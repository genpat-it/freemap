use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use rayon::prelude::*;
use rustc_hash::FxHashMap;

// Prefetch intrinsics for cache optimization
use memmap2::Mmap;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
use std::cell::RefCell;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

// ============================================================================
// Data Structures
// ============================================================================

/// Anchor represents a k-mer match between read and reference.
#[derive(Debug, Clone, Copy, Default)]
struct Anchor {
    read_start: u32,
    ref_start: u32,
    len: u16,
}

#[derive(Debug, Clone)]
struct AlignmentResult {
    ref_start: usize,
    ref_end: usize,
    read_start: usize,
    read_end: usize,
    cigar: String,
    chain_score: i32,
    mapq: u8,
    strand: char,
    _chrom_idx: usize, // Index into chromosome list
}

#[derive(Debug, Clone)]
struct FastqRecord {
    name: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

/// Paired-end read pair
#[derive(Debug, Clone)]
struct ReadPair {
    name: String, // Base name (without /1, /2 suffix)
    r1: FastqRecord,
    r2: FastqRecord,
}


/// A supplementary alignment with clip boundaries in the original read
#[derive(Debug, Clone)]
struct SupplementaryAlignment {
    result: AlignmentResult,
    clip_start: usize, // Start position of the clip in the original forward-strand read
    clip_len: usize,   // Length of the clip portion
}

/// Result containing primary + secondary + supplementary alignments for --multi mode
#[derive(Debug, Clone)]
struct MultiAlignmentResult {
    primary: AlignmentResult,
    secondaries: Vec<AlignmentResult>,          // Alternative placements (SAM flag 0x100)
    supplementaries: Vec<SupplementaryAlignment>, // Split-read portions (SAM flag 0x800)
}

/// Stores chromosome name, start offset in concatenated reference, and length
#[derive(Debug, Clone)]
struct ChromInfo {
    name: String,
    offset: usize, // Start position in concatenated reference
    len: usize,
}

/// Global reference with all chromosomes concatenated
#[derive(Debug)]
pub struct GlobalReference {
    sequence: Vec<u8>,      // Concatenated sequence
    chroms: Vec<ChromInfo>, // Chromosome metadata
}

impl GlobalReference {
    /// Build from FASTA records
    fn from_records(records: Vec<(String, Vec<u8>)>) -> Self {
        let total_len: usize = records.iter().map(|(_, seq)| seq.len()).sum();
        let mut sequence = Vec::with_capacity(total_len);
        let mut chroms = Vec::with_capacity(records.len());

        for (name, seq) in records {
            let offset = sequence.len();
            let len = seq.len();
            chroms.push(ChromInfo { name, offset, len });
            sequence.extend(seq);
        }

        GlobalReference { sequence, chroms }
    }

    /// Convert global position to (chrom_idx, local_position)
    fn global_to_local(&self, global_pos: usize) -> (usize, usize) {
        // Binary search for the chromosome
        let mut lo = 0;
        let mut hi = self.chroms.len();
        while lo < hi {
            let mid = (lo + hi) / 2;
            let chrom = &self.chroms[mid];
            if global_pos < chrom.offset {
                hi = mid;
            } else if global_pos >= chrom.offset + chrom.len {
                lo = mid + 1;
            } else {
                return (mid, global_pos - chrom.offset);
            }
        }
        // Fallback to last chromosome if position is beyond
        let last = self.chroms.len() - 1;
        (last, global_pos.saturating_sub(self.chroms[last].offset))
    }

    /// Get chromosome name by index
    fn chrom_name(&self, idx: usize) -> &str {
        &self.chroms[idx].name
    }

    /// Get chromosome length by index
    fn chrom_len(&self, idx: usize) -> usize {
        self.chroms[idx].len
    }

    /// Check if a chromosome is a primary/canonical contig (not decoy/alt)
    /// Primary: chr1-22, chrX, chrY, chrM (with or without 'chr' prefix)
    /// Also handles NCBI accessions: CM000663.2-CM000686.2, J01415.2 (chrM)
    fn is_primary_contig(&self, idx: usize) -> bool {
        let name = &self.chroms[idx].name;

        // UCSC style: chr1-22, chrX, chrY, chrM
        if name.starts_with("chr") {
            let suffix = &name[3..];
            // Check for numeric chromosomes 1-22
            if let Ok(n) = suffix.parse::<u32>() {
                return n >= 1 && n <= 22;
            }
            // Check for X, Y, M
            return suffix == "X" || suffix == "Y" || suffix == "M";
        }

        // Ensembl/NCBI style without chr prefix: 1-22, X, Y, MT
        if let Ok(n) = name.parse::<u32>() {
            return n >= 1 && n <= 22;
        }
        if name == "X" || name == "Y" || name == "MT" {
            return true;
        }

        // GRCh38 RefSeq accessions: CM000663.2 (chr1) to CM000686.2 (chrY)
        // CM000663.2 = chr1, CM000664.2 = chr2, ..., CM000684.2 = chr22
        // CM000685.2 = chrX, CM000686.2 = chrY
        // J01415.2 = chrM (mitochondria)
        if name.starts_with("CM0006") {
            // CM000663 to CM000686 (covers 63-86 range)
            if let Some(num_str) = name.strip_prefix("CM").and_then(|s| s.split('.').next()) {
                if let Ok(n) = num_str.parse::<u32>() {
                    return n >= 663 && n <= 686;
                }
            }
        }
        if name.starts_with("J01415") {
            return true; // Mitochondria
        }

        // NC_ accessions for primary chromosomes
        if name.starts_with("NC_") {
            // NC_000001 to NC_000024 are primary human chromosomes
            if let Some(num_str) = name.strip_prefix("NC_").and_then(|s| s.split('.').next()) {
                if let Ok(n) = num_str.parse::<u32>() {
                    return n >= 1 && n <= 24;
                }
            }
        }

        // Everything else is considered non-primary (decoy, alt, random, etc.)
        // This includes: chrUn*, *_random, KI*, GL*, HLA-*, chr*_alt, etc.
        false
    }

    /// Check if a global position is on a primary contig
    fn is_primary_position(&self, global_pos: usize) -> bool {
        let (chrom_idx, _) = self.global_to_local(global_pos);
        self.is_primary_contig(chrom_idx)
    }

    /// Check if a contig is blacklisted (known problematic decoys)
    /// GL000209.2 is a mitochondrial decoy that causes 62 misassignments
    fn is_blacklisted_contig(&self, idx: usize) -> bool {
        let name = &self.chroms[idx].name;
        // Mitochondrial decoy - causes confusion with real mitochondria (J01415.2)
        if name.starts_with("GL000209") {
            return true;
        }
        false
    }

    /// Check if a global position is on a blacklisted contig
    #[allow(dead_code)]
    fn is_blacklisted_position(&self, global_pos: usize) -> bool {
        let (chrom_idx, _) = self.global_to_local(global_pos);
        self.is_blacklisted_contig(chrom_idx)
    }

    /// Build blacklist mask for fast position filtering
    /// Returns (start, end) ranges of blacklisted regions
    fn get_blacklisted_ranges(&self) -> Vec<(usize, usize)> {
        let mut ranges = Vec::new();
        for (idx, chrom) in self.chroms.iter().enumerate() {
            if self.is_blacklisted_contig(idx) {
                ranges.push((chrom.offset, chrom.offset + chrom.len));
            }
        }
        ranges
    }
}

// ============================================================================
// FASTA/FASTQ Parsing
// ============================================================================

fn parse_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
    let file = File::open(path).expect("Cannot open FASTA file");
    let reader = BufReader::with_capacity(1 << 20, file);
    let mut records = Vec::new();
    let mut current_name = String::new();
    let mut current_seq = Vec::new();

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if line.starts_with('>') {
            if !current_name.is_empty() {
                records.push((current_name.clone(), current_seq.clone()));
                current_seq.clear();
            }
            current_name = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
        } else {
            current_seq.extend(
                line.trim()
                    .as_bytes()
                    .iter()
                    .map(|&b| b.to_ascii_uppercase()),
            );
        }
    }
    if !current_name.is_empty() {
        records.push((current_name, current_seq));
    }
    records
}

fn parse_fastq(path: &str) -> Vec<FastqRecord> {
    let file = File::open(path).unwrap_or_else(|e| panic!("Cannot open FASTQ file: {}", e));
    // Transparently decompress gzip input. Without this, a .gz file was read as raw bytes,
    // no line started with '@', and parsing silently returned 0 reads (minimap2 reads .gz).
    let mut reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        Box::new(BufReader::with_capacity(
            1 << 20,
            flate2::read::MultiGzDecoder::new(file),
        ))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, file))
    };

    let file_size = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
    let estimated_records = (file_size / 400) as usize;
    let mut records = Vec::with_capacity(estimated_records);

    let mut header_buf = String::with_capacity(256);
    let mut seq_buf = String::with_capacity(256);
    let mut skip_buf = String::with_capacity(256);
    let mut qual_buf = String::with_capacity(256);

    loop {
        header_buf.clear();
        if reader.read_line(&mut header_buf).unwrap_or(0) == 0 {
            break;
        }

        if !header_buf.starts_with('@') {
            continue;
        }

        let header_bytes = header_buf.as_bytes();
        let name_end = header_bytes[1..]
            .iter()
            .position(|&b| b == b' ' || b == b'\t' || b == b'\n' || b == b'\r')
            .map(|p| p + 1)
            .unwrap_or(header_buf.len() - 1);
        let name = header_buf[1..name_end].to_string();

        seq_buf.clear();
        if reader.read_line(&mut seq_buf).unwrap_or(0) == 0 {
            break;
        }
        let seq_len = seq_buf.trim_end().len();
        let seq: Vec<u8> = seq_buf.as_bytes()[..seq_len].to_vec();

        skip_buf.clear();
        let _ = reader.read_line(&mut skip_buf);

        qual_buf.clear();
        let _ = reader.read_line(&mut qual_buf);
        let qual: Vec<u8> = qual_buf.trim_end().as_bytes().to_vec();

        records.push(FastqRecord { name, seq, qual });
    }
    records
}

/// Strip paired-end suffix from read name (/1, /2, .1, .2, _R1, _R2, etc.)
fn strip_pair_suffix(name: &str) -> String {
    // Common suffixes to strip
    if name.ends_with("/1") || name.ends_with("/2") || name.ends_with(".1") || name.ends_with(".2")
    {
        name[..name.len() - 2].to_string()
    } else if name.ends_with("_R1") || name.ends_with("_R2") {
        name[..name.len() - 3].to_string()
    } else {
        name.to_string()
    }
}

/// Load paired-end reads from two FASTQ files and pair them by name
fn load_paired_reads(r1_path: &str, r2_path: &str) -> Vec<ReadPair> {
    let r1_reads = parse_fastq(r1_path);
    let r2_reads = parse_fastq(r2_path);

    // Build a map of R1 reads by normalized name
    let mut r1_map: FxHashMap<String, FastqRecord> = FxHashMap::default();
    for read in r1_reads {
        let base_name = strip_pair_suffix(&read.name);
        r1_map.insert(base_name, read);
    }

    // Match R2 reads to R1 reads
    let mut pairs = Vec::with_capacity(r2_reads.len());
    for r2 in r2_reads {
        let base_name = strip_pair_suffix(&r2.name);
        if let Some(r1) = r1_map.remove(&base_name) {
            pairs.push(ReadPair {
                name: base_name,
                r1,
                r2,
            });
        }
    }

    pairs
}

// ============================================================================
// Reverse Complement
// ============================================================================

#[inline(always)]
fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

// ============================================================================
// HOMOPOLYMER COMPRESSION (DP-free improvement for ONT)
// ============================================================================
//
// ONT errors are dominated by indels in homopolymers (e.g., AAAAA → AAAA or AAAAAA).
// By compressing runs before alignment, we eliminate this systematic noise source.
//
// Example: ACGTTTTTACGAAAA → ACGTACGA (compressed)
//          Position map: [0,1,2,3,8,9,10,11,15] (original positions)
//
// This is strictly DP-free: just O(n) compression and position mapping.

/// Compress homopolymers in a sequence, returning compressed sequence and position map.
/// Position map[i] = original position of compressed base i
fn compress_homopolymers(seq: &[u8]) -> (Vec<u8>, Vec<u32>) {
    if seq.is_empty() {
        return (Vec::new(), Vec::new());
    }

    let mut compressed = Vec::with_capacity(seq.len());
    let mut pos_map = Vec::with_capacity(seq.len());

    let mut prev_base = seq[0];
    compressed.push(prev_base);
    pos_map.push(0u32);

    for (i, &base) in seq.iter().enumerate().skip(1) {
        if base != prev_base {
            compressed.push(base);
            pos_map.push(i as u32);
            prev_base = base;
        }
    }

    (compressed, pos_map)
}

/// Convert position in compressed sequence back to original sequence position
#[inline(always)]
fn decompress_pos(pos: usize, pos_map: &[u32]) -> usize {
    if pos < pos_map.len() {
        pos_map[pos] as usize
    } else {
        // Beyond compressed sequence - extrapolate
        pos_map
            .last()
            .map(|&p| p as usize + (pos - pos_map.len() + 1))
            .unwrap_or(pos)
    }
}

// ============================================================================
// Index Serialization (pre-indexing support)
// ============================================================================

const INDEX_MAGIC: &[u8; 4] = b"GEOI"; // GEO Index magic bytes
const INDEX_VERSION: u8 = 6; // Version 6: embedded reference (v5 + chrom metadata + ref sequence)

// Mmap index format (version 6):
// [Header: 64 bytes]
//   - magic: 4 bytes "GEOI"
//   - version: 1 byte (6)
//   - flags: 1 byte (bit 0 = homopolymer mode)
//   - k: 2 bytes (u16)
//   - w: 2 bytes (u16)
//   - max_freq: 4 bytes (u32)
//   - n_buckets: 4 bytes (u32) - power of 2 for fast modulo
//   - n_positions: 8 bytes (u64) - total positions stored
//   - ref_checksum: 8 bytes (u64)
//   - ref_len: 8 bytes (u64)
//   - compressed_ref_len: 8 bytes (u64) - 0 if not homopolymer
//   - pos_map_len: 8 bytes (u64) - 0 if not homopolymer
//   - n_chroms: 2 bytes (u16) - number of chromosomes (0 = no embedded ref)
//   - chrom_meta_len: 4 bytes (u32) - total bytes of chrom metadata section
// [Bucket table: n_buckets * 16 bytes]
//   - Each bucket: offset (8 bytes) + count (4 bytes) + padding (4 bytes)
// [Positions: variable]
//   - For each entry: hash (8 bytes) + positions (4 bytes each)
// [Compressed ref: compressed_ref_len bytes] (if homopolymer)
// [Pos map: pos_map_len * 4 bytes] (if homopolymer)
// [Chrom metadata: chrom_meta_len bytes] (if n_chroms > 0)
//   - For each chrom: name_len(u16) + name(bytes) + offset(u64) + len(u64)
// [Reference sequence: ref_len bytes] (if n_chroms > 0)

const HEADER_SIZE: usize = 64;
const BUCKET_SIZE: usize = 16; // offset (8) + count (4) + padding (4)

#[repr(C, packed)]
#[derive(Clone, Copy)]
#[allow(dead_code)]
struct MmapHeader {
    magic: [u8; 4],
    version: u8,
    flags: u8,
    k: u16,
    w: u16,
    max_freq: u32,
    n_buckets: u32,
    n_positions: u64,
    ref_checksum: u64,
    ref_len: u64,
    compressed_ref_len: u64,
    pos_map_len: u64,
    n_chroms: u16,        // Number of chromosomes (0 = no embedded ref)
    chrom_meta_len: u32,  // Total bytes of chromosome metadata section
}

#[repr(C, packed)]
#[derive(Clone, Copy, Default)]
#[allow(dead_code)]
struct BucketEntry {
    offset: u64, // Offset into positions section
    count: u32,  // Number of entries in this bucket
    _padding: u32,
}

/// Trait for k-mer index lookup (works with both FxHashMap and MmapIndex)
pub trait KmerIndex: Sync {
    fn get_positions(&self, hash: u64) -> Option<&[u32]>;
    /// Prefetch data for a hash (hint to CPU cache)
    fn prefetch(&self, _hash: u64) {}
}

impl KmerIndex for FxHashMap<u64, Vec<u32>> {
    #[inline]
    fn get_positions(&self, hash: u64) -> Option<&[u32]> {
        self.get(&hash).map(|v| v.as_slice())
    }

    #[inline]
    fn prefetch(&self, _hash: u64) {
        // FxHashMap prefetch is a no-op - we can't access internal buckets
        // The sorted lookup order provides cache locality instead
    }
}

/// Sorted array-based index with binary search lookup
/// More memory-efficient and faster to build than HashMap
#[derive(Debug)]
pub struct SortedIndex {
    // Sorted by hash, stores (hash, start_idx, count) for each unique hash
    // positions[start_idx..start_idx+count] are the reference positions
    hash_table: Vec<(u64, u32, u16)>, // (hash, start_index, count)
    positions: Vec<u32>,
}

impl SortedIndex {
    /// Check if the index is empty
    pub fn is_empty(&self) -> bool {
        self.hash_table.is_empty()
    }

    /// Get number of unique k-mers
    pub fn len(&self) -> usize {
        self.hash_table.len()
    }

    /// Iterate over (hash, positions) pairs
    pub fn iter(&self) -> impl Iterator<Item = (u64, &[u32])> + '_ {
        self.hash_table.iter().map(|&(hash, start, count)| {
            let start_idx = start as usize;
            let end_idx = start_idx + count as usize;
            (hash, &self.positions[start_idx..end_idx])
        })
    }

    /// Convert to FxHashMap for O(1) lookup (useful when lookup dominates)
    pub fn to_hashmap(&self) -> FxHashMap<u64, Vec<u32>> {
        let mut map =
            FxHashMap::with_capacity_and_hasher(self.hash_table.len(), Default::default());
        for &(hash, start, count) in &self.hash_table {
            let start_idx = start as usize;
            let end_idx = start_idx + count as usize;
            map.insert(hash, self.positions[start_idx..end_idx].to_vec());
        }
        map
    }
}

impl KmerIndex for SortedIndex {
    #[inline]
    fn get_positions(&self, hash: u64) -> Option<&[u32]> {
        // Binary search for the full 64-bit hash
        let idx = self
            .hash_table
            .binary_search_by_key(&hash, |&(h, _, _)| h)
            .ok()?;
        let (_, start, count) = self.hash_table[idx];
        let start = start as usize;
        let end = start + count as usize;
        Some(&self.positions[start..end])
    }
}

/// Mmap-backed index for near-instant loading
#[derive(Debug)]
pub struct MmapIndex {
    mmap: Mmap,
    n_buckets: u32,
    bucket_mask: u64,
    positions_offset: usize,
    compressed_ref_offset: usize,
    pos_map_offset: usize,
    compressed_ref_len: usize,
    pos_map_len: usize,
    // Embedded reference (v6)
    n_chroms: usize,
    chrom_meta_offset: usize,
    ref_seq_offset: usize,
    ref_seq_len: usize,
}

impl KmerIndex for MmapIndex {
    /// Look up positions for a k-mer hash
    #[inline]
    fn get_positions(&self, hash: u64) -> Option<&[u32]> {
        let bucket_idx = (hash & self.bucket_mask) as usize;
        let bucket_offset = HEADER_SIZE + bucket_idx * BUCKET_SIZE;

        #[cfg(target_arch = "x86_64")]
        unsafe {
            _mm_prefetch(
                self.mmap.as_ptr().add(bucket_offset) as *const i8,
                _MM_HINT_T0,
            );
        }

        // Read bucket entry
        let bucket_data = &self.mmap[bucket_offset..bucket_offset + BUCKET_SIZE];
        let entry_offset = u64::from_le_bytes(bucket_data[0..8].try_into().unwrap()) as usize;
        let entry_count = u32::from_le_bytes(bucket_data[8..12].try_into().unwrap()) as usize;

        if entry_count == 0 {
            return None;
        }

        // Linear search within bucket (buckets are small due to good hash distribution)
        // Compact format: hash64(8) + n_pos16(2) + padding(2) + positions(4*n)
        let mut pos = self.positions_offset + entry_offset;
        for _ in 0..entry_count {
            #[cfg(target_arch = "x86_64")]
            unsafe {
                _mm_prefetch(self.mmap.as_ptr().add(pos) as *const i8, _MM_HINT_T0);
            }
            let stored_hash64 = u64::from_le_bytes(self.mmap[pos..pos + 8].try_into().unwrap());
            pos += 8;
            let n_pos = u16::from_le_bytes(self.mmap[pos..pos + 2].try_into().unwrap()) as usize;
            pos += 2;
            pos += 2; // Skip padding

            if stored_hash64 == hash {
                // Found it! Return slice of positions (now aligned thanks to padding)
                let positions_ptr = &self.mmap[pos..pos + n_pos * 4];
                let positions = unsafe {
                    std::slice::from_raw_parts(positions_ptr.as_ptr() as *const u32, n_pos)
                };
                return Some(positions);
            }
            pos += n_pos * 4;
        }
        None
    }

    #[inline]
    fn prefetch(&self, hash: u64) {
        let bucket_idx = (hash & self.bucket_mask) as usize;
        let bucket_offset = HEADER_SIZE + bucket_idx * BUCKET_SIZE;
        #[cfg(target_arch = "x86_64")]
        unsafe {
            _mm_prefetch(
                self.mmap.as_ptr().add(bucket_offset) as *const i8,
                _MM_HINT_T0,
            );
        }
    }
}

impl MmapIndex {
    /// Get compressed reference (for homopolymer mode)
    pub fn compressed_ref(&self) -> &[u8] {
        if self.compressed_ref_len == 0 {
            return &[];
        }
        &self.mmap[self.compressed_ref_offset..self.compressed_ref_offset + self.compressed_ref_len]
    }

    /// Get position map (for homopolymer mode)
    pub fn pos_map(&self) -> &[u32] {
        if self.pos_map_len == 0 {
            return &[];
        }
        let data = &self.mmap[self.pos_map_offset..self.pos_map_offset + self.pos_map_len * 4];
        unsafe { std::slice::from_raw_parts(data.as_ptr() as *const u32, self.pos_map_len) }
    }

    /// Check if this index has an embedded reference
    pub fn has_embedded_ref(&self) -> bool {
        self.n_chroms > 0
    }

    /// Extract GlobalReference from embedded data (fast memcpy, no FASTA parsing)
    pub fn extract_reference(&self) -> Option<GlobalReference> {
        if self.n_chroms == 0 {
            return None;
        }

        // Parse chromosome metadata
        let mut chroms = Vec::with_capacity(self.n_chroms);
        let mut pos = self.chrom_meta_offset;

        for _ in 0..self.n_chroms {
            let name_len = u16::from_le_bytes(
                self.mmap[pos..pos + 2].try_into().unwrap(),
            ) as usize;
            pos += 2;
            let name = std::str::from_utf8(&self.mmap[pos..pos + name_len])
                .unwrap_or("unknown")
                .to_string();
            pos += name_len;
            let offset = u64::from_le_bytes(
                self.mmap[pos..pos + 8].try_into().unwrap(),
            ) as usize;
            pos += 8;
            let len = u64::from_le_bytes(
                self.mmap[pos..pos + 8].try_into().unwrap(),
            ) as usize;
            pos += 8;
            chroms.push(ChromInfo { name, offset, len });
        }

        // Copy reference sequence from mmap (fast memcpy ~1s for 3.1GB)
        let sequence = self.mmap[self.ref_seq_offset..self.ref_seq_offset + self.ref_seq_len].to_vec();

        Some(GlobalReference { sequence, chroms })
    }

    /// Get number of entries from header
    pub fn len(&self) -> usize {
        // n_entries is stored at position 18-26 in header (n_positions field actually stores entry count)
        // Actually we need to count unique hashes - use n_buckets as approximation
        // For now just return a reasonable estimate
        self.n_buckets as usize // This is just for display, not critical
    }
}

/// Enum to hold either FxHashMap or MmapIndex
pub enum IndexType {
    HashMap(FxHashMap<u64, Vec<u32>>),
    Sorted(SortedIndex),
    Mmap(MmapIndex),
}

impl KmerIndex for IndexType {
    #[inline]
    fn get_positions(&self, hash: u64) -> Option<&[u32]> {
        match self {
            IndexType::HashMap(m) => m.get_positions(hash),
            IndexType::Sorted(s) => s.get_positions(hash),
            IndexType::Mmap(m) => m.get_positions(hash),
        }
    }
}

impl IndexType {
    pub fn len(&self) -> usize {
        match self {
            IndexType::HashMap(m) => m.len(),
            IndexType::Sorted(s) => s.len(),
            IndexType::Mmap(m) => m.len(),
        }
    }
}

/// Save index to file in mmap-friendly format (version 6: with embedded reference)
fn save_index(
    path: &str,
    index: &SortedIndex,
    reference: &[u8],
    compressed_ref: &[u8],
    pos_map: &[u32],
    k: usize,
    w: usize,
    max_freq: usize,
    homopolymer_mode: bool,
    global_ref: &GlobalReference,
) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    // Calculate number of buckets (power of 2, n/8 = best tradeoff)
    let n_entries = index.hash_table.len();
    let n_buckets = (n_entries / 8).next_power_of_two().max(1024) as u32;
    let bucket_mask = (n_buckets - 1) as u64;

    // Group entries by bucket
    let mut buckets: Vec<Vec<(u64, u32, u16)>> = vec![Vec::new(); n_buckets as usize];
    for &(hash, start, count) in &index.hash_table {
        let bucket_idx = (hash & bucket_mask) as usize;
        buckets[bucket_idx].push((hash, start, count));
    }

    // Calculate total positions for header
    let n_positions: u64 = index.positions.len() as u64;

    // Compute checksum
    let checksum: u64 = reference
        .iter()
        .fold(0u64, |acc, &b| acc.wrapping_mul(31).wrapping_add(b as u64));

    // Write header (64 bytes)
    writer.write_all(INDEX_MAGIC)?; // 4 bytes
    writer.write_all(&[INDEX_VERSION])?; // 1 byte
    writer.write_all(&[if homopolymer_mode { 1u8 } else { 0u8 }])?; // 1 byte flags
    writer.write_all(&(k as u16).to_le_bytes())?; // 2 bytes
    writer.write_all(&(w as u16).to_le_bytes())?; // 2 bytes
    writer.write_all(&(max_freq as u32).to_le_bytes())?; // 4 bytes
    writer.write_all(&n_buckets.to_le_bytes())?; // 4 bytes
    writer.write_all(&n_positions.to_le_bytes())?; // 8 bytes
    writer.write_all(&checksum.to_le_bytes())?; // 8 bytes
    writer.write_all(&(reference.len() as u64).to_le_bytes())?; // 8 bytes
    writer.write_all(&(compressed_ref.len() as u64).to_le_bytes())?; // 8 bytes
    writer.write_all(&(pos_map.len() as u64).to_le_bytes())?; // 8 bytes
    // Compute chrom metadata size
    let chrom_meta_len: u32 = global_ref.chroms.iter()
        .map(|c| 2 + c.name.len() + 8 + 8) // name_len(u16) + name + offset(u64) + len(u64)
        .sum::<usize>() as u32;
    writer.write_all(&(global_ref.chroms.len() as u16).to_le_bytes())?; // 2 bytes n_chroms
    writer.write_all(&chrom_meta_len.to_le_bytes())?; // 4 bytes chrom_meta_len = 64 total

    // Calculate offsets for each bucket and write bucket table
    let mut current_offset: u64 = 0;
    for bucket in &buckets {
        writer.write_all(&current_offset.to_le_bytes())?; // 8 bytes offset
        writer.write_all(&(bucket.len() as u32).to_le_bytes())?; // 4 bytes count
        writer.write_all(&[0u8; 4])?; // 4 bytes padding

        // Calculate size of this bucket's entries
        for &(_, start, count) in bucket {
            current_offset += 8 + 2 + 2 + (count as u64) * 4; // hash64 + n_pos16 + pad + positions
            let _ = start; // used to access positions
        }
    }

    // Write positions data with compact format (64-bit hash)
    for bucket in &buckets {
        for &(hash, start, count) in bucket {
            writer.write_all(&hash.to_le_bytes())?; // 8 bytes (full 64-bit hash)
            writer.write_all(&count.to_le_bytes())?; // 2 bytes
            writer.write_all(&[0u8; 2])?; // 2 bytes padding for alignment
            let start_idx = start as usize;
            let end_idx = start_idx + count as usize;
            for &pos in &index.positions[start_idx..end_idx] {
                writer.write_all(&pos.to_le_bytes())?;
            }
        }
    }

    // Write compressed reference (if homopolymer mode)
    if homopolymer_mode {
        writer.write_all(compressed_ref)?;
        // Pad to 4-byte alignment so pos_map u32 values are aligned for mmap access
        let padding = (4 - (compressed_ref.len() % 4)) % 4;
        if padding > 0 {
            writer.write_all(&vec![0u8; padding])?;
        }
        for &p in pos_map {
            writer.write_all(&p.to_le_bytes())?;
        }
    }

    // Write embedded reference (v6): chromosome metadata + raw sequence
    if !global_ref.chroms.is_empty() {
        // Chromosome metadata: for each chrom: name_len(u16) + name(bytes) + offset(u64) + len(u64)
        for chrom in &global_ref.chroms {
            writer.write_all(&(chrom.name.len() as u16).to_le_bytes())?;
            writer.write_all(chrom.name.as_bytes())?;
            writer.write_all(&(chrom.offset as u64).to_le_bytes())?;
            writer.write_all(&(chrom.len as u64).to_le_bytes())?;
        }
        // Raw reference sequence
        writer.write_all(&global_ref.sequence)?;
    }

    writer.flush()?;
    Ok(())
}

/// Load index from file using mmap for near-instant loading
/// If reference is None, skips checksum verification (only valid for v6 indices with embedded ref).
fn load_index_mmap(
    path: &str,
    reference: Option<&[u8]>,
) -> Result<(MmapIndex, usize, usize, usize, bool), String> {
    let file = File::open(path).map_err(|e| format!("Cannot open index: {}", e))?;

    // Memory-map the file
    let mmap = unsafe { Mmap::map(&file) }.map_err(|e| format!("Cannot mmap index: {}", e))?;

    if mmap.len() < HEADER_SIZE {
        return Err("Index file too small".to_string());
    }

    // Read and verify header
    if &mmap[0..4] != INDEX_MAGIC {
        return Err("Invalid index file (wrong magic)".to_string());
    }

    let version = mmap[4];
    if version != 5 && version != 6 {
        return Err(format!(
            "Unsupported index version: {} (expected 5 or 6)",
            version
        ));
    }

    // Parse header fields
    let flags = mmap[5];
    let homopolymer_mode = (flags & 1) != 0;
    let k = u16::from_le_bytes(mmap[6..8].try_into().unwrap()) as usize;
    let w = u16::from_le_bytes(mmap[8..10].try_into().unwrap()) as usize;
    let max_freq = u32::from_le_bytes(mmap[10..14].try_into().unwrap()) as usize;
    let n_buckets = u32::from_le_bytes(mmap[14..18].try_into().unwrap());
    // n_positions at bytes 18-26 (not needed for loading)
    let stored_checksum = u64::from_le_bytes(mmap[26..34].try_into().unwrap());
    let stored_ref_len = u64::from_le_bytes(mmap[34..42].try_into().unwrap()) as usize;
    let compressed_ref_len = u64::from_le_bytes(mmap[42..50].try_into().unwrap()) as usize;
    let pos_map_len = u64::from_le_bytes(mmap[50..58].try_into().unwrap()) as usize;

    // v6 fields (bytes 58-63); for v5 these were zero padding
    let n_chroms = u16::from_le_bytes(mmap[58..60].try_into().unwrap()) as usize;
    let chrom_meta_len = u32::from_le_bytes(mmap[60..64].try_into().unwrap()) as usize;

    // Verify reference (skip if None - caller will use embedded ref)
    if let Some(reference) = reference {
        let computed_checksum: u64 = reference
            .iter()
            .fold(0u64, |acc, &b| acc.wrapping_mul(31).wrapping_add(b as u64));
        if stored_checksum != computed_checksum {
            return Err("Reference mismatch: index was built for a different reference".to_string());
        }
        if stored_ref_len != reference.len() {
            return Err(format!(
                "Reference length mismatch: {} vs {}",
                stored_ref_len,
                reference.len()
            ));
        }
    } else if n_chroms == 0 {
        return Err("Index has no embedded reference and no external reference provided".to_string());
    }

    // Calculate section offsets
    let bucket_table_size = (n_buckets as usize) * BUCKET_SIZE;
    let positions_offset = HEADER_SIZE + bucket_table_size;

    // Calculate trailing data size (everything after positions section)
    let compressed_ref_padding = if compressed_ref_len > 0 {
        (4 - (compressed_ref_len % 4)) % 4
    } else {
        0
    };
    let hpc_data = compressed_ref_len + compressed_ref_padding + pos_map_len * 4;
    let embedded_ref_data = if n_chroms > 0 {
        chrom_meta_len + stored_ref_len
    } else {
        0
    };
    let trailing_data = hpc_data + embedded_ref_data;

    let file_size = mmap.len();
    let positions_end = file_size - trailing_data;

    let compressed_ref_offset = positions_end;
    let pos_map_offset = compressed_ref_offset + compressed_ref_len + compressed_ref_padding;
    let pos_map_end = pos_map_offset + pos_map_len * 4;

    // Embedded reference offsets
    let chrom_meta_offset = pos_map_end;
    let ref_seq_offset = chrom_meta_offset + chrom_meta_len;
    let ref_seq_len = if n_chroms > 0 { stored_ref_len } else { 0 };

    let bucket_mask = (n_buckets - 1) as u64;

    let index = MmapIndex {
        mmap,
        n_buckets,
        bucket_mask,
        positions_offset,
        compressed_ref_offset,
        pos_map_offset,
        compressed_ref_len,
        pos_map_len,
        n_chroms,
        chrom_meta_offset,
        ref_seq_offset,
        ref_seq_len,
    };

    Ok((index, k, w, max_freq, homopolymer_mode))
}

// ============================================================================
// Indexing with minimizers for speed
// ============================================================================

// Lookup table for base to 2-bit encoding (0-3 for ACGT, 4 for N)
static BASE_TO_BITS: [u8; 256] = {
    let mut table = [4u8; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

#[inline(always)]
fn base_to_bits(b: u8) -> u64 {
    BASE_TO_BITS[b as usize] as u64
}

/// splitmix64 hash function - excellent avalanche properties
/// Used to mix k-mer encoding for uniform distribution in hash tables
#[inline(always)]
fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}

#[inline(always)]
fn kmer_hash(kmer: &[u8]) -> u64 {
    let mut code: u64 = 0;
    for &b in kmer {
        code = (code << 2) | base_to_bits(b);
    }
    splitmix64(code)
}

/// Compute reverse complement of a k-mer encoding using bit-parallel operations
/// O(1) instead of O(k) - uses parallel bit swaps
/// Encoding: 2 bits per base, LSB = last base, MSB = first base
/// Complement: A(00)↔T(11), C(01)↔G(10) → XOR with 0b11 per base
#[inline(always)]
fn revcomp_encoding(enc: u64, k: usize) -> u64 {
    // Step 1: Complement all bases (XOR each 2-bit pair with 0b11)
    let mask = (1u64 << (2 * k)) - 1;
    let comp = enc ^ mask; // All bases complemented

    // Step 2: Reverse the order of 2-bit pairs using parallel swaps
    // This is a divide-and-conquer bit reversal for 2-bit units
    let mut x = comp;

    // Swap adjacent 2-bit pairs (positions 0-1 with 2-3, 4-5 with 6-7, etc.)
    x = ((x >> 2) & 0x3333_3333_3333_3333) | ((x & 0x3333_3333_3333_3333) << 2);
    // Swap adjacent 4-bit groups
    x = ((x >> 4) & 0x0f0f_0f0f_0f0f_0f0f) | ((x & 0x0f0f_0f0f_0f0f_0f0f) << 4);
    // Swap adjacent bytes
    x = ((x >> 8) & 0x00ff_00ff_00ff_00ff) | ((x & 0x00ff_00ff_00ff_00ff) << 8);
    // Swap adjacent 16-bit groups
    x = ((x >> 16) & 0x0000_ffff_0000_ffff) | ((x & 0x0000_ffff_0000_ffff) << 16);
    // Swap 32-bit halves
    x = (x >> 32) | (x << 32);

    // Step 3: Shift right to align result (we reversed all 32 positions, but only k are valid)
    x >> (64 - 2 * k)
}

/// Get canonical k-mer encoding (lexicographically smaller of forward/revcomp)
/// This ensures forward and reverse complement k-mers map to the same key
#[inline(always)]
fn canonical_encoding(enc: u64, k: usize) -> u64 {
    let rev_enc = revcomp_encoding(enc, k);
    if enc <= rev_enc {
        enc
    } else {
        rev_enc
    }
}

/// Compute canonical hash for a k-mer encoding
/// 1. Canonicalize encoding (min of fwd/revcomp)
/// 2. Apply splitmix64 hash (only ONE hash, not two!)
#[inline(always)]
fn canonical_hash(enc: u64, k: usize) -> u64 {
    let canon = canonical_encoding(enc, k);
    splitmix64(canon)
}

/// Scan minimizers with callback - core logic used by 2-pass index building
/// Uses canonical k-mers with splitmix64 hash (like minimap2)
#[inline]
fn scan_minimizers<F: FnMut(u64, u32)>(reference: &[u8], k: usize, w: usize, mut emit: F) {
    if reference.len() < k {
        return;
    }

    let mask: u64 = (1u64 << (2 * k)) - 1;
    let mut deque_buf: Vec<(u64, u32)> = Vec::with_capacity(w + 1);
    let mut deque_front: usize = 0;
    let mut last_min_pos: u32 = u32::MAX;

    let mut kmer_enc: u64 = 0; // Raw k-mer encoding
    let mut valid_bases: usize = 0;

    // Initialize first k-mer
    for i in 0..k.min(reference.len()) {
        let bits = BASE_TO_BITS[reference[i] as usize];
        if bits > 3 {
            valid_bases = 0;
            kmer_enc = 0;
        } else {
            kmer_enc = ((kmer_enc << 2) | bits as u64) & mask;
            valid_bases += 1;
        }
    }

    if valid_bases == k {
        let hash = canonical_hash(kmer_enc, k); // Canonical hash for minimizer selection
        deque_buf.push((hash, 0));
    }

    for i in 1..=(reference.len() - k) {
        let old_base = BASE_TO_BITS[reference[i - 1] as usize];
        let new_base = BASE_TO_BITS[reference[i + k - 1] as usize];

        if new_base > 3 {
            valid_bases = 0;
            kmer_enc = 0;
            deque_buf.clear();
            deque_front = 0;
            last_min_pos = u32::MAX;
            continue;
        }

        if old_base > 3 {
            kmer_enc = 0;
            valid_bases = 0;
            for j in i..(i + k) {
                let bits = BASE_TO_BITS[reference[j] as usize];
                if bits > 3 {
                    valid_bases = 0;
                    kmer_enc = 0;
                    break;
                }
                kmer_enc = ((kmer_enc << 2) | bits as u64) & mask;
                valid_bases += 1;
            }
            if valid_bases < k {
                deque_buf.clear();
                deque_front = 0;
                last_min_pos = u32::MAX;
                continue;
            }
        } else {
            kmer_enc = ((kmer_enc << 2) | new_base as u64) & mask;
            if valid_bases < k {
                valid_bases += 1;
            }
        }

        if valid_bases < k {
            continue;
        }

        let i32 = i as u32;
        let hash = canonical_hash(kmer_enc, k); // Canonical hash for minimizer selection

        // Monotonic deque: remove elements >= current hash from back
        while deque_buf.len() > deque_front {
            if deque_buf[deque_buf.len() - 1].0 >= hash {
                deque_buf.pop();
            } else {
                break;
            }
        }
        deque_buf.push((hash, i32));

        // Remove elements outside window from front
        let window_start = if i >= w { (i - w + 1) as u32 } else { 0 };
        while deque_front < deque_buf.len() && deque_buf[deque_front].1 < window_start {
            deque_front += 1;
        }

        // Compact buffer periodically
        if deque_front > w * 4 {
            let remaining = deque_buf.len() - deque_front;
            for idx in 0..remaining {
                deque_buf[idx] = deque_buf[deque_front + idx];
            }
            deque_buf.truncate(remaining);
            deque_front = 0;
        }

        // Emit minimizer (hash value, not raw encoding)
        if i >= w - 1 && deque_front < deque_buf.len() {
            let (min_hash, min_pos) = deque_buf[deque_front];
            if min_pos != last_min_pos {
                last_min_pos = min_pos;
                emit(min_hash, min_pos);
            }
        }
    }
}

/// Build index: collect minimizers, sort by hash, build flat index
/// Uses full 64-bit k-mer encoding to avoid collisions for any k value
/// OPTIMIZED: Parallel collection + radix sort for O(n) sorting
fn build_kmer_index(reference: &[u8], k: usize, w: usize, max_hits: usize) -> SortedIndex {
    if reference.len() < k {
        return SortedIndex {
            hash_table: Vec::new(),
            positions: Vec::new(),
        };
    }

    let n_threads = rayon::current_num_threads();

    // For small references or single thread, use sequential approach
    if reference.len() < 1_000_000 || n_threads == 1 {
        return build_kmer_index_sequential(reference, k, w, max_hits);
    }

    // =========================================================================
    // PARALLEL COLLECTION: Split reference into chunks, collect in parallel
    // =========================================================================
    let n_chunks = n_threads;
    let chunk_size = reference.len() / n_chunks;
    let overlap = k + w; // Overlap to avoid missing minimizers at boundaries

    let chunk_results: Vec<Vec<(u64, u32)>> = (0..n_chunks)
        .into_par_iter()
        .map(|chunk_id| {
            let start = chunk_id * chunk_size;
            let end = if chunk_id == n_chunks - 1 {
                reference.len()
            } else {
                ((chunk_id + 1) * chunk_size + overlap).min(reference.len())
            };

            if start >= reference.len() {
                return Vec::new();
            }

            let chunk = &reference[start..end];
            let estimated = chunk.len() / w + 100;
            let mut entries: Vec<(u64, u32)> = Vec::with_capacity(estimated);

            scan_minimizers(chunk, k, w, |hash, local_pos| {
                let global_pos = start as u32 + local_pos;
                // Only emit if position belongs to this chunk (avoid duplicates from overlap)
                let chunk_boundary = if chunk_id == n_chunks - 1 {
                    reference.len() as u32
                } else {
                    ((chunk_id + 1) * chunk_size) as u32
                };
                if global_pos < chunk_boundary {
                    entries.push((hash, global_pos));
                }
            });

            entries
        })
        .collect();

    // Concatenate all chunks
    let total_count: usize = chunk_results.iter().map(|v| v.len()).sum();
    let mut entries: Vec<(u64, u32)> = Vec::with_capacity(total_count);
    for chunk in chunk_results {
        entries.extend(chunk);
    }

    if entries.is_empty() {
        return SortedIndex {
            hash_table: Vec::new(),
            positions: Vec::new(),
        };
    }

    // =========================================================================
    // PARALLEL SORT: Uses rayon's work-stealing parallel sort
    // =========================================================================
    entries.par_sort_unstable_by_key(|&(h, _)| h);

    // Build SortedIndex from sorted entries
    build_sorted_index_from_entries(&entries, max_hits)
}

/// Sequential index building for small references or single-threaded mode
fn build_kmer_index_sequential(
    reference: &[u8],
    k: usize,
    w: usize,
    max_hits: usize,
) -> SortedIndex {
    let estimated = reference.len() / w + 1000;
    let mut entries: Vec<(u64, u32)> = Vec::with_capacity(estimated);

    scan_minimizers(reference, k, w, |h, pos| {
        entries.push((h, pos));
    });

    if entries.is_empty() {
        return SortedIndex {
            hash_table: Vec::new(),
            positions: Vec::new(),
        };
    }

    // Standard sort for sequential case
    entries.sort_unstable_by_key(|&(h, _)| h);

    build_sorted_index_from_entries(&entries, max_hits)
}

/// Build SortedIndex from sorted (hash64, pos32) entries
fn build_sorted_index_from_entries(entries: &[(u64, u32)], max_hits: usize) -> SortedIndex {
    if entries.is_empty() {
        return SortedIndex {
            hash_table: Vec::new(),
            positions: Vec::new(),
        };
    }

    // Count unique hashes (full 64-bit)
    let mut unique_count = 1usize;
    let mut prev_hash = entries[0].0;
    for &(h, _) in &entries[1..] {
        if h != prev_hash {
            unique_count += 1;
            prev_hash = h;
        }
    }

    let mut hash_table: Vec<(u64, u32, u16)> = Vec::with_capacity(unique_count);
    let mut positions: Vec<u32> = Vec::with_capacity(entries.len());

    let mut i = 0;
    while i < entries.len() {
        let current_hash = entries[i].0;
        let start_idx = positions.len() as u32;
        let group_start = i;

        // Find end of this hash group and collect positions
        while i < entries.len() && entries[i].0 == current_hash {
            i += 1;
        }

        let count = i - group_start;
        if count < max_hits && count <= u16::MAX as usize {
            // Add all positions for this hash
            for j in group_start..i {
                positions.push(entries[j].1);
            }
            // Store full 64-bit hash
            hash_table.push((current_hash, start_idx, count as u16));
        }
    }

    SortedIndex {
        hash_table,
        positions,
    }
}

/// Count 32-bit hash collisions in minimizer index
/// A collision occurs when two DIFFERENT 64-bit hashes map to the same 32-bit hash
/// Returns (total_minimizers, unique_hash64, unique_hash32, colliding_buckets, collision_positions)
fn count_hash_collisions(
    reference: &[u8],
    k: usize,
    w: usize,
) -> (usize, usize, usize, usize, usize) {
    use std::collections::HashMap;

    // Map: hash32 -> set of distinct hash64 values
    let mut hash32_to_hash64s: HashMap<u32, Vec<u64>> = HashMap::new();
    let mut total_minimizers = 0usize;

    scan_minimizers(reference, k, w, |h64, _pos| {
        total_minimizers += 1;
        let h32 = h64 as u32;
        hash32_to_hash64s.entry(h32).or_default().push(h64);
    });

    // Count unique hash64 and hash32
    let unique_hash32 = hash32_to_hash64s.len();

    // For each hash32 bucket, count distinct hash64 values
    let mut unique_hash64 = 0usize;
    let mut colliding_buckets = 0usize;
    let mut collision_positions = 0usize;

    for (_h32, hash64_list) in &hash32_to_hash64s {
        // Deduplicate hash64 values in this bucket
        let mut unique_in_bucket: Vec<u64> = hash64_list.clone();
        unique_in_bucket.sort_unstable();
        unique_in_bucket.dedup();

        unique_hash64 += unique_in_bucket.len();

        // If more than one distinct hash64 maps to this hash32, it's a collision
        if unique_in_bucket.len() > 1 {
            colliding_buckets += 1;
            // Count how many positions are affected by collisions
            collision_positions += hash64_list.len();
        }
    }

    (
        total_minimizers,
        unique_hash64,
        unique_hash32,
        colliding_buckets,
        collision_positions,
    )
}

/// Print detailed hash collision statistics
fn print_collision_stats(reference: &[u8], k: usize, w: usize) {
    eprintln!("\n[Collision Analysis] Minimizer index collision statistics");
    eprintln!("============================================================");
    eprintln!(
        "Reference size: {} bp ({:.2} Gbp)",
        reference.len(),
        reference.len() as f64 / 1e9
    );
    eprintln!("Parameters: k={}, w={}", k, w);
    eprintln!();

    // geomap now uses full 64-bit k-mer encoding - NO collisions for k <= 32
    let kmer_bits = 2 * k;
    eprintln!("K-mer encoding: {} bits", kmer_bits);
    eprintln!("Index storage: 64 bits (full encoding preserved)");
    eprintln!();
    if kmer_bits <= 64 {
        eprintln!("✓ ZERO COLLISIONS: full 64-bit encoding used");
        eprintln!("  Each unique k-mer maps to exactly one index entry");
        eprintln!();
    }

    // Show what collisions WOULD have been with 32-bit truncation (for reference)
    eprintln!("--- Comparison: if 32-bit truncation were used ---");

    let (total, unique64, unique32, colliding, collision_pos) =
        count_hash_collisions(reference, k, w);

    let collision_rate = if unique64 > 0 {
        (unique64 - unique32) as f64 / unique64 as f64 * 100.0
    } else {
        0.0
    };

    let affected_rate = if total > 0 {
        collision_pos as f64 / total as f64 * 100.0
    } else {
        0.0
    };

    eprintln!("Total minimizer positions:     {:>12}", total);
    eprintln!("Unique 64-bit hashes:          {:>12}", unique64);
    eprintln!("Unique 32-bit hashes:          {:>12}", unique32);
    eprintln!(
        "Hash64 lost to truncation:     {:>12} ({:.4}%)",
        unique64 - unique32,
        collision_rate
    );
    eprintln!();
    eprintln!("Buckets with collisions:       {:>12}", colliding);
    eprintln!(
        "Positions affected:            {:>12} ({:.4}%)",
        collision_pos, affected_rate
    );
    eprintln!();

    // Theoretical analysis
    let n = unique64 as f64;
    let m = (1u64 << 32) as f64;
    let expected_collisions = n * (1.0 - (-n / m).exp());
    eprintln!(
        "Birthday paradox prediction:   {:>12.0} collisions",
        expected_collisions
    );
    eprintln!();

    if collision_rate < 0.01 {
        eprintln!("(32-bit would have: negligible collision rate)");
    } else if collision_rate < 1.0 {
        eprintln!(
            "(32-bit would have: {} collisions - chaining would filter)",
            unique64 - unique32
        );
    } else {
        eprintln!(
            "(32-bit would have: {} collisions - SIGNIFICANT impact)",
            unique64 - unique32
        );
    }
    eprintln!();
    eprintln!("✓ Using 64-bit: ALL {} unique k-mers preserved", unique64);
    eprintln!("============================================================\n");
}

// ============================================================================
// Seeding with minimizers
// ============================================================================

/// Ring buffer for monotonic deque - avoids VecDeque allocation overhead
/// Fixed size 64 is enough for any reasonable window size (w <= 50)
const RING_SIZE: usize = 64;
const RING_MASK: usize = RING_SIZE - 1;

/// Radix sort for (u64, u32) tuples, sorting by the u64 key
/// Uses 8-bit radix (256 buckets) with 4 passes on the lower 32 bits
/// This is faster than comparison sort for large arrays of integers
#[inline]
fn radix_sort_by_hash(data: &mut [(u64, u32)]) {
    if data.len() < 2 {
        return;
    }

    let mut aux = vec![(0u64, 0u32); data.len()];
    let mut counts = [0usize; 256];

    // 4 passes on lower 32 bits (usually enough for good distribution)
    for pass in 0..4 {
        let shift = pass * 8;

        // Count occurrences
        counts.fill(0);
        for &(hash, _) in data.iter() {
            let bucket = ((hash >> shift) & 0xFF) as usize;
            counts[bucket] += 1;
        }

        // Compute prefix sums (starting positions)
        let mut sum = 0;
        for c in counts.iter_mut() {
            let tmp = *c;
            *c = sum;
            sum += tmp;
        }

        // Scatter to aux array
        for &item in data.iter() {
            let bucket = ((item.0 >> shift) & 0xFF) as usize;
            aux[counts[bucket]] = item;
            counts[bucket] += 1;
        }

        // Swap data and aux
        data.copy_from_slice(&aux);
    }
}

/// Complement a single 2-bit encoded base: A(00)↔T(11), C(01)↔G(10)
#[inline(always)]
fn complement_bits(b: u64) -> u64 {
    b ^ 3
}

/// Per-platform seeding limits to balance speed/coverage
#[derive(Clone, Copy)]
struct SeedingLimits {
    read_hash_occ_limit: usize,
    hash_emit_cap: usize,
    max_anchors_per_kb: usize,
    min_dynamic_max_occ: usize,
}

#[inline]
fn select_seeding_limits(
    k: usize,
    w: usize,
    max_occ: usize,
    ultralong: bool,
    short_read_mode: bool,
    recovery_mode: bool,
) -> SeedingLimits {
    let _ = w;
    if short_read_mode {
        return SeedingLimits {
            read_hash_occ_limit: 6,
            hash_emit_cap: 12,
            max_anchors_per_kb: 800,
            min_dynamic_max_occ: max_occ.max(200),
        };
    }

    if k <= 16 {
        if recovery_mode {
            return SeedingLimits {
                read_hash_occ_limit: if ultralong { 256 } else { 192 },
                hash_emit_cap: if ultralong { 512 } else { 384 },
                max_anchors_per_kb: if ultralong { 900 } else { 720 },
                min_dynamic_max_occ: max_occ.max(400),
            };
        } else {
            return SeedingLimits {
                read_hash_occ_limit: if ultralong { 96 } else { 72 },
                hash_emit_cap: if ultralong { 256 } else { 192 },
                max_anchors_per_kb: if ultralong { 600 } else { 480 },
                min_dynamic_max_occ: max_occ.max(200),
            };
        }
    }

    if k == 17 {
        return SeedingLimits {
            read_hash_occ_limit: 20,
            hash_emit_cap: 48,
            max_anchors_per_kb: 220,
            min_dynamic_max_occ: max_occ.max(48),
        };
    }

    SeedingLimits {
        read_hash_occ_limit: 12,
        hash_emit_cap: 32,
        max_anchors_per_kb: 160,
        min_dynamic_max_occ: max_occ.max(32),
    }
}

/// Seed read using canonical k-mers (same as index)
/// TWO-PHASE VERSION: Collect minimizers, sort by hash, then batch lookup
/// This improves cache locality by accessing hash table in sorted order
/// max_occ controls filtering of repetitive k-mers (lower = more aggressive filtering)
#[inline]
fn seed_read<I: KmerIndex>(
    read: &[u8],
    ref_index: &I,
    k: usize,
    w: usize,
    max_occ: usize,
    limits: &SeedingLimits,
) -> Vec<Anchor> {
    if read.len() < k {
        return Vec::new();
    }

    // PHASE 1: Collect minimizers (hash, read_pos)
    let mut minimizers: Vec<(u64, u32)> = Vec::with_capacity(read.len() / w + 1);

    // Ring buffer for monotonic deque
    let mut ring: [(u64, u32); RING_SIZE] = [(0, 0); RING_SIZE];
    let mut head: usize = 0;
    let mut tail: usize = 0;
    let mut last_min_pos: usize = usize::MAX;
    let mut n_count = 0usize;
    let mask = (1u64 << (2 * k)) - 1;
    let shift = 2 * (k - 1);

    // Pre-count N's in first k-mer and compute initial encodings
    let mut kmer_fwd: u64 = 0;
    let mut kmer_rev: u64 = 0;
    for i in 0..k.min(read.len()) {
        if read[i] == b'N' {
            n_count += 1;
        }
        let base = base_to_bits(read[i]);
        kmer_fwd = (kmer_fwd << 2) | base;
        kmer_rev = (kmer_rev >> 2) | (complement_bits(base) << shift);
    }

    // Scan for minimizers
    for j in 0..=(read.len() - k) {
        if j > 0 {
            let old_base = read[j - 1];
            let new_base = read[j + k - 1];

            if old_base == b'N' {
                n_count -= 1;
            }
            if new_base == b'N' {
                n_count += 1;
            }

            let new_bits = base_to_bits(new_base);
            kmer_fwd = ((kmer_fwd << 2) | new_bits) & mask;
            kmer_rev = (kmer_rev >> 2) | (complement_bits(new_bits) << shift);
        }

        if n_count > 0 {
            head = 0;
            tail = 0;
            last_min_pos = usize::MAX;
            continue;
        }

        let j32 = j as u32;
        let canon = if kmer_fwd <= kmer_rev {
            kmer_fwd
        } else {
            kmer_rev
        };
        let hash = splitmix64(canon);

        // Monotonic deque operations
        while head != tail {
            let back_idx = (tail - 1) & RING_MASK;
            if ring[back_idx].0 >= hash {
                tail = back_idx;
            } else {
                break;
            }
        }

        ring[tail] = (hash, j32);
        tail = (tail + 1) & RING_MASK;

        let window_start = if j >= w { (j - w + 1) as u32 } else { 0 };
        while head != tail && ring[head].1 < window_start {
            head = (head + 1) & RING_MASK;
        }

        // Collect minimizer
        if j >= w - 1 && head != tail {
            let (min_hash, min_pos32) = ring[head];
            let min_pos = min_pos32 as usize;

            if min_pos != last_min_pos {
                last_min_pos = min_pos;
                minimizers.push((min_hash, min_pos32));
            }
        }
    }

    // PHASE 2: Sort minimizers by hash for cache-friendly lookup
    if minimizers.len() > 256 {
        radix_sort_by_hash(&mut minimizers);
    } else {
        minimizers.sort_unstable_by_key(|&(h, _)| h);
    }

    // PHASE 3: Batch lookup in sorted order with prefetching
    let mut anchors = Vec::with_capacity(minimizers.len() * 2);
    let k16 = k as u16;
    let n = minimizers.len();
    let mut dynamic_max_occ = max_occ;

    // Soft budget: avoid runaway anchor counts in repetitive reads
    let anchor_budget = ((read.len() + 999) / 1000) * limits.max_anchors_per_kb + 256;

    // Track per-hash reuse to avoid per-read repeats
    let mut current_hash: u64 = u64::MAX;
    let mut hash_run: usize = 0;
    let mut emitted_for_hash: usize = 0;

    // Process with lookahead prefetching (prefetch 8 entries ahead)
    const PREFETCH_DISTANCE: usize = 8;

    for i in 0..n {
        #[cfg(target_arch = "x86_64")]
        if i + PREFETCH_DISTANCE < n {
            ref_index.prefetch(minimizers[i + PREFETCH_DISTANCE].0);
        }

        let (hash, read_pos) = minimizers[i];

        if hash == current_hash {
            hash_run += 1;
        } else {
            current_hash = hash;
            hash_run = 0;
            emitted_for_hash = 0;
        }

        if hash_run >= limits.read_hash_occ_limit {
            continue;
        }

        if let Some(ref_positions) = ref_index.get_positions(hash) {
            let freq = ref_positions.len();
            if freq <= dynamic_max_occ {
                let stride = (freq / limits.hash_emit_cap).max(1);
                for &ref_pos in ref_positions.iter().step_by(stride) {
                    if emitted_for_hash >= limits.hash_emit_cap {
                        break;
                    }
                    anchors.push(Anchor {
                        read_start: read_pos,
                        ref_start: ref_pos,
                        len: k16,
                    });
                    emitted_for_hash += 1;

                    if anchors.len() >= anchor_budget
                        && dynamic_max_occ > limits.min_dynamic_max_occ
                    {
                        dynamic_max_occ = (dynamic_max_occ / 2).max(limits.min_dynamic_max_occ);
                    }
                }
            }
        }
    }

    anchors
}

/// Filter anchors to remove those landing on blacklisted regions
/// This removes anchors from problematic decoys like mitochondrial decoy GL000209.2
#[inline]
fn filter_blacklisted_anchors(
    anchors: Vec<Anchor>,
    blacklist_ranges: &[(usize, usize)],
) -> Vec<Anchor> {
    if blacklist_ranges.is_empty() {
        return anchors;
    }
    anchors
        .into_iter()
        .filter(|a| {
            let pos = a.ref_start as usize;
            !blacklist_ranges
                .iter()
                .any(|(start, end)| pos >= *start && pos < *end)
        })
        .collect()
}

// ============================================================================
// Optimized Chaining with Top-K extraction for MAPQ
// ============================================================================

/// Result of chaining: best chain, its score, and alternative chains for MAPQ
#[derive(Debug, Clone)]
struct ChainResult {
    best_chain: Vec<Anchor>,
    best_score: i32,
    /// Alternative chains (second-best, etc.) for MAPQ calculation
    alternatives: Vec<(Vec<Anchor>, i32, i32)>, // (chain, score, diagonal)
}

/// Apply primary chromosome preference as a tie-breaker.
/// When the best chain is on a decoy/alt contig and there's an alternative
/// on a primary chromosome with a similar score (within epsilon), prefer primary.
/// This is a "safe" tie-breaker that only affects ambiguous mappings.
///
/// For compressed mode, pass ref_pos_map to convert compressed positions to original.
/// Returns (best_chain, best_score, alternatives) - potentially swapped
fn apply_primary_contig_tiebreaker(
    mut chain_result: ChainResult,
    global_ref: &GlobalReference,
    ref_pos_map: Option<&[u32]>, // For compressed mode: maps compressed pos -> original pos
    score_epsilon: f32,          // e.g., 0.05 for 5% tolerance
) -> ChainResult {
    // If no alternatives or best chain is empty, nothing to do
    if chain_result.best_chain.is_empty() || chain_result.alternatives.is_empty() {
        return chain_result;
    }

    // Helper to convert position (handles compressed vs standard mode)
    let to_original_pos = |compressed_pos: usize| -> usize {
        match ref_pos_map {
            Some(map) if compressed_pos < map.len() => map[compressed_pos] as usize,
            _ => compressed_pos, // Standard mode: positions are already original
        }
    };

    // Check if best chain is on a primary contig
    let best_ref_pos = to_original_pos(chain_result.best_chain[0].ref_start as usize);
    let best_is_primary = global_ref.is_primary_position(best_ref_pos);

    // If already on primary, no need to change
    if best_is_primary {
        return chain_result;
    }

    // Best is on decoy/alt - look for a primary alternative with similar score
    let score_threshold = (chain_result.best_score as f32 * (1.0 - score_epsilon)) as i32;

    // Find best alternative on a primary contig
    let mut best_primary_idx: Option<usize> = None;
    let mut best_primary_score: i32 = 0;

    for (idx, (alt_chain, alt_score, _diag)) in chain_result.alternatives.iter().enumerate() {
        if alt_chain.is_empty() {
            continue;
        }

        // Check if this alternative is on a primary contig
        let alt_ref_pos = to_original_pos(alt_chain[0].ref_start as usize);
        let alt_is_primary = global_ref.is_primary_position(alt_ref_pos);

        if alt_is_primary && *alt_score >= score_threshold && *alt_score > best_primary_score {
            best_primary_idx = Some(idx);
            best_primary_score = *alt_score;
        }
    }

    // If we found a good primary alternative, swap it with the best
    if let Some(idx) = best_primary_idx {
        // Swap: move current best to alternatives, promote primary alternative to best
        let old_best_chain = std::mem::take(&mut chain_result.best_chain);
        let old_best_score = chain_result.best_score;
        let old_best_diag = if !old_best_chain.is_empty() {
            old_best_chain[0].ref_start as i32 - old_best_chain[0].read_start as i32
        } else {
            0
        };

        // Get the primary alternative
        let (new_best_chain, new_best_score, _) = chain_result.alternatives.remove(idx);

        // Update chain_result
        chain_result.best_chain = new_best_chain;
        chain_result.best_score = new_best_score;

        // Add old best to alternatives
        chain_result
            .alternatives
            .push((old_best_chain, old_best_score, old_best_diag));
    }

    chain_result
}

/// Filter anchors using diagonal histogram gating for ultralong reads.
/// Returns anchors near the top-K diagonal peaks only.
fn filter_anchors_by_diagonal_peaks(
    anchors: &[Anchor],
    bin_size: i32,
    num_peaks: usize,
    band: i32,
) -> Vec<Anchor> {
    if anchors.is_empty() {
        return Vec::new();
    }

    // Calculate diagonal for each anchor (in bp)
    let diags: Vec<i64> = anchors
        .iter()
        .map(|a| a.ref_start as i64 - a.read_start as i64)
        .collect();

    // Find min/max diagonal to size the histogram
    let min_diag = *diags.iter().min().unwrap();
    let max_diag = *diags.iter().max().unwrap();
    let diag_range = max_diag - min_diag;

    // If range is very small, keep all anchors
    if diag_range < bin_size as i64 * 2 {
        return anchors.to_vec();
    }

    // For human genome scale: use larger bins proportional to diagonal range
    // min 500bp, max 50000bp bins
    let effective_bin_size = if diag_range > 100_000_000 {
        50000i64 // Very large range: use 50kb bins
    } else if diag_range > 10_000_000 {
        10000i64 // Large range: use 10kb bins
    } else if diag_range > 1_000_000 {
        5000i64 // Medium range: use 5kb bins
    } else {
        bin_size as i64 // Small range: use default
    };

    // Build histogram
    let num_bins = (diag_range / effective_bin_size + 1) as usize;
    let mut histogram: Vec<usize> = vec![0; num_bins.min(100_000)]; // Cap at 100k bins

    for &diag in &diags {
        let bin = ((diag - min_diag) / effective_bin_size) as usize;
        if bin < histogram.len() {
            histogram[bin] += 1;
        }
    }

    // Find top-K peaks (bins with most anchors)
    let top_peaks: Vec<i64> = if num_peaks == 1 {
        // Fast path: O(n) linear scan for single peak (common case)
        let mut best_bin = 0usize;
        let mut best_count = 0usize;
        for (i, &count) in histogram.iter().enumerate() {
            if count > best_count {
                best_count = count;
                best_bin = i;
            }
        }
        if best_count > 0 {
            vec![min_diag + (best_bin as i64 * effective_bin_size) + (effective_bin_size / 2)]
        } else {
            Vec::new()
        }
    } else {
        // General case: partial sort for top-K peaks
        let mut bin_counts: Vec<(usize, usize)> = histogram
            .iter()
            .enumerate()
            .map(|(i, &count)| (i, count))
            .collect();
        if num_peaks < bin_counts.len() {
            bin_counts.select_nth_unstable_by(num_peaks - 1, |a, b| b.1.cmp(&a.1));
            bin_counts.truncate(num_peaks);
        }
        bin_counts
            .iter()
            .filter(|(_, count)| *count > 0)
            .map(|(bin, _)| min_diag + (*bin as i64 * effective_bin_size) + (effective_bin_size / 2))
            .collect()
    };

    if top_peaks.is_empty() {
        return anchors.to_vec();
    }

    // Use MUCH larger band for human genome scale
    // Goal: keep ~50-80% of anchors concentrated around top peaks
    // With 50kb bins, need at least 500kb band to capture enough anchors
    let effective_band = if diag_range > 100_000_000 {
        500_000i64 // Very large genome: 500kb band
    } else if diag_range > 10_000_000 {
        200_000i64 // Large range: 200kb band
    } else {
        (effective_bin_size * 10).max(band as i64) // 10 bins worth
    };

    // Filter anchors: keep only those within 'band' of any top peak
    let filtered: Vec<Anchor> = anchors
        .iter()
        .zip(diags.iter())
        .filter(|(_, &diag)| {
            top_peaks
                .iter()
                .any(|&peak| (diag - peak).abs() <= effective_band)
        })
        .map(|(anchor, _)| *anchor)
        .collect();

    // If filtering removed too many anchors (>90%), fall back to original
    // This handles cases where anchors are evenly distributed (no clear peaks)
    if filtered.len() * 10 < anchors.len() {
        return anchors.to_vec();
    }

    filtered
}

// ============================================================================
// Thread-local buffer pool for chain_anchors (eliminates repeated allocations)
// ============================================================================

struct ChainBufferPool {
    best_score: Vec<i32>,
    prev: Vec<i16>,
    endings: Vec<(usize, i32, i32)>,
    chain_indices: Vec<usize>,
}

impl ChainBufferPool {
    fn new() -> Self {
        Self {
            best_score: Vec::with_capacity(8192),
            prev: Vec::with_capacity(8192),
            endings: Vec::with_capacity(256),
            chain_indices: Vec::with_capacity(64),
        }
    }

    #[inline(always)]
    fn prepare(&mut self, n: usize) {
        // Clear and resize buffers - reuses existing capacity
        self.best_score.clear();
        self.prev.clear();
        self.endings.clear();

        // Ensure capacity without reallocating if already sufficient
        if self.best_score.capacity() < n {
            self.best_score.reserve(n - self.best_score.capacity());
        }
        if self.prev.capacity() < n {
            self.prev.reserve(n - self.prev.capacity());
        }
    }
}

thread_local! {
    static CHAIN_BUFFERS: RefCell<ChainBufferPool> = RefCell::new(ChainBufferPool::new());
}

/// Chain anchors and return top-K chains for topological MAPQ calculation
/// If ultralong=true, uses diagonal histogram gating and relaxed chaining
/// Uses thread-local buffer pools to avoid repeated allocations
fn chain_anchors_topk(
    anchors: &mut [Anchor],
    band_width: i32,
    top_k: usize,
    max_lookback: usize,
    gap_max: i32,
    gap_scale: i32,
    ultralong: bool,
) -> ChainResult {
    if anchors.is_empty() {
        return ChainResult {
            best_chain: Vec::new(),
            best_score: 0,
            alternatives: Vec::new(),
        };
    }

    // ULTRALONG MODE: Apply diagonal histogram gating FIRST
    let working_anchors: Vec<Anchor>;
    let anchors_slice: &mut [Anchor] = if ultralong {
        // Filter to keep only anchors near THE TOP diagonal peak (most anchors)
        // This ensures all anchors are collinear (same chromosome/region)
        // bin_size=500bp, num_peaks=1, band=2000bp around peak
        working_anchors = filter_anchors_by_diagonal_peaks(anchors, 500, 1, 2000);
        if working_anchors.is_empty() {
            return ChainResult {
                best_chain: Vec::new(),
                best_score: 0,
                alternatives: Vec::new(),
            };
        }
        // Can't mutate working_anchors directly, need to copy back
        anchors[..working_anchors.len()].copy_from_slice(&working_anchors);
        &mut anchors[..working_anchors.len()]
    } else {
        anchors
    };

    // Use shift instead of division for band calculation
    let band_shift = band_width.trailing_zeros();

    if ultralong {
        // Ultralong mode: sort by COARSE diagonal band (1MB) then read_start
        // This keeps anchors from same chromosome/region together
        let coarse_shift = 20; // 2^20 = ~1MB bands
        anchors_slice.sort_unstable_by(|a, b| {
            let diag_a = (a.ref_start as i64 - a.read_start as i64) >> coarse_shift;
            let diag_b = (b.ref_start as i64 - b.read_start as i64) >> coarse_shift;
            diag_a.cmp(&diag_b).then(a.read_start.cmp(&b.read_start))
        });
    } else {
        // Default mode: sort by (diagonal band, read_start) for strict chaining
        anchors_slice.sort_unstable_by(|a, b| {
            let diag_a = (a.ref_start as i32 - a.read_start as i32) >> band_shift;
            let diag_b = (b.ref_start as i32 - b.read_start as i32) >> band_shift;
            diag_a.cmp(&diag_b).then(a.read_start.cmp(&b.read_start))
        });
    }

    let n = anchors_slice.len();

    // Use thread-local buffer pool to avoid repeated allocations
    CHAIN_BUFFERS.with(|pool| {
        let mut pool = pool.borrow_mut();
        pool.prepare(n);

        // Initialize best_score with anchor lengths
        pool.best_score
            .extend(anchors_slice.iter().map(|a| a.len as i32));
        pool.prev.resize(n, -1);

        // Optimized DP with limited lookback
        // Compute diagonals on-the-fly instead of pre-allocating Vec
        for i in 1..n {
            let ai = &anchors_slice[i];
            let ai_diag = (ai.ref_start as i32 - ai.read_start as i32) >> band_shift;
            let ai_read_start = ai.read_start;
            let ai_ref_start = ai.ref_start;
            let ai_len = ai.len as i32;

            let start_j = if i > max_lookback {
                i - max_lookback
            } else {
                0
            };

            for j in (start_j..i).rev() {
                let aj = &anchors_slice[j];

                if ultralong {
                    // Ultralong mode: allow chaining within same 1MB diagonal band
                    let coarse_shift = 20;
                    let diag_i = (ai_ref_start as i64 - ai_read_start as i64) >> coarse_shift;
                    let diag_j = (aj.ref_start as i64 - aj.read_start as i64) >> coarse_shift;
                    if diag_i != diag_j {
                        break; // Sorted by diagonal, early termination
                    }
                } else {
                    // Strict mode: only chain within same diagonal band
                    let aj_diag = (aj.ref_start as i32 - aj.read_start as i32) >> band_shift;
                    if aj_diag != ai_diag {
                        continue;
                    }
                }

                let j_end_read = aj.read_start + aj.len as u32;
                let j_end_ref = aj.ref_start + aj.len as u32;

                if j_end_read <= ai_read_start && j_end_ref <= ai_ref_start {
                    let gap_read = ai_read_start - j_end_read;
                    let gap_ref = ai_ref_start - j_end_ref;

                    // For wide-band chaining (SV-aware mode), reject connections across
                    // extremely long reference gaps. 50kb max covers the vast majority
                    // of real SVs while preventing false chains across distant regions.
                    if band_width >= 8192 && gap_ref > 50_000 {
                        continue;
                    }

                    let gap_diff = (gap_read as i32 - gap_ref as i32).unsigned_abs() as i32;
                    // For wide-band mode, use sqrt-scaled penalty so large diagonal shifts
                    // (SVs) are properly penalized while small SVs remain cheap.
                    // 100bp→5, 500bp→11, 1kb→16, 5kb→35, 10kb→50, 50kb→112
                    let gap_penalty = if band_width >= 8192 && gap_diff > 100 {
                        ((gap_diff as f64).sqrt() / 2.0) as i32
                    } else {
                        gap_diff.min(gap_max) / gap_scale
                    };

                    let candidate = pool.best_score[j] + ai_len - gap_penalty;
                    if candidate > pool.best_score[i] {
                        pool.best_score[i] = candidate;
                        pool.prev[i] = j as i16;
                    }
                }
            }
        }

        // Find top-K ending positions (different diagonals = different chains)
        // Use index-based iteration to avoid borrow conflicts
        for i in 0..n {
            let s = pool.best_score[i];
            let diag = (anchors_slice[i].ref_start as i32 - anchors_slice[i].read_start as i32)
                / band_width;
            pool.endings.push((i, s, diag));
        }

        // Use partial sort - only need top elements, not full sort
        let k_limit = (top_k * 3).min(pool.endings.len());
        if k_limit > 0 && k_limit < pool.endings.len() {
            pool.endings
                .select_nth_unstable_by(k_limit - 1, |a, b| b.1.cmp(&a.1));
            pool.endings.truncate(k_limit);
        }
        pool.endings.sort_by(|a, b| b.1.cmp(&a.1));

        // Copy endings to local vec to avoid borrow conflicts when using chain_indices
        let endings_copy: Vec<(usize, i32, i32)> = pool.endings.clone();

        // Extract top-K chains from different diagonals
        // Also track reference positions to detect segmental duplications (>1MB apart)
        let mut chains: Vec<(Vec<Anchor>, i32, i32)> = Vec::with_capacity(top_k + 5);
        let mut used_diagonals: Vec<i32> = Vec::with_capacity(top_k + 5);
        let mut used_ref_starts: Vec<u32> = Vec::with_capacity(top_k + 5);
        const SEGDUP_THRESHOLD: u32 = 1_000_000; // 1MB

        for &(end_idx, score, diag) in &endings_copy {
            // Backtrack to get the chain (reuse chain_indices buffer)
            pool.chain_indices.clear();
            let mut current = end_idx as i16;
            while current >= 0 {
                pool.chain_indices.push(current as usize);
                current = pool.prev[current as usize];
            }
            pool.chain_indices.reverse();

            let chain: Vec<Anchor> = pool
                .chain_indices
                .iter()
                .map(|&i| anchors_slice[i])
                .collect();

            if chain.is_empty() {
                continue;
            }

            let chain_ref_start = chain[0].ref_start;

            // Check if dominated by existing chain
            // A chain is dominated ONLY if:
            // 1. It's on a similar diagonal (within 2 bands), AND
            // 2. It's at a similar reference position (within 1MB)
            // This ensures segmental duplications are kept as alternatives
            let dominated = used_diagonals
                .iter()
                .zip(used_ref_starts.iter())
                .any(|(&d, &r)| {
                    let diag_close = (d - diag).abs() <= 2;
                    let ref_close = if chain_ref_start > r {
                        chain_ref_start - r < SEGDUP_THRESHOLD
                    } else {
                        r - chain_ref_start < SEGDUP_THRESHOLD
                    };
                    diag_close && ref_close
                });

            if dominated && !chains.is_empty() {
                continue;
            }

            chains.push((chain, score, diag));
            used_diagonals.push(diag);
            used_ref_starts.push(chain_ref_start);

            if chains.len() >= top_k + 5 {
                break;
            }
        }

        if chains.is_empty() {
            return ChainResult {
                best_chain: Vec::new(),
                best_score: 0,
                alternatives: Vec::new(),
            };
        }

        // TIE-BREAKER: When top chains have similar scores (within 5%),
        // prefer the one with more anchors (better coverage of the read)
        // This helps select the correct copy in segmental duplications
        let top_score = chains[0].1;
        let score_threshold = top_score * 95 / 100; // 95% of top score

        // Find the chain with most anchors among those with similar scores
        let mut best_idx = 0;
        let mut best_anchor_count = chains[0].0.len();

        for (idx, (chain, score, _)) in chains.iter().enumerate().skip(1) {
            if *score < score_threshold {
                break; // Chains are sorted by score, so we can stop early
            }
            if chain.len() > best_anchor_count {
                best_idx = idx;
                best_anchor_count = chain.len();
            }
        }

        // If best_idx changed, swap the chains
        if best_idx > 0 {
            chains.swap(0, best_idx);
        }

        let (best_chain, best_score, _best_diag) = chains.remove(0);

        ChainResult {
            best_chain,
            best_score,
            alternatives: chains,
        }
    })
}

// ============================================================================
// TOPOLOGICAL MAPQ CALCULATION
// ============================================================================
//
// Novel approach inspired by ALTA: Calculate MAPQ from information-theoretic
// limits rather than arbitrary score differences.
//
// For k topologically equivalent alignments (e.g., tandem repeats):
// - P(correct) = 1/k (best any algorithm can do)
// - P(wrong) = (k-1)/k
// - MAPQ = -10 * log10(P(wrong))
//
// This correctly assigns low MAPQ to ambiguous tandem repeat regions where
// no algorithm can determine the true position.

/// Classification of alignment ambiguity
#[derive(Debug, Clone, Copy, PartialEq)]
enum TopologicalClass {
    /// Single unambiguous alignment
    Unique,
    /// Best is clearly better than alternatives (score ratio > 1.5)
    WeaklyAmbiguous,
    /// Multiple nearly-equivalent alignments (tandem repeat scenario)
    StronglyAmbiguous { n_equivalent: usize },
}

/// Check if two chains are topologically equivalent (same shape, different position)
/// This indicates a tandem repeat where the read could map equally well to either copy
///
/// STRICT criteria for true topological equivalence:
/// 1. Scores must be nearly identical (within 5%)
/// 2. Anchor counts must be similar
/// 3. Read positions must be nearly identical (same path shape)
/// 4. Reference shift must be constant across all anchors (pure translation)
fn chains_are_equivalent(chain1: &[Anchor], chain2: &[Anchor], score1: i32, score2: i32) -> bool {
    // Scores must be very similar (within 5%) - stricter than before
    let score_ratio = score1.max(score2) as f64 / score1.min(score2).max(1) as f64;
    if score_ratio > 1.05 {
        return false;
    }

    // Chains must have similar number of anchors (within 20%)
    if chain1.len() < 2 || chain2.len() < 2 {
        return false; // Need at least 2 anchors to verify pattern
    }

    let len_ratio = chain1.len().max(chain2.len()) as f64 / chain1.len().min(chain2.len()) as f64;
    if len_ratio > 1.2 {
        return false;
    }

    // Check if it's a pure translation (constant reference shift)
    // This is the hallmark of tandem repeats

    // Calculate reference shift between first anchors
    let ref_shift = chain2[0].ref_start as i64 - chain1[0].ref_start as i64;

    // Reference shift must be non-trivial (otherwise they're the same chain)
    if ref_shift.abs() < 50 {
        return false;
    }

    // For tandem repeats, read positions should be nearly identical,
    // ref positions should be shifted by a constant amount
    let n_check = chain1.len().min(chain2.len()).min(10);
    let mut n_consistent = 0;

    for i in 0..n_check {
        // Read positions should be nearly identical (same path shape)
        // Allow only small variation due to k-mer offsets
        let read_diff = (chain2[i].read_start as i64 - chain1[i].read_start as i64).abs();
        if read_diff > 20 {
            continue; // This anchor pair doesn't match
        }

        // Reference shift should be constant (pure translation)
        let this_ref_shift = chain2[i].ref_start as i64 - chain1[i].ref_start as i64;
        if (this_ref_shift - ref_shift).abs() > 30 {
            continue; // This anchor pair doesn't match
        }

        n_consistent += 1;
    }

    // At least 60% of checked anchors must show consistent translation pattern
    n_consistent as f64 / n_check as f64 >= 0.6
}

/// Classify the topological ambiguity of an alignment
fn classify_ambiguity(chain_result: &ChainResult, read_len: usize) -> TopologicalClass {
    if chain_result.alternatives.is_empty() {
        return TopologicalClass::Unique;
    }

    let best_score = chain_result.best_score;
    let best_chain = &chain_result.best_chain;

    // Count how many alternatives are topologically equivalent (true tandem repeats)
    let mut n_equivalent = 1; // Include the best chain itself

    for (alt_chain, alt_score, _alt_diag) in &chain_result.alternatives {
        if chains_are_equivalent(best_chain, alt_chain, best_score, *alt_score) {
            n_equivalent += 1;
        }
    }

    if n_equivalent > 1 {
        return TopologicalClass::StronglyAmbiguous { n_equivalent };
    }

    // No true tandem repeat detected - check score difference
    let second_best_score = chain_result
        .alternatives
        .first()
        .map(|(_, s, _)| *s)
        .unwrap_or(0);

    // If second best is much worse, we're confident
    // Use ratio of (best - second) / best as confidence metric
    let score_margin = (best_score - second_best_score) as f64 / best_score.max(1) as f64;

    // Require larger margin for short reads (score differences are smaller, fewer anchors)
    let unique_threshold = if read_len < 300 { 0.30 } else { 0.20 };

    if score_margin >= unique_threshold {
        TopologicalClass::Unique
    } else {
        TopologicalClass::WeaklyAmbiguous
    }
}

/// Calculate MAPQ using topological analysis
///
/// Key insight: For k topologically equivalent positions (tandem repeats),
/// no algorithm can achieve better than 1/k accuracy. We calculate MAPQ
/// from this fundamental limit, not from arbitrary score differences.
///
/// MAPQ formula follows SAM spec: MAPQ = -10 * log10(P(wrong))
/// - Unique: high MAPQ based on coverage (up to 60)
/// - WeaklyAmbiguous: moderate MAPQ based on score margin (10-30)
/// - StronglyAmbiguous: low MAPQ from information-theoretic limit (0-10)
fn calculate_topological_mapq(chain_result: &ChainResult, read_len: usize) -> u8 {
    let class = classify_ambiguity(chain_result, read_len);
    let is_short = read_len < 300;

    match class {
        TopologicalClass::Unique => {
            // High confidence - use coverage-based MAPQ
            // Coverage of chain score relative to read length
            let coverage = chain_result.best_score as f64 / read_len as f64;

            // Scale to max 60 (50 for short reads), with reasonable floor at ~30 for good alignments
            let max_mapq = if is_short { 50.0 } else { 60.0 };
            let base_mapq = 30.0 + (coverage * 30.0);
            base_mapq.min(max_mapq) as u8
        }

        TopologicalClass::WeaklyAmbiguous => {
            // Second best exists and is competitive but not equivalent
            let second_score = chain_result
                .alternatives
                .first()
                .map(|(_, s, _)| *s)
                .unwrap_or(0);

            // Score margin determines confidence
            let margin = (chain_result.best_score - second_score) as f64
                / chain_result.best_score.max(1) as f64;

            // Map margin to MAPQ (more conservative for short reads)
            let (scale, max_mapq) = if is_short {
                (0.3, 25.0) // margin 0-0.3 → MAPQ 10-25
            } else {
                (0.2, 30.0) // margin 0-0.2 → MAPQ 10-30
            };
            let base_mapq = 10.0 + (margin / scale) * (max_mapq - 10.0);
            base_mapq.min(max_mapq) as u8
        }

        TopologicalClass::StronglyAmbiguous { n_equivalent } => {
            // CRITICAL: For k equivalent alignments, P(wrong) = (k-1)/k
            // MAPQ = -10 * log10(P(wrong))
            //
            // k=2: P(wrong)=0.5  → MAPQ = 3
            // k=3: P(wrong)=0.67 → MAPQ = 1.8
            // k=5: P(wrong)=0.8  → MAPQ = 1.0
            // k=10: P(wrong)=0.9 → MAPQ = 0.5
            let k = n_equivalent as f64;
            let p_wrong = (k - 1.0) / k;

            if p_wrong >= 1.0 {
                return 0;
            }

            let mapq = -10.0 * p_wrong.log10();
            mapq.max(0.0).min(10.0) as u8 // Cap at 10 for tandem repeats
        }
    }
}

/// Apply MAPQ penalty for multi-mapping scenarios:
/// 1. Cross-chromosome alternatives (different chromosomes)
/// 2. Same-chromosome segmental duplications (>1MB apart on same chromosome)
///
/// For CLR reads especially, segmental duplications can cause high-scoring
/// alternative alignments on the same chromosome at very different positions.
fn apply_cross_chromosome_penalty(
    mapq: u8,
    chain_result: &ChainResult,
    global_ref: &GlobalReference,
    ref_pos_map: Option<&[u32]>, // For compressed mode
) -> u8 {
    if chain_result.alternatives.is_empty() {
        return mapq;
    }

    let best_score = chain_result.best_score;
    if best_score == 0 {
        return mapq;
    }

    // Get chromosome and position of best chain
    let best_ref_pos = chain_result
        .best_chain
        .first()
        .map(|a| a.ref_start as usize)
        .unwrap_or(0);
    let best_original_pos = if let Some(pos_map) = ref_pos_map {
        if best_ref_pos < pos_map.len() {
            pos_map[best_ref_pos] as usize
        } else {
            best_ref_pos
        }
    } else {
        best_ref_pos
    };
    let (best_chrom_idx, best_local_pos) = global_ref.global_to_local(best_original_pos);

    // Check for multi-mapping scenarios:
    // 1. Cross-chromosome alternatives
    // 2. Same-chromosome alternatives at very different positions (segmental duplications)
    let mut best_ambiguous_score = 0i32;
    const SEGDUP_DISTANCE: usize = 1_000_000; // 1MB threshold for segmental duplications

    for (alt_chain, alt_score, _) in &chain_result.alternatives {
        // Score must be competitive (within 30% of best)
        // Lower threshold to catch segmental duplications with degraded scores
        if *alt_score < best_score * 3 / 10 {
            continue;
        }

        let alt_ref_pos = alt_chain.first().map(|a| a.ref_start as usize).unwrap_or(0);
        let alt_original_pos = if let Some(pos_map) = ref_pos_map {
            if alt_ref_pos < pos_map.len() {
                pos_map[alt_ref_pos] as usize
            } else {
                alt_ref_pos
            }
        } else {
            alt_ref_pos
        };
        let (alt_chrom_idx, alt_local_pos) = global_ref.global_to_local(alt_original_pos);

        // Case 1: Different chromosome
        if alt_chrom_idx != best_chrom_idx {
            best_ambiguous_score = best_ambiguous_score.max(*alt_score);
            continue;
        }

        // Case 2: Same chromosome but >1MB apart (segmental duplication)
        let distance = if alt_local_pos > best_local_pos {
            alt_local_pos - best_local_pos
        } else {
            best_local_pos - alt_local_pos
        };

        if distance > SEGDUP_DISTANCE {
            // This is a segmental duplication scenario - treat similarly to cross-chrom
            best_ambiguous_score = best_ambiguous_score.max(*alt_score);
        }
    }

    if best_ambiguous_score == 0 {
        return mapq;
    }

    // Ambiguous alternative exists - apply penalty based on score ratio
    let score_ratio = best_ambiguous_score as f64 / best_score as f64;

    // More aggressive MAPQ capping for ambiguous alignments
    // If alternative is very close (>80% of best), cap MAPQ at 3 (essentially random)
    // If alternative is close (>60% of best), cap MAPQ at 10
    // If alternative is moderate (>40% of best), cap MAPQ at 20
    // If alternative is weak (>30% of best), cap MAPQ at 25
    let max_mapq = if score_ratio > 0.8 {
        3
    } else if score_ratio > 0.6 {
        10
    } else if score_ratio > 0.4 {
        20
    } else {
        25
    };

    mapq.min(max_mapq)
}

// ============================================================================
// Banded alignment for CIGAR generation
// ============================================================================

/// Fast gap alignment for CIGAR generation
/// Uses simplified approach: for small gaps, compute edit operations directly
/// For larger gaps, use diagonal-following heuristic
#[inline]
fn fast_gap_align(query: &[u8], target: &[u8]) -> Vec<(char, usize)> {
    let n = query.len();
    let m = target.len();

    // Empty sequences
    if n == 0 && m == 0 {
        return vec![];
    }
    if n == 0 {
        return vec![('D', m)];
    }
    if m == 0 {
        return vec![('I', n)];
    }

    // For very small gaps, use simple linear comparison
    if n <= 10 && m <= 10 {
        return simple_align(query, target);
    }

    // For medium gaps, use diagonal-following heuristic
    if n <= 500 && m <= 500 {
        return diagonal_align(query, target);
    }

    // For large gaps, just report the length difference as indel + matches
    let min_len = n.min(m);
    let mut ops = Vec::with_capacity(3);
    if n > m {
        ops.push(('M', min_len));
        ops.push(('I', n - m));
    } else if m > n {
        ops.push(('M', min_len));
        ops.push(('D', m - n));
    } else {
        ops.push(('M', n));
    }
    ops
}

/// Simple alignment for very short sequences - just count differences
#[inline]
fn simple_align(query: &[u8], target: &[u8]) -> Vec<(char, usize)> {
    let n = query.len();
    let m = target.len();

    if n == m {
        // Same length: just report as matches
        return vec![('M', n)];
    }

    // Different lengths: report shorter as matches, difference as indel
    let min_len = n.min(m);
    let mut ops = Vec::with_capacity(2);
    ops.push(('M', min_len));
    if n > m {
        ops.push(('I', n - m));
    } else {
        ops.push(('D', m - n));
    }
    ops
}

/// Diagonal-following alignment - O(n) heuristic
/// Follows the diagonal, handling indels when sequences diverge
#[inline]
fn diagonal_align(query: &[u8], target: &[u8]) -> Vec<(char, usize)> {
    let n = query.len();
    let m = target.len();

    let mut ops: Vec<(char, usize)> = Vec::with_capacity(16);
    let mut qi = 0usize;
    let mut ti = 0usize;

    while qi < n && ti < m {
        // Count matching bases
        let mut match_len = 0;
        while qi + match_len < n
            && ti + match_len < m
            && query[qi + match_len] == target[ti + match_len]
        {
            match_len += 1;
        }

        if match_len > 0 {
            add_cigar_op_fast(&mut ops, 'M', match_len);
            qi += match_len;
            ti += match_len;
            continue;
        }

        // Mismatch: look ahead to find best recovery
        // Try: 1 mismatch, 1 insertion, 1 deletion
        let remaining_q = n - qi;
        let remaining_t = m - ti;

        if remaining_q == 0 || remaining_t == 0 {
            break;
        }

        // Simple heuristic: if one sequence is longer, prefer that indel
        if remaining_q > remaining_t + 2 {
            // Query is longer - insertion
            add_cigar_op_fast(&mut ops, 'I', 1);
            qi += 1;
        } else if remaining_t > remaining_q + 2 {
            // Target is longer - deletion
            add_cigar_op_fast(&mut ops, 'D', 1);
            ti += 1;
        } else {
            // Similar lengths - treat as mismatch
            add_cigar_op_fast(&mut ops, 'M', 1);
            qi += 1;
            ti += 1;
        }
    }

    // Handle remaining bases
    if qi < n {
        add_cigar_op_fast(&mut ops, 'I', n - qi);
    }
    if ti < m {
        add_cigar_op_fast(&mut ops, 'D', m - ti);
    }

    ops
}

/// Refined diagonal alignment with lookahead-based indel detection - O(n), DP-free
///
/// Like `diagonal_align()` but uses a short lookahead window to detect small indels
/// (1-4bp) even in equal-length segments. When a mismatch is encountered, it tries:
/// - Mismatch (substitution): skip 1 base in both sequences
/// - Insertion of 1-4bp: skip bases in query, check if sequences realign
/// - Deletion of 1-4bp: skip bases in target, check if sequences realign
/// Picks the option with the most matching bases in the lookahead window.
///
/// Complexity: O(n * MAX_TRY * LOOKAHEAD) = O(n) with small constants
#[inline]
fn diagonal_align_refined(query: &[u8], target: &[u8], max_indel: usize) -> Vec<(char, usize)> {
    let n = query.len();
    let m = target.len();

    if n == 0 && m == 0 {
        return vec![];
    }
    if n == 0 {
        return vec![('D', m)];
    }
    if m == 0 {
        return vec![('I', n)];
    }

    // For very short sequences, use simple comparison
    if n <= 3 && m <= 3 {
        return simple_align(query, target);
    }

    const LOOKAHEAD: usize = 8; // Check next 8 bases for match quality

    let mut ops: Vec<(char, usize)> = Vec::with_capacity(32);
    let mut qi = 0usize;
    let mut ti = 0usize;

    while qi < n && ti < m {
        // Count exact matches on the diagonal
        let mut match_len = 0;
        while qi + match_len < n
            && ti + match_len < m
            && query[qi + match_len] == target[ti + match_len]
        {
            match_len += 1;
        }

        if match_len > 0 {
            add_cigar_op_fast(&mut ops, 'M', match_len);
            qi += match_len;
            ti += match_len;
            continue;
        }

        // Mismatch: evaluate candidates using lookahead
        // Score = number of matching bases in lookahead window after the operation
        // We want to pick the operation that best re-synchronizes the sequences

        // Candidate 1: Mismatch (substitution) - skip 1 in both
        let (mut best_score, mut best_op, mut best_qi_adv, mut best_ti_adv) = {
            let mut score = 0i32;
            for k in 0..LOOKAHEAD {
                if qi + 1 + k < n && ti + 1 + k < m
                    && query[qi + 1 + k] == target[ti + 1 + k]
                {
                    score += 1;
                }
            }
            (score, b'M', 1usize, 1usize)
        };

        // Candidate 2: Insertions (1..MAX_TRY bp in query)
        for gap in 1..=max_indel {
            if qi + gap >= n {
                break;
            }
            let mut score = 0i32;
            for k in 0..LOOKAHEAD {
                if qi + gap + k < n && ti + k < m
                    && query[qi + gap + k] == target[ti + k]
                {
                    score += 1;
                }
            }
            // Prefer smaller gaps when scores are equal (penalize by gap size)
            if score > best_score + (gap as i32 - 1) {
                best_score = score;
                best_op = b'I';
                best_qi_adv = gap;
                best_ti_adv = 0;
            }
        }

        // Candidate 3: Deletions (1..MAX_TRY bp in target)
        for gap in 1..=max_indel {
            if ti + gap >= m {
                break;
            }
            let mut score = 0i32;
            for k in 0..LOOKAHEAD {
                if qi + k < n && ti + gap + k < m
                    && query[qi + k] == target[ti + gap + k]
                {
                    score += 1;
                }
            }
            if score > best_score + (gap as i32 - 1) {
                best_score = score;
                best_op = b'D';
                best_qi_adv = 0;
                best_ti_adv = gap;
            }
        }

        // Apply the best operation
        match best_op {
            b'M' => {
                add_cigar_op_fast(&mut ops, 'M', 1);
                qi += best_qi_adv;
                ti += best_ti_adv;
            }
            b'I' => {
                add_cigar_op_fast(&mut ops, 'I', best_qi_adv);
                qi += best_qi_adv;
            }
            b'D' => {
                add_cigar_op_fast(&mut ops, 'D', best_ti_adv);
                ti += best_ti_adv;
            }
            _ => unreachable!(),
        }
    }

    // Handle remaining bases
    if qi < n {
        add_cigar_op_fast(&mut ops, 'I', n - qi);
    }
    if ti < m {
        add_cigar_op_fast(&mut ops, 'D', m - ti);
    }

    ops
}

#[inline]
fn add_cigar_op_fast(ops: &mut Vec<(char, usize)>, op: char, len: usize) {
    if len == 0 {
        return;
    }
    if let Some(last) = ops.last_mut() {
        if last.0 == op {
            last.1 += len;
            return;
        }
    }
    ops.push((op, len));
}

/// Optimal affine-gap (Gotoh) global alignment of two SHORT segments, returning CIGAR ops.
/// Used by the optional DP-polish mode (`--dp`) in place of the greedy `diagonal_align_refined`,
/// to place indels/mismatches optimally between anchors. DP is confined to short inter-anchor
/// gaps (dense anchors → cheap), keeping overall alignment fast. Falls back to the greedy path
/// for over-long segments (caller-enforced). Costs: match 0, mismatch 4, gap-open 6, extend 1.
fn nw_affine_cigar(query: &[u8], target: &[u8]) -> Vec<(char, usize)> {
    let n = query.len();
    let m = target.len();
    if n == 0 && m == 0 {
        return Vec::new();
    }
    if n == 0 {
        return vec![('D', m)];
    }
    if m == 0 {
        return vec![('I', n)];
    }
    const MM: i32 = 4; // mismatch cost
    const GO: i32 = 6; // gap open (added on the first gap base)
    const GE: i32 = 1; // gap extend (per gap base)
    const INF: i32 = 1 << 28;
    let w = m + 1;
    let idx = |i: usize, j: usize| i * w + j;
    // mm = best ending in match/mismatch; ix = ending in insertion (gap in target, consumes
    // query); iy = ending in deletion (gap in query, consumes target).
    let mut mm = vec![INF; (n + 1) * w];
    let mut ix = vec![INF; (n + 1) * w];
    let mut iy = vec![INF; (n + 1) * w];
    mm[idx(0, 0)] = 0;
    for i in 1..=n {
        ix[idx(i, 0)] = GO + GE * i as i32;
    }
    for j in 1..=m {
        iy[idx(0, j)] = GO + GE * j as i32;
    }
    for i in 1..=n {
        for j in 1..=m {
            let sub = if query[i - 1] == target[j - 1] { 0 } else { MM };
            let diag = mm[idx(i - 1, j - 1)]
                .min(ix[idx(i - 1, j - 1)])
                .min(iy[idx(i - 1, j - 1)]);
            mm[idx(i, j)] = diag.saturating_add(sub);
            ix[idx(i, j)] = (mm[idx(i - 1, j)] + GO + GE).min(ix[idx(i - 1, j)] + GE);
            iy[idx(i, j)] = (mm[idx(i, j - 1)] + GO + GE).min(iy[idx(i, j - 1)] + GE);
        }
    }
    // Traceback from the cheapest of the three final states.
    let mut state = {
        let (mut s, mut best) = (0u8, mm[idx(n, m)]);
        if ix[idx(n, m)] < best {
            s = 1;
            best = ix[idx(n, m)];
        }
        if iy[idx(n, m)] < best {
            s = 2;
        }
        s
    };
    let (mut i, mut j) = (n, m);
    let mut rev: Vec<(char, usize)> = Vec::with_capacity(16);
    while i > 0 || j > 0 {
        match state {
            0 => {
                // came into mm[i][j] from diag (min of the three at i-1,j-1)
                let sub = if query[i - 1] == target[j - 1] { 0 } else { MM };
                let target_cost = mm[idx(i, j)] - sub;
                add_cigar_op_fast(&mut rev, 'M', 1);
                let (pi, pj) = (i - 1, j - 1);
                state = if mm[idx(pi, pj)] == target_cost {
                    0
                } else if ix[idx(pi, pj)] == target_cost {
                    1
                } else {
                    2
                };
                i = pi;
                j = pj;
            }
            1 => {
                // insertion: consumes query
                let from_open = mm[idx(i - 1, j)] + GO + GE;
                add_cigar_op_fast(&mut rev, 'I', 1);
                state = if ix[idx(i, j)] == from_open { 0 } else { 1 };
                i -= 1;
            }
            _ => {
                // deletion: consumes target
                let from_open = mm[idx(i, j - 1)] + GO + GE;
                add_cigar_op_fast(&mut rev, 'D', 1);
                state = if iy[idx(i, j)] == from_open { 0 } else { 2 };
                j -= 1;
            }
        }
    }
    rev.reverse();
    // merge adjacent same-op runs (add_cigar_op_fast already merged consecutive pushes,
    // but reversal can juxtapose equal ops across the boundary of runs)
    let mut out: Vec<(char, usize)> = Vec::with_capacity(rev.len());
    for (op, len) in rev {
        if let Some(last) = out.last_mut() {
            if last.0 == op {
                last.1 += len;
                continue;
            }
        }
        out.push((op, len));
    }
    out
}

/// Global flag: when set (`--dp`), the polish gap aligner uses optimal affine-gap DP
/// (`nw_affine_cigar`) on short segments instead of the greedy heuristic. Set once in main().
static DP_POLISH: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

/// Align a gap/segment for CIGAR polishing: optimal affine-gap DP when `--dp` is enabled and
/// the segment is short enough (DP confined to short inter-anchor gaps keeps alignment fast),
/// otherwise the greedy heuristic (default; behaviour unchanged when `--dp` is off).
#[inline]
fn polish_gap(query: &[u8], target: &[u8], max_indel: usize) -> Vec<(char, usize)> {
    if DP_POLISH.load(std::sync::atomic::Ordering::Relaxed)
        && query.len() <= 512
        && target.len() <= 512
    {
        nw_affine_cigar(query, target)
    } else {
        diagonal_align_refined(query, target, max_indel)
    }
}

#[inline]
fn add_cigar_op(ops: &mut Vec<(char, usize)>, op: char, len: usize) {
    if let Some(last) = ops.last_mut() {
        if last.0 == op {
            last.1 += len;
            return;
        }
    }
    ops.push((op, len));
}

/// Build CIGAR string from ops using itoa for fast integer conversion
/// 1.5-2x faster than write! macro for CIGAR generation
#[inline]
fn cigar_to_string(ops: &[(char, usize)]) -> String {
    // Calculate exact capacity: number of digits + number of ops
    let capacity: usize = ops
        .iter()
        .map(|(_, len)| {
            if *len < 10 {
                2
            } else if *len < 100 {
                3
            } else if *len < 1000 {
                4
            } else if *len < 10000 {
                5
            } else {
                6
            }
        })
        .sum();

    let mut s = String::with_capacity(capacity);
    let mut buf = itoa::Buffer::new();

    for &(op, len) in ops {
        s.push_str(buf.format(len));
        s.push(op);
    }
    s
}

/// Align between anchors and build full CIGAR
/// Now also extends alignment to read boundaries
fn align_with_anchors_cigar(
    anchors: &[Anchor],
    read: &[u8],
    reference: &[u8],
    chain_score: i32,
    mapq: u8,
    strand: char,
    do_cigar: bool,
) -> Option<AlignmentResult> {
    if anchors.is_empty() {
        return None;
    }

    let first = anchors.first().unwrap();
    let last = anchors.last().unwrap();

    let anchor_ref_start = first.ref_start as usize;
    let anchor_ref_end = (last.ref_start + last.len as u32) as usize;
    let anchor_read_start = first.read_start as usize;
    let anchor_read_end = (last.read_start + last.len as u32) as usize;

    // Extend alignment to read boundaries
    // Left extension: how much read sequence before first anchor?
    let left_extend_read = anchor_read_start;
    let left_extend_ref = left_extend_read.min(anchor_ref_start); // don't go past reference start

    // Right extension: how much read sequence after last anchor?
    let right_extend_read = read.len().saturating_sub(anchor_read_end);
    let right_extend_ref = right_extend_read.min(reference.len().saturating_sub(anchor_ref_end));

    // Extended coordinates
    let ref_start = anchor_ref_start.saturating_sub(left_extend_ref);
    let ref_end = (anchor_ref_end + right_extend_ref).min(reference.len());
    let read_start = anchor_read_start.saturating_sub(left_extend_read);
    let read_end = (anchor_read_end + right_extend_read).min(read.len());

    let cigar = if do_cigar && ref_end <= reference.len() && read_end <= read.len() {
        let mut all_ops: Vec<(char, usize)> = Vec::new();

        // Left extension (before first anchor)
        if left_extend_read > 0 || left_extend_ref > 0 {
            let ext_read = &read[read_start..anchor_read_start];
            let ext_ref = &reference[ref_start..anchor_ref_start];
            if !ext_read.is_empty() || !ext_ref.is_empty() {
                let ext_ops = fast_gap_align(ext_read, ext_ref);
                for (op, len) in ext_ops {
                    add_cigar_op(&mut all_ops, op, len);
                }
            }
        }

        // Process each anchor and gap between them
        let mut prev_read_end = anchor_read_start;
        let mut prev_ref_end = anchor_ref_start;

        for anchor in anchors {
            let a_read_start = anchor.read_start as usize;
            let a_ref_start = anchor.ref_start as usize;
            let a_len = anchor.len as usize;

            // Align gap before this anchor
            if a_read_start > prev_read_end || a_ref_start > prev_ref_end {
                let gap_read = &read[prev_read_end..a_read_start];
                let gap_ref = &reference[prev_ref_end..a_ref_start];

                if !gap_read.is_empty() || !gap_ref.is_empty() {
                    let gap_ops = fast_gap_align(gap_read, gap_ref);
                    for (op, len) in gap_ops {
                        add_cigar_op(&mut all_ops, op, len);
                    }
                }
            }

            // Add anchor match (k-mer match, so it's exact)
            add_cigar_op(&mut all_ops, 'M', a_len);

            prev_read_end = a_read_start + a_len;
            prev_ref_end = a_ref_start + a_len;
        }

        // Right extension (after last anchor)
        if right_extend_read > 0 || right_extend_ref > 0 {
            let ext_read = &read[anchor_read_end..read_end];
            let ext_ref = &reference[anchor_ref_end..ref_end];
            if !ext_read.is_empty() || !ext_ref.is_empty() {
                let ext_ops = fast_gap_align(ext_read, ext_ref);
                for (op, len) in ext_ops {
                    add_cigar_op(&mut all_ops, op, len);
                }
            }
        }

        cigar_to_string(&all_ops)
    } else {
        // Simple CIGAR: just count total matched bases
        let total_match: usize = anchors.iter().map(|a| a.len as usize).sum();
        format!("{}M", total_match)
    };

    Some(AlignmentResult {
        ref_start,
        ref_end,
        read_start,
        read_end,
        cigar,
        chain_score,
        mapq,
        strand,
        _chrom_idx: 0, // Will be set by caller based on global position
    })
}

// ============================================================================
// TRAJECTORY-BASED ALIGNMENT (NO DP)
// Novel approach: CIGAR derived purely from geometric gap analysis
// ============================================================================

/// Micro-anchor boundary refinement - DP-free local refinement
///
/// When the first anchor is far from read start (common in noisy ONT reads),
/// we search for micro-anchors (smaller k-mers) in a window near the expected
/// start position to find a more accurate boundary.
///
/// Algorithm (O(W) where W = window size, typically 500-1000bp):
/// 1. Take first W bp of read
/// 2. Find all k-mer matches (k=8) in reference near expected diagonal
/// 3. Filter to matches consistent with global trajectory
/// 4. Use leftmost valid micro-anchor as refined start position
///
/// This is strictly DP-free: only hashing, filtering, and min operations.
/// Adaptive k-mer shrinking for boundary refinement
///
/// Key insight: ONT errors break long k-mers (k=15), but shorter k-mers (k=7-9)
/// can still find matches even in noisy regions. We try progressively smaller k
/// until we find consistent anchors near the read boundary.
///
/// Algorithm:
/// 1. Start with k=13, try to find anchors in first W bp of read
/// 2. If not enough consistent anchors → reduce to k=11
/// 3. If still not enough → k=9
/// 4. If still not enough → k=7
/// 5. Stop at first k that gives ≥3 consistent anchors
///
/// This is like mini-SW without any DP - just reducing seed stringency!
fn refine_boundary_with_microanchors(
    read: &[u8],
    reference: &[u8],
    first_anchor: &Anchor,
    slope: f64,
    window_size: usize,
    _micro_k: usize, // ignored, we use adaptive k
) -> Option<usize> {
    // Only refine if first anchor is far from read start
    if first_anchor.read_start < 50 {
        return None; // Already close enough
    }

    let window = window_size
        .min(read.len())
        .min(first_anchor.read_start as usize + 200);

    // Expected diagonal from first anchor
    let expected_diag = first_anchor.ref_start as i64 - first_anchor.read_start as i64;

    // Search range in reference (around expected start)
    let expected_ref_start = (expected_diag.max(0) as usize).saturating_sub(500);
    let expected_ref_end = ((expected_diag + window as i64 + 500) as usize).min(reference.len());

    // Trajectory estimate as fallback reference
    let trajectory_estimate =
        (first_anchor.ref_start as f64 - first_anchor.read_start as f64 * slope) as i64;

    // ADAPTIVE K-MER SHRINKING: try k=13, 11, 9, 7
    // Key: require MORE anchors for smaller k (to reduce false positives)
    let k_configs: [(usize, usize, i64); 4] = [
        (13, 3, 30), // (k, min_anchors, bandwidth)
        (11, 5, 25), // More anchors needed for smaller k
        (9, 8, 20),
        (7, 12, 15), // Very strict for k=7
    ];

    for (k, min_anchors, bandwidth) in k_configs {
        if window < k || expected_ref_end <= expected_ref_start + k {
            continue;
        }

        // Build k-mer index for reference window at this k
        let mut ref_kmers: FxHashMap<u64, Vec<u32>> = FxHashMap::default();
        for i in expected_ref_start..expected_ref_end.saturating_sub(k) {
            let kmer = &reference[i..i + k];
            if kmer
                .iter()
                .all(|&b| b == b'A' || b == b'C' || b == b'G' || b == b'T')
            {
                let hash = kmer_hash(kmer);
                ref_kmers.entry(hash).or_default().push(i as u32);
            }
        }

        // Find micro-anchors in read window
        let mut micro_anchors: Vec<(u32, u32)> = Vec::new();

        for read_pos in 0..window.saturating_sub(k) {
            let kmer = &read[read_pos..read_pos + k];
            if kmer
                .iter()
                .all(|&b| b == b'A' || b == b'C' || b == b'G' || b == b'T')
            {
                let hash = kmer_hash(kmer);
                if let Some(ref_positions) = ref_kmers.get(&hash) {
                    for &ref_pos in ref_positions {
                        let diag = ref_pos as i64 - read_pos as i64;
                        if (diag - expected_diag).abs() <= bandwidth {
                            micro_anchors.push((read_pos as u32, ref_pos));
                        }
                    }
                }
            }
        }

        // Need minimum anchors (more for smaller k)
        if micro_anchors.len() < min_anchors {
            continue; // Try smaller k
        }

        // Found enough anchors! Compute median estimate
        micro_anchors.sort_by_key(|(r, _)| *r);

        let mut estimates: Vec<i64> = micro_anchors
            .iter()
            .map(|(r, f)| (*f as i64) - ((*r as f64) * slope) as i64)
            .collect();
        estimates.sort();

        let median_ref_start = estimates[estimates.len() / 2];

        // Sanity check: stricter for smaller k
        let max_deviation = if k >= 11 { 150 } else { 100 };
        if (median_ref_start - trajectory_estimate).abs() > max_deviation {
            continue; // This k gave bad results, try smaller
        }

        // Success! Return the refined position
        return Some(median_ref_start.max(0) as usize);
    }

    // No k worked, fallback to trajectory extrapolation
    None
}

/// RANSAC (Random Sample Consensus) regression estimator
///
/// More robust than Theil-Sen: 50-60% breakdown point (can handle up to 50-60% outliers).
///
/// Algorithm:
/// 1. Randomly sample 2 points to define a line
/// 2. Count inliers (points within threshold distance)
/// 3. Repeat N iterations, keep best model
/// 4. Refit on all inliers of best model
///
/// For alignment: slope should be ~1.0, intercept = ref_start estimate
fn ransac_regression(points: &[(f64, f64)], iterations: usize, threshold: f64) -> (f64, f64) {
    if points.len() < 2 {
        return (1.0, points.first().map(|(x, y)| y - x).unwrap_or(0.0));
    }
    if points.len() == 2 {
        let dx = points[1].0 - points[0].0;
        if dx.abs() < 1e-10 {
            return (1.0, points[0].1 - points[0].0);
        }
        let slope = (points[1].1 - points[0].1) / dx;
        let intercept = points[0].1 - slope * points[0].0;
        return (slope, intercept);
    }

    // Simple pseudo-random using point coordinates as seed
    let mut best_inliers = 0;
    let mut best_slope = 1.0;
    let mut best_intercept = points[0].1 - points[0].0;

    let n = points.len();

    for iter in 0..iterations {
        // Pseudo-random selection based on iteration
        let i = iter % n;
        let j = (iter * 7 + 3) % n;
        if i == j {
            continue;
        }

        let (x1, y1) = points[i];
        let (x2, y2) = points[j];

        let dx = x2 - x1;
        if dx.abs() < 1e-10 {
            continue;
        }

        let slope = (y2 - y1) / dx;
        let intercept = y1 - slope * x1;

        // Skip unreasonable slopes (for alignment, slope should be ~1.0)
        if slope < 0.8 || slope > 1.2 {
            continue;
        }

        // Count inliers
        let mut inliers = 0;
        for &(x, y) in points {
            let predicted = slope * x + intercept;
            let error = (y - predicted).abs();
            if error <= threshold {
                inliers += 1;
            }
        }

        if inliers > best_inliers {
            best_inliers = inliers;
            best_slope = slope;
            best_intercept = intercept;
        }
    }

    // Refit on all inliers of best model
    if best_inliers >= 3 {
        let inlier_points: Vec<(f64, f64)> = points
            .iter()
            .filter(|&&(x, y)| {
                let predicted = best_slope * x + best_intercept;
                (y - predicted).abs() <= threshold
            })
            .copied()
            .collect();

        if inlier_points.len() >= 2 {
            // Simple least squares on inliers
            let n = inlier_points.len() as f64;
            let sum_x: f64 = inlier_points.iter().map(|(x, _)| x).sum();
            let sum_y: f64 = inlier_points.iter().map(|(_, y)| y).sum();
            let sum_xy: f64 = inlier_points.iter().map(|(x, y)| x * y).sum();
            let sum_x2: f64 = inlier_points.iter().map(|(x, _)| x * x).sum();

            let denom = n * sum_x2 - sum_x * sum_x;
            if denom.abs() > 1e-10 {
                best_slope = (n * sum_xy - sum_x * sum_y) / denom;
                best_intercept = (sum_y - best_slope * sum_x) / n;
            }
        }
    }

    (best_slope, best_intercept)
}

/// Theil-Sen robust regression estimator
///
/// Given a set of (x, y) points, estimates slope and intercept robustly.
/// - Slope: median of all pairwise slopes (y_j - y_i) / (x_j - x_i)
/// - Intercept: median of (y_i - slope * x_i)

// ============================================================================
// Trajectory quality scoring for chain selection
// ============================================================================

/// Calculate trajectory quality score for a chain
/// Returns a score 0.0-1.0 based on how well anchors fit a linear trajectory
fn calculate_trajectory_quality(anchors: &[Anchor]) -> f64 {
    if anchors.len() < 3 {
        return 0.0;
    }

    let first = anchors.first().unwrap();
    let last = anchors.last().unwrap();

    if last.read_start <= first.read_start {
        return 0.0;
    }

    // Calculate linear fit
    let slope = (last.ref_start as f64 - first.ref_start as f64)
        / (last.read_start as f64 - first.read_start as f64);
    let intercept = first.ref_start as f64 - slope * first.read_start as f64;

    // Count inliers with progressive thresholds
    let mut tight_inliers = 0; // Within 20bp
    let mut medium_inliers = 0; // Within 50bp
    let mut loose_inliers = 0; // Within 100bp

    for a in anchors {
        let expected_ref = slope * a.read_start as f64 + intercept;
        let residual = (a.ref_start as f64 - expected_ref).abs();
        if residual < 20.0 {
            tight_inliers += 1;
        }
        if residual < 50.0 {
            medium_inliers += 1;
        }
        if residual < 100.0 {
            loose_inliers += 1;
        }
    }

    let n = anchors.len() as f64;

    // Combined quality score: weighted inlier ratios
    let tight_ratio = tight_inliers as f64 / n;
    let medium_ratio = medium_inliers as f64 / n;
    let loose_ratio = loose_inliers as f64 / n;

    // Weight tight inliers more heavily
    0.5 * tight_ratio + 0.3 * medium_ratio + 0.2 * loose_ratio
}

/// TRAJECTORY-GUIDED CHAIN SELECTION
/// Use trajectory quality as TIE-BREAKER when chains have very similar scores
/// Only swap to alternative if: (1) scores within 5% AND (2) alternative has MUCH better trajectory
fn select_best_chain_by_trajectory(
    chain_result: &ChainResult,
    _read_len: usize,
) -> (Vec<Anchor>, i32, f64) {
    // Returns (best_chain, score, trajectory_quality)
    let best_traj_quality = calculate_trajectory_quality(&chain_result.best_chain);

    // Only consider alternatives within 5% of best score (tight tie-breaker)
    let score_threshold = (chain_result.best_score as f64 * 0.95) as i32;

    let mut best_alternative: Option<(&Vec<Anchor>, i32, f64)> = None;

    for (alt_chain, alt_score, _diag) in &chain_result.alternatives {
        if *alt_score >= score_threshold && !alt_chain.is_empty() {
            let alt_quality = calculate_trajectory_quality(alt_chain);

            // Only consider if trajectory quality is SIGNIFICANTLY better (>0.15 improvement)
            if alt_quality > best_traj_quality + 0.15 {
                match &best_alternative {
                    None => best_alternative = Some((alt_chain, *alt_score, alt_quality)),
                    Some((_, _, prev_quality)) => {
                        if alt_quality > *prev_quality {
                            best_alternative = Some((alt_chain, *alt_score, alt_quality));
                        }
                    }
                }
            }
        }
    }

    // Return alternative only if found one with significantly better trajectory
    match best_alternative {
        Some((chain, score, quality)) => (chain.clone(), score, quality),
        None => (
            chain_result.best_chain.clone(),
            chain_result.best_score,
            best_traj_quality,
        ),
    }
}

/// Trajectory-based CIGAR generation - NO DP whatsoever
///
/// Key insight: The ratio (ref_gap / read_gap) between consecutive anchors
/// directly tells us about indels:
/// - ratio = 1: match
/// - ratio > 1: deletion (ref advanced more)
/// - ratio < 1: insertion (read advanced more)
///
/// This is fundamentally different from DP-based approaches!
fn align_trajectory_based(
    anchors: &[Anchor],
    read: &[u8],
    reference: &[u8],
    chain_score: i32,
    mapq: u8,
    strand: char,
    refine_boundaries: bool,
    polish: bool,
    ransac_threshold: f64,
    polish_max_indel: usize,
) -> Option<AlignmentResult> {
    if anchors.len() < 2 {
        return None;
    }

    let first = anchors.first().unwrap();
    let last = anchors.last().unwrap();

    // =========================================================================
    // SEGMENTED RANSAC ROBUST REGRESSION for trajectory estimation
    //
    // Key insight for ONT reads: the first ~500-1000bp are often soft-clipped
    // or contain high error rates, making anchors in this region noisy/unreliable.
    //
    // Solution: fit RANSAC regression on anchors, then extrapolate back to read_pos=0.
    // RANSAC has 50-60% breakdown point (robust to ~50% outliers).
    //
    // Mathematically: intercept = ref_pos at read_pos=0 (extrapolated)
    // =========================================================================

    // Use anchors for RANSAC regression with segmented fit
    let n_for_regression = anchors.len().min(25);
    let points: Vec<(f64, f64)> = anchors[..n_for_regression]
        .iter()
        .map(|a| (a.read_start as f64, a.ref_start as f64))
        .collect();

    // Segmented fit: use first third and last third of anchors
    // This helps handle large indels in the middle of reads
    let (slope, intercept) = if points.len() >= 9 {
        let third = points.len() / 3;

        // Fit on first segment (determines ref_start)
        let first_seg = &points[..third.max(3)];
        let (s1, i1) = ransac_regression(first_seg, 30, ransac_threshold);

        // Fit on last segment
        let last_seg = &points[points.len() - third.max(3)..];
        let (s2, _i2) = ransac_regression(last_seg, 30, ransac_threshold);

        // Average slopes, use first segment's intercept
        let avg_slope = (s1 + s2) / 2.0;

        // If slopes differ significantly, use first segment only
        if (s1 - s2).abs() > 0.1 {
            (s1, i1)
        } else {
            // FIX: Use RANSAC intercept (i1) directly instead of min()
            // The min() caused systematic negative bias by always choosing leftward position
            (avg_slope, i1)
        }
    } else {
        ransac_regression(&points, 50, ransac_threshold * 1.2)
    };

    // =========================================================================
    // POSITION CALCULATION
    //
    // Use extrapolated position at read_pos=0 (intercept from regression).
    // This is a "synthetic extension" that estimates where the alignment would
    // start if we extended back to the beginning of the read.
    //
    // Rationale: minimap2 uses DP to extend alignment back to read_pos=0.
    // freemap uses trajectory extrapolation instead, which gives comparable
    // positions without the DP overhead. The soft-clip in CIGAR indicates
    // that these bases weren't verified by anchors.
    //
    // Validation (GIAB HG002 HiFi, 10k reads, Jan 2026):
    // - Mean position error vs minimap2: only 4.8 bp (excellent!)
    // - 93% concordance within 100bp vs minimap2's 94% self-concordance
    // - Remaining ~7% discordance is from repetitive regions, not extrapolation
    // - Offset analysis confirms current implementation is optimal (offset=0)
    // =========================================================================

    // POS = reference position at read_pos=0, from trajectory extrapolation (intercept).
    // The leading region [0, anchor_read_start) is ALIGNED (below), not soft-clipped, so
    // that POS is consistent with the CIGAR and read-start bases contribute to coverage /
    // variant calling — matching minimap2's DP-extension behaviour.
    //
    // BUGFIX: previously this region was soft-clipped while POS was set to the intercept.
    // A soft-clip consumes QUERY only (no reference), so the first 'M' base actually landed
    // ~anchor_read_start bases to the RIGHT of POS: every 'M' base was frame-shifted against
    // the reference and base-level identity collapsed to ~random (~28% on ONT), even though
    // coarse metrics (placement within 500bp, coverage in kb bins, indel-event counts) all
    // passed and hid the defect. Aligning the leading region restores ~98% identity.
    let mut ref_start = intercept.max(0.0) as usize;

    // Optionally refine boundary using micro-anchors (smaller k-mers)
    // This can help find matches closer to the extrapolated start
    if refine_boundaries {
        if let Some(refined) = refine_boundary_with_microanchors(
            read, reference, first, slope, 500, // window_size
            11,  // micro_k (ignored, function uses adaptive k)
        ) {
            // Use refined position if it's reasonable (between extrapolated and first anchor)
            let first_anchor_pos = first.ref_start as usize;
            if refined >= ref_start && refined <= first_anchor_pos {
                ref_start = refined;
            }
        }
    }

    // For the end, extrapolate from last anchor position rather than regression intercept.
    // This is more accurate when the chain spans an SV (large diagonal shift in the middle),
    // because the regression intercept reflects the first segment's diagonal.
    let last_anchor_ref_end = (last.ref_start as usize) + (last.len as usize);
    let right_read_remaining = read.len().saturating_sub(last.read_start as usize + last.len as usize);
    let ref_end = (last_anchor_ref_end + (right_read_remaining as f64 * slope) as usize)
        .min(reference.len());

    // Keep original anchor positions for CIGAR generation
    let anchor_read_start = first.read_start as usize;
    let anchor_read_end = (last.read_start + last.len as u32) as usize;
    let read_start = 0; // Always start from 0 for read
    let read_end = read.len();

    let mut all_ops: Vec<(char, usize)> = Vec::with_capacity(anchors.len() * 2);

    // Leading region (read[0..anchor_read_start]): ALIGN it against the reference window
    // [ref_start, first.ref_start) instead of soft-clipping, so the read-start bases
    // contribute to coverage / variant calling (like minimap2's extension). Guards:
    //  - cap the extended length (BOUNDARY_EXTEND_CAP): long unanchored prefixes are
    //    usually adapters / N-runs / junk that minimap2 also soft-clips (z-drop);
    //  - require read and reference spans to be comparable (slope ~ 1) — otherwise the
    //    extrapolation is unreliable and we soft-clip, anchoring POS at the first anchor.
    const BOUNDARY_EXTEND_CAP: usize = 1000;
    if anchor_read_start > 0 {
        let lead_ref_end = first.ref_start as usize;
        let lead_ref_len = lead_ref_end.saturating_sub(ref_start);
        let spans_ok = lead_ref_len > 0
            && lead_ref_end <= reference.len()
            && lead_ref_len <= anchor_read_start.saturating_mul(2) + 16
            && anchor_read_start <= lead_ref_len.saturating_mul(2) + 16;
        if anchor_read_start <= BOUNDARY_EXTEND_CAP && spans_ok {
            let lead_read = &read[0..anchor_read_start];
            let lead_ref = &reference[ref_start..lead_ref_end];
            let ops = if polish && anchor_read_start.abs_diff(lead_ref_len) <= polish_max_indel {
                polish_gap(lead_read, lead_ref, polish_max_indel)
            } else {
                fast_gap_align(lead_read, lead_ref)
            };
            for (op, len) in ops {
                add_cigar_op(&mut all_ops, op, len);
            }
        } else {
            // Unreliable / over-long prefix: soft-clip it and anchor POS at the first anchor.
            add_cigar_op(&mut all_ops, 'S', anchor_read_start);
            ref_start = lead_ref_end;
        }
    }

    // Process consecutive anchor pairs - THIS IS THE KEY TRAJECTORY ANALYSIS
    for i in 0..anchors.len() - 1 {
        let curr = &anchors[i];
        let next = &anchors[i + 1];

        let curr_read_end = curr.read_start as usize + curr.len as usize;
        let curr_ref_end = curr.ref_start as usize + curr.len as usize;
        let next_read_start = next.read_start as usize;
        let next_ref_start = next.ref_start as usize;

        // Gap between anchors
        let read_gap = next_read_start.saturating_sub(curr_read_end);
        let ref_gap = next_ref_start.saturating_sub(curr_ref_end);

        // Current anchor contributes matches
        if i == 0 {
            add_cigar_op(&mut all_ops, 'M', curr.len as usize);
        }

        // Analyze gap between anchors
        const MAX_INDEL_SIZE: usize = 50000; // 50kb max indel
        let indel_size = if ref_gap > read_gap {
            ref_gap - read_gap
        } else {
            read_gap - ref_gap
        };

        if indel_size > MAX_INDEL_SIZE {
            // Skip this anchor pair - gap too large to be real
            add_cigar_op(&mut all_ops, 'M', read_gap);
            add_cigar_op(&mut all_ops, 'M', next.len as usize);
        } else if read_gap == 0 && ref_gap == 0 {
            // Overlapping anchors - just add next anchor's contribution
            add_cigar_op(&mut all_ops, 'M', next.len as usize);
        } else if polish && (read_gap > 0 || ref_gap > 0)
            && curr_ref_end + ref_gap <= reference.len()
            && curr_read_end + read_gap <= read.len()
        {
            // POLISH MODE: use actual sequence comparison for gap regions
            let gap_read = &read[curr_read_end..curr_read_end + read_gap];
            let gap_ref = &reference[curr_ref_end..curr_ref_end + ref_gap];

            if gap_read.is_empty() && !gap_ref.is_empty() {
                add_cigar_op(&mut all_ops, 'D', ref_gap);
            } else if !gap_read.is_empty() && gap_ref.is_empty() {
                add_cigar_op(&mut all_ops, 'I', read_gap);
            } else if read_gap.abs_diff(ref_gap) <= polish_max_indel {
                // Small indel: diagonal_align_refined can detect it directly
                let gap_ops = polish_gap(gap_read, gap_ref, polish_max_indel);
                for (op, len) in gap_ops {
                    add_cigar_op(&mut all_ops, op, len);
                }
            } else {
                // Large indel: hybrid approach
                // 1. Polish the overlapping (equal-length) portion for base-level detail
                // 2. Emit the structural indel (excess) geometrically
                let overlap = std::cmp::min(read_gap, ref_gap);
                if overlap > 0 {
                    let gap_ops = polish_gap(&gap_read[..overlap], &gap_ref[..overlap], polish_max_indel);
                    for (op, len) in gap_ops {
                        add_cigar_op(&mut all_ops, op, len);
                    }
                }
                if read_gap > ref_gap {
                    add_cigar_op(&mut all_ops, 'I', read_gap - ref_gap);
                } else {
                    add_cigar_op(&mut all_ops, 'D', ref_gap - read_gap);
                }
            }
            add_cigar_op(&mut all_ops, 'M', next.len as usize);
        } else if read_gap == ref_gap {
            // Perfect diagonal - pure matches in gap (geometric estimate)
            add_cigar_op(&mut all_ops, 'M', read_gap);
            add_cigar_op(&mut all_ops, 'M', next.len as usize);
        } else if ref_gap > read_gap {
            // DELETION detected: reference advanced more than read
            add_cigar_op(&mut all_ops, 'M', read_gap);
            add_cigar_op(&mut all_ops, 'D', ref_gap - read_gap);
            add_cigar_op(&mut all_ops, 'M', next.len as usize);
        } else {
            // INSERTION detected: read advanced more than reference
            add_cigar_op(&mut all_ops, 'M', ref_gap);
            add_cigar_op(&mut all_ops, 'I', read_gap - ref_gap);
            add_cigar_op(&mut all_ops, 'M', next.len as usize);
        }
    }

    // Trailing region (read[anchor_read_end..]): ALIGN it against the reference window
    // [last_anchor_ref_end, ref_end) instead of soft-clipping (symmetric to the leading
    // region above; same cap and span-ratio guards).
    let remaining = read.len().saturating_sub(anchor_read_end);
    if remaining > 0 {
        let trail_ref_start = last_anchor_ref_end;
        let trail_ref_end = ref_end.min(reference.len());
        let trail_ref_len = trail_ref_end.saturating_sub(trail_ref_start);
        let spans_ok = trail_ref_len > 0
            && trail_ref_len <= remaining.saturating_mul(2) + 16
            && remaining <= trail_ref_len.saturating_mul(2) + 16;
        if remaining <= BOUNDARY_EXTEND_CAP && spans_ok {
            let trail_read = &read[anchor_read_end..read_end];
            let trail_ref = &reference[trail_ref_start..trail_ref_end];
            let ops = if polish && remaining.abs_diff(trail_ref_len) <= polish_max_indel {
                polish_gap(trail_read, trail_ref, polish_max_indel)
            } else {
                fast_gap_align(trail_read, trail_ref)
            };
            for (op, len) in ops {
                add_cigar_op(&mut all_ops, op, len);
            }
        } else {
            add_cigar_op(&mut all_ops, 'S', remaining);
        }
    }

    // Compact CIGAR
    let mut compact: Vec<(char, usize)> = Vec::new();
    for (op, len) in all_ops {
        if len == 0 {
            continue;
        }
        if let Some(last) = compact.last_mut() {
            if last.0 == op {
                last.1 += len;
                continue;
            }
        }
        compact.push((op, len));
    }

    let cigar = cigar_to_string(&compact);

    Some(AlignmentResult {
        ref_start,
        ref_end,
        read_start,
        read_end,
        cigar,
        chain_score,
        mapq,
        strand,
        _chrom_idx: 0, // Will be set by caller based on global position
    })
}

// ============================================================================
// Full alignment pipeline
// ============================================================================

fn align_read_single_strand<I: KmerIndex>(
    read: &[u8],
    reference: &[u8],
    index: &I,
    k: usize,
    w: usize,
    band_width: i32,
    strand: char,
    do_cigar: bool,
    trajectory_mode: bool,
    refine_boundaries: bool,
    max_lookback: usize,
    gap_max: i32,
    gap_scale: i32,
    ultralong: bool,
    max_occ: usize, // Maximum k-mer occurrence threshold for seeding
    global_ref: Option<&GlobalReference>, // For primary contig tie-breaker
    seeding_limits: &SeedingLimits,
    polish: bool, // Refine gap regions with sequence-level alignment (DP-free)
    ransac_threshold: f64,
    polish_max_indel: usize,
    no_tiebreaker: bool, // Disable primary-chromosome tiebreaker
) -> Option<(AlignmentResult, i32, Vec<(Vec<Anchor>, i32, i32)>)> {
    // Use provided max_occ, but override to 20 for ultralong mode
    // (ultralong needs more aggressive filtering due to repetitive genomes)
    let effective_max_occ = if ultralong { max_occ.min(20) } else { max_occ };
    let mut anchors = seed_read(read, index, k, w, effective_max_occ, seeding_limits);
    if anchors.is_empty() {
        return None;
    }

    // Filter out anchors on blacklisted regions (e.g., mitochondrial decoy GL000209.2)
    if let Some(gref) = global_ref {
        let blacklist = gref.get_blacklisted_ranges();
        if !blacklist.is_empty() {
            anchors = filter_blacklisted_anchors(anchors, &blacklist);
            if anchors.is_empty() {
                return None;
            }
        }
    }

    // Use top-K chaining for topological MAPQ calculation
    // Collect more alternatives (10) to better detect segmental duplications
    let mut chain_result = chain_anchors_topk(
        &mut anchors,
        band_width,
        10,
        max_lookback,
        gap_max,
        gap_scale,
        ultralong,
    );
    if chain_result.best_chain.is_empty() || chain_result.best_score < (k as i32) {
        return None;
    }

    // NOVEL: Topological MAPQ calculation
    // This correctly identifies tandem repeats and assigns low MAPQ
    // Calculate BEFORE tie-breaker to use as gating condition
    let mut mapq = calculate_topological_mapq(&chain_result, read.len());

    // Apply cross-chromosome penalty if alternatives exist on different chromosomes
    if let Some(gref) = global_ref {
        mapq = apply_cross_chromosome_penalty(mapq, &chain_result, gref, None);
    }

    // Apply primary contig tie-breaker if GlobalReference is available
    // AGGRESSIVE: Apply whenever best is on decoy (not gated by MAPQ)
    // This fixes chromosome errors where reads map to unplaced contigs (KI270xxx)
    // instead of the correct main chromosome with slightly lower score
    const TIEBREAKER_SCORE_EPSILON: f32 = 0.10; // 10% - allow primary with 90%+ of decoy score
    if !no_tiebreaker {
        if let Some(gref) = global_ref {
            let best_ref_pos = chain_result
                .best_chain
                .first()
                .map(|a| a.ref_start as usize)
                .unwrap_or(0);
            let best_is_primary = gref.is_primary_position(best_ref_pos);
            // Only apply tiebreaker when best is on decoy/alt contig
            if !best_is_primary {
                chain_result =
                    apply_primary_contig_tiebreaker(chain_result, gref, None, TIEBREAKER_SCORE_EPSILON);
            }
        }
    }

    // TRAJECTORY-GUIDED CHAIN SELECTION (key innovation!)
    // Instead of just using highest chain_score, pick the chain with best trajectory fit
    let (best_chain, score, _traj_quality) =
        select_best_chain_by_trajectory(&chain_result, read.len());
    let alternatives = chain_result.alternatives;

    // Choose alignment method
    let result = if trajectory_mode {
        // Pure trajectory (DP-free)
        align_trajectory_based(
            &best_chain,
            read,
            reference,
            score,
            mapq,
            strand,
            refine_boundaries,
            polish,
            ransac_threshold,
            polish_max_indel,
        )?
    } else {
        // Standard: heuristic gap alignment
        align_with_anchors_cigar(&best_chain, read, reference, score, mapq, strand, do_cigar)?
    };

    Some((result, score, alternatives))
}

fn align_read<I: KmerIndex>(
    read: &[u8],
    reference: &[u8],
    index: &I,
    k: usize,
    w: usize,
    band_width: i32,
    do_cigar: bool,
    trajectory_mode: bool,
    refine_boundaries: bool,
    short_read_mode: bool,
    max_lookback: usize,
    gap_max: i32,
    gap_scale: i32,
    ultralong: bool,
    max_occ: usize, // Maximum k-mer occurrence threshold for seeding
    global_ref: Option<&GlobalReference>, // For primary contig tie-breaker
    recovery_mode: bool,
    on_the_fly_hpc: bool,
    polish: bool,
    ransac_threshold: f64,
    polish_max_indel: usize,
    no_tiebreaker: bool, // Disable primary-chromosome tiebreaker
) -> Option<AlignmentResult> {
    let seeding_limits =
        select_seeding_limits(k, w, max_occ, ultralong, short_read_mode, recovery_mode);

    // Optional homopolymer compression on-the-fly (read side only)
    let (read_seq, read_len_effective) = if on_the_fly_hpc {
        let comp = compress_homopolymers(read).0;
        let len = comp.len();
        (comp, len)
    } else {
        (read.to_vec(), read.len())
    };

    // Try forward
    let fwd = align_read_single_strand(
        &read_seq,
        reference,
        index,
        k,
        w,
        band_width,
        '+',
        do_cigar,
        trajectory_mode,
        refine_boundaries,
        max_lookback,
        gap_max,
        gap_scale,
        ultralong,
        max_occ,
        global_ref,
        &seeding_limits,
        polish,
        ransac_threshold,
        polish_max_indel,
        no_tiebreaker,
    );

    // Short read optimization: skip reverse complement if forward alignment is strong
    if short_read_mode {
        if let Some((ref fwd_res, fwd_score, _)) = &fwd {
            let coverage = *fwd_score as f64 / read_len_effective as f64;
            if coverage >= 0.75 && fwd_res.mapq >= 30 {
                return fwd.map(|(mut res, _, _)| {
                    adjust_mapq_by_geometry(&mut res, read_len_effective);
                    res
                });
            }
        }
    }

    // Try reverse complement
    let rc_read = reverse_complement(&read_seq);
    let rev = align_read_single_strand(
        &rc_read,
        reference,
        index,
        k,
        w,
        band_width,
        '-',
        do_cigar,
        trajectory_mode,
        refine_boundaries,
        max_lookback,
        gap_max,
        gap_scale,
        ultralong,
        max_occ,
        global_ref,
        &seeding_limits,
        polish,
        ransac_threshold,
        polish_max_indel,
        no_tiebreaker,
    );

    // Select best alignment and apply identity-based MAPQ adjustment
    let mut result = match (fwd, rev) {
        (Some((fwd_res, fwd_score, _)), Some((rev_res, rev_score, _))) => {
            if fwd_score >= rev_score {
                Some(fwd_res)
            } else {
                Some(rev_res)
            }
        }
        (Some((res, _, _)), None) => Some(res),
        (None, Some((res, _, _))) => Some(res),
        (None, None) => None,
    };

    if let Some(ref mut res) = result {
        adjust_mapq_by_geometry(res, read_len_effective);
    }

    result
}

// ============================================================================
// MULTI-ALIGNMENT MODE (--multi)
// ============================================================================
//
// Output secondary (0x100) and supplementary (0x800) alignments for SV callers.

/// Find supplementary alignments by re-aligning large soft-clipped portions of the read.
/// This detects split reads spanning SV breakpoints.
fn find_supplementary_alignments<I: KmerIndex>(
    read: &[u8],
    primary: &AlignmentResult,
    reference: &[u8],
    index: &I,
    k: usize,
    w: usize,
    band_width: i32,
    do_cigar: bool,
    trajectory_mode: bool,
    refine_boundaries: bool,
    max_lookback: usize,
    gap_max: i32,
    gap_scale: i32,
    max_occ: usize,
    global_ref: Option<&GlobalReference>,
    seeding_limits: &SeedingLimits,
    polish: bool,
    ransac_threshold: f64,
    polish_max_indel: usize,
    no_tiebreaker: bool, // Disable primary-chromosome tiebreaker
) -> Vec<SupplementaryAlignment> {
    let min_clip = 500; // Minimum clip size to attempt supplementary alignment
    let mut supps = Vec::new();

    // Parse CIGAR to find clip sizes
    let cigar = &primary.cigar;
    let (left_clip, right_clip) = parse_clip_sizes(cigar);

    // Check left clip
    // In SAM, left clip of CIGAR corresponds to:
    // - Forward primary: start of original read (bases 0..left_clip)
    // - Reverse primary: end of original read (bases read_len-left_clip..read_len)
    if left_clip >= min_clip {
        let (clip_start, clip_len) = if primary.strand == '-' {
            (read.len() - left_clip, left_clip)
        } else {
            (0, left_clip)
        };
        let clip_seq = &read[clip_start..clip_start + clip_len];

        for &clip_strand in &['+', '-'] {
            let seq = if clip_strand == '-' {
                reverse_complement(clip_seq)
            } else {
                clip_seq.to_vec()
            };
            if let Some((mut result, _, _)) = align_read_single_strand(
                &seq, reference, index, k, w, band_width,
                clip_strand, do_cigar, trajectory_mode, refine_boundaries,
                max_lookback, gap_max, gap_scale, false, max_occ,
                global_ref, seeding_limits, polish,
                ransac_threshold, polish_max_indel, no_tiebreaker,
            ) {
                result.mapq = result.mapq.min(primary.mapq);
                supps.push(SupplementaryAlignment {
                    result,
                    clip_start,
                    clip_len,
                });
                break;
            }
        }
    }

    // Check right clip
    // Right clip of CIGAR corresponds to:
    // - Forward primary: end of original read (bases read_len-right_clip..read_len)
    // - Reverse primary: start of original read (bases 0..right_clip)
    if right_clip >= min_clip {
        let (clip_start, clip_len) = if primary.strand == '-' {
            (0, right_clip)
        } else {
            (read.len() - right_clip, right_clip)
        };
        let clip_seq = &read[clip_start..clip_start + clip_len];

        for &clip_strand in &['+', '-'] {
            let seq = if clip_strand == '-' {
                reverse_complement(clip_seq)
            } else {
                clip_seq.to_vec()
            };
            if let Some((mut result, _, _)) = align_read_single_strand(
                &seq, reference, index, k, w, band_width,
                clip_strand, do_cigar, trajectory_mode, refine_boundaries,
                max_lookback, gap_max, gap_scale, false, max_occ,
                global_ref, seeding_limits, polish,
                ransac_threshold, polish_max_indel, no_tiebreaker,
            ) {
                result.mapq = result.mapq.min(primary.mapq);
                supps.push(SupplementaryAlignment {
                    result,
                    clip_start,
                    clip_len,
                });
                break;
            }
        }
    }

    // Filter out supplementary alignments that overlap with the primary on the reference.
    // These are false supplementaries where the clip re-aligns to the same genomic region
    // as the primary, providing no useful split-read signal for SV callers.
    if let Some(gref) = global_ref {
        let (pri_chrom, _) = gref.global_to_local(primary.ref_start);
        supps.retain(|sup| {
            let (sup_chrom, _) = gref.global_to_local(sup.result.ref_start);
            if sup_chrom != pri_chrom {
                return true; // Different chromosome — keep (translocation evidence)
            }
            // Same chromosome: check for reference overlap
            let pri_start = primary.ref_start;
            let pri_end = primary.ref_end;
            let sup_start = sup.result.ref_start;
            let sup_end = sup.result.ref_end;
            // Compute overlap fraction relative to the supplementary
            let overlap_start = pri_start.max(sup_start);
            let overlap_end = pri_end.min(sup_end);
            if overlap_start < overlap_end {
                let overlap_len = overlap_end - overlap_start;
                let sup_len = sup_end.saturating_sub(sup_start).max(1);
                let overlap_frac = overlap_len as f64 / sup_len as f64;
                overlap_frac < 0.5 // Keep only if <50% overlap
            } else {
                true // No overlap — keep
            }
        });
    }

    supps
}

/// Parse leading and trailing soft-clip sizes from a CIGAR string.
fn parse_clip_sizes(cigar: &str) -> (usize, usize) {
    let mut left_clip = 0usize;
    let mut right_clip = 0usize;

    // Parse left clip
    let mut num = 0usize;
    let mut chars = cigar.chars().peekable();
    while let Some(&c) = chars.peek() {
        if c.is_ascii_digit() {
            num = num * 10 + (c as usize - '0' as usize);
            chars.next();
        } else {
            if c == 'S' || c == 'H' {
                left_clip = num;
            }
            break;
        }
    }

    // Parse right clip - find last operation
    let mut num = 0usize;
    let mut last_num = 0usize;
    let mut last_op = ' ';
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num = num * 10 + (c as usize - '0' as usize);
        } else {
            last_num = num;
            last_op = c;
            num = 0;
        }
    }
    if last_op == 'S' || last_op == 'H' {
        right_clip = last_num;
    }

    (left_clip, right_clip)
}

/// Align a read in multi-mode: produces primary + secondary + supplementary alignments.
fn align_read_multi<I: KmerIndex>(
    read: &[u8],
    reference: &[u8],
    index: &I,
    k: usize,
    w: usize,
    band_width: i32,
    do_cigar: bool,
    trajectory_mode: bool,
    refine_boundaries: bool,
    short_read_mode: bool,
    max_lookback: usize,
    gap_max: i32,
    gap_scale: i32,
    ultralong: bool,
    max_occ: usize,
    global_ref: Option<&GlobalReference>,
    recovery_mode: bool,
    on_the_fly_hpc: bool,
    polish: bool,
    ransac_threshold: f64,
    polish_max_indel: usize,
    max_secondary: usize,
    no_tiebreaker: bool, // Disable primary-chromosome tiebreaker
) -> Option<MultiAlignmentResult> {
    let seeding_limits =
        select_seeding_limits(k, w, max_occ, ultralong, short_read_mode, recovery_mode);

    let (read_seq, _read_len_effective) = if on_the_fly_hpc {
        let comp = compress_homopolymers(read).0;
        let len = comp.len();
        (comp, len)
    } else {
        (read.to_vec(), read.len())
    };

    // SV-aware chaining: use much wider diagonal band and larger lookback
    // to allow chaining across structural variants (deletions, insertions up to ~64kb).
    // With band_width=65536 (2^16), anchors whose diagonal differs by <65536 are
    // in the same band and can be chained. The gap penalty still correctly penalizes
    // diagonal shifts, but doesn't prevent cross-SV connections.
    let sv_band = 65536i32.max(band_width);
    let sv_lookback = 64usize.max(max_lookback);

    // Align both strands, keeping alternatives
    let fwd = align_read_single_strand(
        &read_seq, reference, index, k, w, sv_band, '+',
        do_cigar, trajectory_mode, refine_boundaries, sv_lookback,
        gap_max, gap_scale, ultralong, max_occ, global_ref,
        &seeding_limits, polish, ransac_threshold, polish_max_indel,
        no_tiebreaker,
    );

    let rc_read = reverse_complement(&read_seq);
    let rev = align_read_single_strand(
        &rc_read, reference, index, k, w, sv_band, '-',
        do_cigar, trajectory_mode, refine_boundaries, sv_lookback,
        gap_max, gap_scale, ultralong, max_occ, global_ref,
        &seeding_limits, polish, ransac_threshold, polish_max_indel,
        no_tiebreaker,
    );

    // Pick the best primary and get its alternatives
    let (mut primary, _score, winning_alts) = match (fwd, rev) {
        (Some((fwd_res, fwd_score, fwd_alts)), Some((rev_res, rev_score, rev_alts))) => {
            if fwd_score >= rev_score {
                (fwd_res, fwd_score, fwd_alts)
            } else {
                (rev_res, rev_score, rev_alts)
            }
        }
        (Some((res, score, alts)), None) => (res, score, alts),
        (None, Some((res, score, alts))) => (res, score, alts),
        (None, None) => return None,
    };

    adjust_mapq_by_geometry(&mut primary, read.len());

    // Build secondary alignments from the winning strand's alternatives
    let mut secondaries = Vec::new();
    let winning_strand = primary.strand;
    let winning_read = if winning_strand == '-' { &rc_read } else { &read_seq };

    for (alt_chain, alt_score, _diag) in winning_alts.iter().take(max_secondary) {
        if alt_chain.len() < 2 {
            continue;
        }
        // Skip alternatives with very low scores relative to primary
        if (*alt_score as f64) < (primary.chain_score as f64 * 0.5) {
            continue;
        }
        // Generate CIGAR for this alternative chain
        if let Some(mut result) = if trajectory_mode {
            align_trajectory_based(
                alt_chain, winning_read, reference, *alt_score,
                0, // MAPQ = 0 for secondaries
                winning_strand, refine_boundaries, polish, ransac_threshold, polish_max_indel,
            )
        } else {
            align_with_anchors_cigar(alt_chain, winning_read, reference, *alt_score, 0, winning_strand, do_cigar)
        } {
            result.mapq = 0; // Secondaries always MAPQ 0
            secondaries.push(result);
        }
    }

    // Find supplementary alignments from large soft-clips
    let supplementaries = find_supplementary_alignments(
        read, &primary, reference, index, k, w, band_width,
        do_cigar, trajectory_mode, refine_boundaries, max_lookback,
        gap_max, gap_scale, max_occ, global_ref,
        &seeding_limits, polish, ransac_threshold, polish_max_indel,
        no_tiebreaker,
    );

    Some(MultiAlignmentResult {
        primary,
        secondaries,
        supplementaries,
    })
}

// ============================================================================
// HOMOPOLYMER-COMPRESSED ALIGNMENT
// ============================================================================
//
// Align in homopolymer-compressed space, then map positions back to original.
// This eliminates ONT's systematic homopolymer indel errors.

/// Scan for exact k-mer match in original reference around estimated position.
/// Tries multiple k-mer sizes (15, 13, 11, 9) and multiple positions in the read.
/// Returns the refined ref_start position or the original estimate if no match found.
///
/// FIX: Now finds the match CLOSEST to estimated position, not the first (leftmost) match.
/// This eliminates systematic negative bias in position reporting.
fn refine_ref_start_kmer_scan(
    read: &[u8],
    original_ref: &[u8],
    estimated_ref_start: usize,
    _scan_kmer_len: usize, // ignored, we try multiple sizes
    scan_window: usize,
) -> usize {
    // Try progressively shorter k-mers (more tolerant to errors)
    let kmer_sizes = [15, 13, 11, 9];

    // Try k-mers at different positions in the read (in case start has errors)
    let read_offsets = [0, 10, 20, 50];

    for &kmer_len in &kmer_sizes {
        for &read_offset in &read_offsets {
            if read_offset + kmer_len > read.len() {
                continue;
            }

            let read_kmer = &read[read_offset..read_offset + kmer_len];

            // Adjust estimated ref position for read offset
            let est_ref_for_offset = estimated_ref_start + read_offset;

            // Define search window
            let window_start = est_ref_for_offset.saturating_sub(scan_window);
            let window_end =
                (est_ref_for_offset + scan_window).min(original_ref.len().saturating_sub(kmer_len));

            if window_start >= window_end {
                continue;
            }

            // Find the match CLOSEST to estimated position (not first/leftmost)
            let mut best_match: Option<usize> = None;
            let mut best_distance = usize::MAX;

            for pos in window_start..window_end {
                if pos + kmer_len > original_ref.len() {
                    break;
                }
                let ref_kmer = &original_ref[pos..pos + kmer_len];
                if ref_kmer == read_kmer {
                    let distance = if pos >= est_ref_for_offset {
                        pos - est_ref_for_offset
                    } else {
                        est_ref_for_offset - pos
                    };
                    if distance < best_distance {
                        best_distance = distance;
                        best_match = Some(pos);
                    }
                }
            }

            if let Some(pos) = best_match {
                // Found match! Calculate ref_start (accounting for read_offset)
                return pos.saturating_sub(read_offset);
            }
        }
    }

    // No match found, return original estimate
    estimated_ref_start
}

fn align_read_single_strand_compressed<I: KmerIndex>(
    read: &[u8],
    _compressed_ref: &[u8],
    ref_pos_map: &[u32],
    original_ref: &[u8],
    index: &I,
    k: usize,
    w: usize,
    band_width: i32,
    strand: char,
    original_ref_len: usize,
    max_lookback: usize,
    gap_max: i32,
    gap_scale: i32,
    ultralong: bool,
    max_occ: usize, // Maximum k-mer occurrence threshold for seeding
    global_ref: Option<&GlobalReference>, // For primary contig tie-breaker
    seeding_limits: &SeedingLimits,
    polish: bool,
    polish_max_indel: usize,
    no_tiebreaker: bool, // Disable primary-chromosome tiebreaker
) -> Option<(AlignmentResult, i32)> {
    // Compress the read
    let (compressed_read, read_pos_map) = compress_homopolymers(read);

    // Use provided max_occ, but cap at 20 for ultralong mode
    let effective_max_occ = if ultralong { max_occ.min(20) } else { max_occ };

    // Seed and chain in compressed space with top-K for MAPQ
    let mut anchors = seed_read(
        &compressed_read,
        index,
        k,
        w,
        effective_max_occ,
        seeding_limits,
    );
    if anchors.is_empty() {
        return None;
    }

    // Filter out anchors on blacklisted regions (convert compressed pos to original)
    if let Some(gref) = global_ref {
        let blacklist = gref.get_blacklisted_ranges();
        if !blacklist.is_empty() {
            anchors.retain(|a| {
                let compressed_pos = a.ref_start as usize;
                let original_pos = if compressed_pos < ref_pos_map.len() {
                    ref_pos_map[compressed_pos] as usize
                } else {
                    compressed_pos
                };
                !blacklist
                    .iter()
                    .any(|(start, end)| original_pos >= *start && original_pos < *end)
            });
            if anchors.is_empty() {
                return None;
            }
        }
    }

    // Collect more alternatives (10) to better detect segmental duplications
    let mut chain_result = chain_anchors_topk(
        &mut anchors,
        band_width,
        10,
        max_lookback,
        gap_max,
        gap_scale,
        ultralong,
    );
    if chain_result.best_chain.is_empty() || chain_result.best_score < (k as i32) {
        return None;
    }

    // NOVEL: Topological MAPQ calculation
    // Calculate BEFORE tie-breaker to use as gating condition
    let mut mapq = calculate_topological_mapq(&chain_result, compressed_read.len());

    // Apply cross-chromosome penalty if alternatives exist on different chromosomes
    // Compressed mode: pass ref_pos_map to convert compressed positions to original
    if let Some(gref) = global_ref {
        mapq = apply_cross_chromosome_penalty(mapq, &chain_result, gref, Some(ref_pos_map));
    }

    // Apply primary contig tie-breaker if GlobalReference is available
    // AGGRESSIVE: Apply whenever best is on decoy (not gated by MAPQ)
    // This fixes chromosome errors where reads map to unplaced contigs (KI270xxx)
    const TIEBREAKER_SCORE_EPSILON: f32 = 0.10; // 10% - allow primary with 90%+ of decoy score
    if !no_tiebreaker {
        if let Some(gref) = global_ref {
            let best_compressed_pos = chain_result
                .best_chain
                .first()
                .map(|a| a.ref_start as usize)
                .unwrap_or(0);
            let best_original_pos = if best_compressed_pos < ref_pos_map.len() {
                ref_pos_map[best_compressed_pos] as usize
            } else {
                best_compressed_pos
            };
            let best_is_primary = gref.is_primary_position(best_original_pos);
            // Only apply tiebreaker when best is on decoy/alt contig
            if !best_is_primary {
                chain_result = apply_primary_contig_tiebreaker(
                    chain_result,
                    gref,
                    Some(ref_pos_map),
                    TIEBREAKER_SCORE_EPSILON,
                );
            }
        }
    }

    // TRAJECTORY-GUIDED CHAIN SELECTION (key innovation!)
    // Instead of just using highest chain_score, pick the chain with best trajectory fit
    let (best_chain, score, _traj_quality) =
        select_best_chain_by_trajectory(&chain_result, read.len());
    let _alternatives = chain_result.alternatives;

    // =========================================================================
    // TRAJECTORY-BASED CIGAR GENERATION FOR COMPRESSED ALIGNMENT
    // Convert anchors to original coordinates and analyze gaps
    // =========================================================================

    let (cigar, cigar_ref_consumed, cigar_ref_start) = if best_chain.len() >= 2 {
        // Convert all anchors to original coordinates
        let orig_anchors: Vec<Anchor> = best_chain
            .iter()
            .map(|a| {
                let orig_read_start = decompress_pos(a.read_start as usize, &read_pos_map);
                let orig_ref_start = decompress_pos(a.ref_start as usize, ref_pos_map);
                // Anchor length in original space: estimate from position mapping
                let orig_len = if (a.read_start as usize + a.len as usize) < read_pos_map.len() {
                    let end_pos =
                        decompress_pos(a.read_start as usize + a.len as usize, &read_pos_map);
                    (end_pos - orig_read_start).max(a.len as usize)
                } else {
                    a.len as usize
                };
                Anchor {
                    read_start: orig_read_start as u32,
                    ref_start: orig_ref_start as u32,
                    len: orig_len as u16,
                }
            })
            .collect();

        // Generate the CIGAR by BASE-ALIGNING in ORIGINAL coordinates, segment by segment
        // between consecutive decompressed anchor start positions.
        //
        // The anchors are exact matches in COMPRESSED space, but in original space read and
        // reference homopolymer run lengths differ, so emitting anchor / equal-gap regions as
        // plain 'M' (as the previous version did) left ~18% residual base error at homopolymer
        // boundaries. Polishing each segment with diagonal_align_refined resolves those bases,
        // while the dense HPC anchors act as re-sync points -> accurate original-space alignment.
        // POS is the first decompressed anchor's reference position (leading region soft-clipped).
        let mut all_ops: Vec<(char, usize)> = Vec::with_capacity(orig_anchors.len() * 2);
        let first = &orig_anchors[0];
        let last = &orig_anchors[orig_anchors.len() - 1];
        let anchor_read_start = first.read_start as usize;
        let anchor_read_end = last.read_start as usize + last.len as usize;
        let mut total_ref_consumed: usize = 0;

        // Leading region (before first anchor) - soft-clip (not verified by anchors)
        if anchor_read_start > 0 {
            add_cigar_op(&mut all_ops, 'S', anchor_read_start);
        }

        // Segment boundaries: each anchor's start position, plus the final aligned end.
        let mut read_bounds: Vec<usize> =
            orig_anchors.iter().map(|a| a.read_start as usize).collect();
        read_bounds.push(anchor_read_end);
        let mut ref_bounds: Vec<usize> =
            orig_anchors.iter().map(|a| a.ref_start as usize).collect();
        ref_bounds.push(last.ref_start as usize + last.len as usize);

        for s in 0..read_bounds.len() - 1 {
            let r0 = read_bounds[s];
            let r1 = read_bounds[s + 1].max(r0);
            let f0 = ref_bounds[s];
            let f1 = ref_bounds[s + 1].max(f0);
            if r1 > read.len() || f1 > original_ref.len() {
                break; // safety: implausible decompressed coordinates
            }
            let rseg = &read[r0..r1];
            let fseg = &original_ref[f0..f1];
            let seg_ops: Vec<(char, usize)> = if rseg.is_empty() && fseg.is_empty() {
                Vec::new()
            } else if rseg.is_empty() {
                vec![('D', fseg.len())]
            } else if fseg.is_empty() {
                vec![('I', rseg.len())]
            } else if polish {
                polish_gap(rseg, fseg, polish_max_indel)
            } else {
                fast_gap_align(rseg, fseg)
            };
            for (op, len) in seg_ops {
                add_cigar_op(&mut all_ops, op, len);
                if op == 'M' || op == 'D' || op == 'N' || op == '=' || op == 'X' {
                    total_ref_consumed += len;
                }
            }
        }

        // Trailing region (after last anchor) - soft-clip (not verified by anchors)
        let remaining = read.len().saturating_sub(anchor_read_end);
        if remaining > 0 {
            add_cigar_op(&mut all_ops, 'S', remaining);
        }

        // Compact CIGAR
        let mut compact: Vec<(char, usize)> = Vec::new();
        for (op, len) in all_ops {
            if len == 0 {
                continue;
            }
            if let Some(last) = compact.last_mut() {
                if last.0 == op {
                    last.1 += len;
                    continue;
                }
            }
            compact.push((op, len));
        }

        (
            cigar_to_string(&compact),
            total_ref_consumed,
            orig_anchors[0].ref_start as usize,
        )
    } else {
        // Single anchor - simple CIGAR
        // Decompress anchor position for ref_start calculation
        let first_anchor = &best_chain[0];
        let anchor_ref_start = decompress_pos(first_anchor.ref_start as usize, ref_pos_map);
        let anchor_read_start = decompress_pos(first_anchor.read_start as usize, &read_pos_map);
        let estimated_ref_start = anchor_ref_start.saturating_sub(anchor_read_start);

        // K-MER SCAN: Refine ref_start with larger window
        let ref_start_calc =
            refine_ref_start_kmer_scan(read, original_ref, estimated_ref_start, 15, 2000);

        (format!("{}M", read.len()), read.len(), ref_start_calc)
    };

    // Calculate ref_end from CIGAR-consumed reference bases
    let actual_ref_start = cigar_ref_start;
    let actual_ref_end = (actual_ref_start + cigar_ref_consumed).min(original_ref_len);

    Some((
        AlignmentResult {
            ref_start: actual_ref_start,
            ref_end: actual_ref_end,
            read_start: 0,
            read_end: read.len(),
            cigar,
            chain_score: score,
            mapq,
            strand,
            _chrom_idx: 0, // Will be set by caller based on global position
        },
        score,
    ))
}

fn align_read_compressed<I: KmerIndex>(
    read: &[u8],
    compressed_ref: &[u8],
    ref_pos_map: &[u32],
    original_ref: &[u8],
    index: &I,
    k: usize,
    w: usize,
    band_width: i32,
    original_ref_len: usize,
    max_lookback: usize,
    gap_max: i32,
    gap_scale: i32,
    ultralong: bool,
    max_occ: usize, // Maximum k-mer occurrence threshold for seeding
    global_ref: Option<&GlobalReference>, // For primary contig tie-breaker
    recovery_mode: bool,
    _on_the_fly_hpc: bool,
    polish: bool,
    polish_max_indel: usize,
    no_tiebreaker: bool, // Disable primary-chromosome tiebreaker
) -> Option<AlignmentResult> {
    let seeding_limits = select_seeding_limits(k, w, max_occ, ultralong, false, recovery_mode);
    // Try forward
    let fwd = align_read_single_strand_compressed(
        read,
        compressed_ref,
        ref_pos_map,
        original_ref,
        index,
        k,
        w,
        band_width,
        '+',
        original_ref_len,
        max_lookback,
        gap_max,
        gap_scale,
        ultralong,
        max_occ,
        global_ref,
        &seeding_limits,
        polish,
        polish_max_indel,
        no_tiebreaker,
    );

    // Try reverse complement
    let rc_read = reverse_complement(read);
    let rev = align_read_single_strand_compressed(
        &rc_read,
        compressed_ref,
        ref_pos_map,
        original_ref,
        index,
        k,
        w,
        band_width,
        '-',
        original_ref_len,
        max_lookback,
        gap_max,
        gap_scale,
        ultralong,
        max_occ,
        global_ref,
        &seeding_limits,
        polish,
        polish_max_indel,
        no_tiebreaker,
    );

    // Select best alignment and apply identity-based MAPQ adjustment
    let mut result = match (fwd, rev) {
        (Some((fwd_res, fwd_score)), Some((rev_res, rev_score))) => {
            if fwd_score >= rev_score {
                Some(fwd_res)
            } else {
                Some(rev_res)
            }
        }
        (Some((res, _)), None) => Some(res),
        (None, Some((res, _))) => Some(res),
        (None, None) => None,
    };

    // Post-hoc MAPQ adjustment based on alignment identity
    if let Some(ref mut res) = result {
        adjust_mapq_by_geometry(res, read.len());
    }

    result
}

// ============================================================================
// PAF output
// ============================================================================

/// Calculate number of matching bases from CIGAR string
/// Counts M (match/mismatch) and = (sequence match) operations
fn calculate_nmatch_from_cigar(cigar: &str) -> usize {
    let mut nmatch = 0;
    let mut num = 0;
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num = num * 10 + (c as usize - '0' as usize);
        } else {
            if c == 'M' || c == '=' {
                nmatch += num;
            }
            num = 0;
        }
    }
    nmatch
}

/// Calculate geometric signals from CIGAR for MAPQ calibration
/// Returns (indel_burden, max_gap_ratio, aligned_len)
/// - indel_burden: (insertions + deletions) / aligned_length
/// - max_gap_ratio: largest single indel / read_length
fn calculate_cigar_geometric_signals(cigar: &str, read_len: usize) -> (f64, f64, usize) {
    let mut aligned_len = 0usize;
    let mut total_indels = 0usize;
    let mut max_indel = 0usize;
    let mut num = 0usize;

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num = num * 10 + (c as usize - '0' as usize);
        } else {
            match c {
                'M' | '=' | 'X' => {
                    aligned_len += num;
                }
                'I' | 'D' => {
                    aligned_len += num;
                    total_indels += num;
                    max_indel = max_indel.max(num);
                }
                _ => {} // S, H, N, P - don't count
            }
            num = 0;
        }
    }

    let indel_burden = if aligned_len > 0 {
        total_indels as f64 / aligned_len as f64
    } else {
        0.0
    };

    let max_gap_ratio = if read_len > 0 {
        max_indel as f64 / read_len as f64
    } else {
        0.0
    };

    (indel_burden, max_gap_ratio, aligned_len)
}

/// Adjust MAPQ based on geometric signals derived from trajectory alignment
///
/// This is a post-hoc geometric calibration that caps MAPQ when alignment quality
/// indicators suggest unreliable mapping. Unlike DP-based MAPQ, this uses only
/// signals derivable from the trajectory representation:
///
/// 1. Indel burden: high (I+D)/aligned_len indicates messy alignment
/// 2. Max gap ratio: large single indel relative to read suggests structural issue
/// 3. Aligned length ratio: partial alignments are less reliable
///
/// This approach is consistent with our DP-free philosophy while significantly
/// reducing wrong-but-confident mappings.
fn adjust_mapq_by_geometry(result: &mut AlignmentResult, read_len: usize) {
    if result.cigar.is_empty() || read_len == 0 {
        return;
    }

    let (indel_burden, max_gap_ratio, aligned_len) =
        calculate_cigar_geometric_signals(&result.cigar, read_len);

    // Aligned length ratio (how much of the read is aligned)
    let aligned_ratio = aligned_len as f64 / read_len as f64;

    let is_short = read_len < 300;

    // Start with maximum MAPQ cap
    let mut max_mapq: u8 = 60;

    // Rule 1: High indel burden indicates messy alignment
    // Short reads (Illumina ~0.1% error): much tighter thresholds
    // Long reads (HiFi 1-5%, ONT 5-15%): more permissive
    if is_short {
        if indel_burden > 0.10 {
            max_mapq = max_mapq.min(3);
        } else if indel_burden > 0.05 {
            max_mapq = max_mapq.min(10);
        } else if indel_burden > 0.02 {
            max_mapq = max_mapq.min(20);
        }
    } else {
        if indel_burden > 0.25 {
            max_mapq = max_mapq.min(3);
        } else if indel_burden > 0.20 {
            max_mapq = max_mapq.min(10);
        } else if indel_burden > 0.15 {
            max_mapq = max_mapq.min(20);
        }
    }

    // Rule 2: Large single gap suggests structural problem
    if is_short {
        if max_gap_ratio > 0.15 {
            max_mapq = max_mapq.min(3);
        } else if max_gap_ratio > 0.08 {
            max_mapq = max_mapq.min(10);
        } else if max_gap_ratio > 0.04 {
            max_mapq = max_mapq.min(20);
        }
    } else {
        if max_gap_ratio > 0.30 {
            max_mapq = max_mapq.min(3);
        } else if max_gap_ratio > 0.20 {
            max_mapq = max_mapq.min(10);
        } else if max_gap_ratio > 0.10 {
            max_mapq = max_mapq.min(20);
        }
    }

    // Rule 3: Partial alignments are less reliable
    if is_short {
        if aligned_ratio < 0.70 {
            max_mapq = max_mapq.min(3);
        } else if aligned_ratio < 0.85 {
            max_mapq = max_mapq.min(10);
        } else if aligned_ratio < 0.92 {
            max_mapq = max_mapq.min(20);
        }
    } else {
        if aligned_ratio < 0.50 {
            max_mapq = max_mapq.min(3);
        } else if aligned_ratio < 0.70 {
            max_mapq = max_mapq.min(10);
        } else if aligned_ratio < 0.85 {
            max_mapq = max_mapq.min(20);
        }
    }

    result.mapq = result.mapq.min(max_mapq);
}

fn format_paf_line(
    read_name: &str,
    read_len: usize,
    ref_name: &str,
    ref_len: usize,
    result: &AlignmentResult,
) -> String {
    let alen = result.ref_end - result.ref_start;
    let nmatch = calculate_nmatch_from_cigar(&result.cigar);
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tcg:Z:{}",
        read_name,
        read_len,
        result.read_start,
        result.read_end,
        result.strand,
        ref_name,
        ref_len,
        result.ref_start,
        result.ref_end,
        nmatch,
        alen,
        result.mapq,
        result.cigar
    )
}

fn format_sam_line(
    read_name: &str,
    read_seq: &[u8],
    qual: &[u8],
    ref_name: &str,
    result: &AlignmentResult,
) -> String {
    // Pre-allocate: name + seq + qual + cigar + fixed fields
    let estimated = read_name.len() + read_seq.len() * 2 + result.cigar.len() + ref_name.len() + 64;
    let mut buf = String::with_capacity(estimated);
    let mut ibuf = itoa::Buffer::new();

    // Read name
    buf.push_str(read_name);
    buf.push('\t');

    // FLAG
    let flag: u16 = if result.strand == '-' { 16 } else { 0 };
    buf.push_str(ibuf.format(flag));
    buf.push('\t');

    // RNAME
    buf.push_str(ref_name);
    buf.push('\t');

    // POS (1-based)
    buf.push_str(ibuf.format(result.ref_start + 1));
    buf.push('\t');

    // MAPQ
    buf.push_str(ibuf.format(result.mapq));
    buf.push('\t');

    // CIGAR
    if result.cigar.is_empty() {
        buf.push_str(ibuf.format(read_seq.len()));
        buf.push('M');
    } else {
        buf.push_str(&result.cigar);
    }

    // RNEXT, PNEXT, TLEN
    buf.push_str("\t*\t0\t0\t");

    // SEQ - write directly, no intermediate allocation
    if result.strand == '-' {
        for &b in read_seq.iter().rev() {
            buf.push(complement(b) as char);
        }
    } else {
        // SAFETY: FASTQ bases are ASCII, valid UTF-8
        let vec = unsafe { buf.as_mut_vec() };
        vec.extend_from_slice(read_seq);
    }
    buf.push('\t');

    // QUAL - write directly, no intermediate allocation
    if qual.is_empty() {
        buf.push('*');
    } else if result.strand == '-' {
        let vec = unsafe { buf.as_mut_vec() };
        vec.extend(qual.iter().rev());
    } else {
        let vec = unsafe { buf.as_mut_vec() };
        vec.extend_from_slice(qual);
    }

    buf
}

/// Generate SA:Z tag value for a single alignment entry.
/// Format: rname,pos,strand,CIGAR,mapQ,NM;
fn format_sa_entry(ref_name: &str, result: &AlignmentResult) -> String {
    format!(
        "{},{},{},{},{},0",
        ref_name,
        result.ref_start + 1,
        result.strand,
        if result.cigar.is_empty() { "*" } else { &result.cigar },
        result.mapq,
    )
}

/// Format a SAM line for multi-mode output (secondary/supplementary alignments).
/// `extra_flags`: 0 for primary, 0x100 for secondary, 0x800 for supplementary.
/// `sa_tag`: Optional SA:Z tag string (for primary and supplementary records).
fn format_sam_line_multi(
    read_name: &str,
    read_seq: &[u8],
    qual: &[u8],
    ref_name: &str,
    result: &AlignmentResult,
    extra_flags: u16,
    sa_tag: Option<&str>,
) -> String {
    let estimated = read_name.len() + read_seq.len() * 2 + result.cigar.len() + ref_name.len() + 128;
    let mut buf = String::with_capacity(estimated);
    let mut ibuf = itoa::Buffer::new();

    // Read name
    buf.push_str(read_name);
    buf.push('\t');

    // FLAG
    let mut flag: u16 = extra_flags;
    if result.strand == '-' {
        flag |= 16;
    }
    // Secondary reads should not be marked as primary
    if extra_flags & 0x100 != 0 {
        flag |= 0x100; // already set
    }
    buf.push_str(ibuf.format(flag));
    buf.push('\t');

    // RNAME
    buf.push_str(ref_name);
    buf.push('\t');

    // POS (1-based)
    buf.push_str(ibuf.format(result.ref_start + 1));
    buf.push('\t');

    // MAPQ
    buf.push_str(ibuf.format(result.mapq));
    buf.push('\t');

    // CIGAR - keep soft clips for all record types (full SEQ is always emitted)
    if result.cigar.is_empty() {
        buf.push_str(ibuf.format(read_seq.len()));
        buf.push('M');
    } else {
        buf.push_str(&result.cigar);
    }

    // RNEXT, PNEXT, TLEN
    buf.push_str("\t*\t0\t0\t");

    // SEQ - output full sequence for all records (sniffles2 and other SV callers need it)
    if result.strand == '-' {
        for &b in read_seq.iter().rev() {
            buf.push(complement(b) as char);
        }
    } else {
        let vec = unsafe { buf.as_mut_vec() };
        vec.extend_from_slice(read_seq);
    }
    buf.push('\t');

    // QUAL
    if qual.is_empty() {
        buf.push('*');
    } else if result.strand == '-' {
        let vec = unsafe { buf.as_mut_vec() };
        vec.extend(qual.iter().rev());
    } else {
        let vec = unsafe { buf.as_mut_vec() };
        vec.extend_from_slice(qual);
    }

    // SA tag
    if let Some(sa) = sa_tag {
        buf.push_str("\tSA:Z:");
        buf.push_str(sa);
    }

    buf
}

/// Format a SAM line for a supplementary alignment (0x800).
/// The clip portion of the original read is output as SEQ, with hard clips
/// for the rest of the read.
fn format_sam_line_supplementary(
    read_name: &str,
    full_read_seq: &[u8],
    full_qual: &[u8],
    ref_name: &str,
    result: &AlignmentResult,
    clip_start: usize,
    clip_len: usize,
    sa_tag: Option<&str>,
) -> String {
    let read_len = full_read_seq.len();
    let clip_seq = &full_read_seq[clip_start..clip_start + clip_len];
    let clip_qual = if full_qual.is_empty() {
        &[] as &[u8]
    } else {
        &full_qual[clip_start..clip_start + clip_len]
    };

    // Build CIGAR: prepend/append hard clips for non-clip portions
    let hard_left = clip_start;
    let hard_right = read_len - (clip_start + clip_len);

    let mut cigar_buf = String::with_capacity(result.cigar.len() + 20);
    let mut cibuf = itoa::Buffer::new();
    if hard_left > 0 {
        cigar_buf.push_str(cibuf.format(hard_left));
        cigar_buf.push('H');
    }
    if result.cigar.is_empty() {
        cigar_buf.push_str(cibuf.format(clip_len));
        cigar_buf.push('M');
    } else {
        cigar_buf.push_str(&result.cigar);
    }
    if hard_right > 0 {
        cigar_buf.push_str(cibuf.format(hard_right));
        cigar_buf.push('H');
    }

    let estimated = read_name.len() + clip_len * 2 + cigar_buf.len() + ref_name.len() + 128;
    let mut buf = String::with_capacity(estimated);
    let mut ibuf = itoa::Buffer::new();

    // Read name
    buf.push_str(read_name);
    buf.push('\t');

    // FLAG: 0x800 (supplementary) + 0x10 if reverse
    let flag: u16 = 0x800 | if result.strand == '-' { 16 } else { 0 };
    buf.push_str(ibuf.format(flag));
    buf.push('\t');

    // RNAME
    buf.push_str(ref_name);
    buf.push('\t');

    // POS (1-based)
    buf.push_str(ibuf.format(result.ref_start + 1));
    buf.push('\t');

    // MAPQ
    buf.push_str(ibuf.format(result.mapq));
    buf.push('\t');

    // CIGAR
    buf.push_str(&cigar_buf);

    // RNEXT, PNEXT, TLEN
    buf.push_str("\t*\t0\t0\t");

    // SEQ - only the clip portion
    if result.strand == '-' {
        for &b in clip_seq.iter().rev() {
            buf.push(complement(b) as char);
        }
    } else {
        let vec = unsafe { buf.as_mut_vec() };
        vec.extend_from_slice(clip_seq);
    }
    buf.push('\t');

    // QUAL - only the clip portion
    if clip_qual.is_empty() {
        buf.push('*');
    } else if result.strand == '-' {
        let vec = unsafe { buf.as_mut_vec() };
        vec.extend(clip_qual.iter().rev());
    } else {
        let vec = unsafe { buf.as_mut_vec() };
        vec.extend_from_slice(clip_qual);
    }

    // SA tag
    if let Some(sa) = sa_tag {
        buf.push_str("\tSA:Z:");
        buf.push_str(sa);
    }

    buf
}

/// Format a paired-end SAM line with proper FLAGS, RNEXT, PNEXT, TLEN
fn format_sam_line_paired(
    read_name: &str,
    read_seq: &[u8],
    qual: &[u8],
    ref_name: &str,
    result: &AlignmentResult,
    mate_ref_name: Option<&str>,
    mate_pos: Option<usize>,
    is_first_in_pair: bool, // true for R1, false for R2
    is_proper_pair: bool,
    tlen: i64,
) -> String {
    // Build SAM FLAG
    let mut flag: u16 = 0x1; // paired

    // Is this a proper pair?
    if is_proper_pair {
        flag |= 0x2; // proper pair
    }

    // Mate unmapped?
    if mate_pos.is_none() {
        flag |= 0x8; // mate unmapped
    }

    // Read strand
    if result.strand == '-' {
        flag |= 0x10; // read reverse strand
    }

    // Mate strand (for proper pair: R1 forward, R2 reverse)
    // If proper pair and we are R1, mate (R2) should be on reverse strand
    // If proper pair and we are R2, mate (R1) should be on forward strand
    if is_proper_pair {
        if is_first_in_pair {
            flag |= 0x20; // mate on reverse strand (R2 is reverse)
        }
        // If R2, mate is forward, so no flag needed
    }

    // First or second in pair
    if is_first_in_pair {
        flag |= 0x40; // first in pair
    } else {
        flag |= 0x80; // second in pair
    }

    // Position (1-based)
    let pos = result.ref_start + 1;

    // Sequence (reverse complement if on reverse strand)
    let seq_str = if result.strand == '-' {
        reverse_complement(read_seq)
            .into_iter()
            .map(|b| b as char)
            .collect::<String>()
    } else {
        read_seq.iter().map(|&b| b as char).collect::<String>()
    };

    // CIGAR
    let cigar = if result.cigar.is_empty() {
        format!("{}M", read_seq.len())
    } else {
        result.cigar.clone()
    };

    // RNEXT (mate reference name)
    let rnext = match mate_ref_name {
        Some(name) if name == ref_name => "=".to_string(),
        Some(name) => name.to_string(),
        None => "*".to_string(),
    };

    // PNEXT (mate position, 1-based)
    let pnext = mate_pos.map(|p| p + 1).unwrap_or(0);

    // TLEN (template length)
    let tlen_str = if is_proper_pair { tlen } else { 0 };

    // Quality string: reverse for reverse strand, * if empty (FASTA input)
    let qual_str = if qual.is_empty() {
        "*".to_string()
    } else if result.strand == '-' {
        qual.iter().rev().map(|&b| b as char).collect::<String>()
    } else {
        qual.iter().map(|&b| b as char).collect::<String>()
    };

    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        read_name, flag, ref_name, pos, result.mapq, cigar, rnext, pnext, tlen_str, seq_str, qual_str
    )
}

/// Format an unmapped paired-end SAM line
fn format_sam_line_unmapped_paired(
    read_name: &str,
    read_seq: &[u8],
    qual: &[u8],
    mate_ref_name: Option<&str>,
    mate_pos: Option<usize>,
    is_first_in_pair: bool,
) -> String {
    let mut flag: u16 = 0x1 | 0x4; // paired + unmapped

    if mate_pos.is_none() {
        flag |= 0x8; // mate also unmapped
    }

    if is_first_in_pair {
        flag |= 0x40;
    } else {
        flag |= 0x80;
    }

    let rnext = mate_ref_name.unwrap_or("*");
    let pnext = mate_pos.map(|p| p + 1).unwrap_or(0);
    let seq_str: String = read_seq.iter().map(|&b| b as char).collect();

    // Unmapped reads: no strand reversal needed, * if empty
    let qual_str = if qual.is_empty() {
        "*".to_string()
    } else {
        qual.iter().map(|&b| b as char).collect::<String>()
    };

    format!(
        "{}\t{}\t*\t0\t0\t*\t{}\t{}\t0\t{}\t{}",
        read_name, flag, rnext, pnext, seq_str, qual_str
    )
}

#[cfg(test)]
fn write_sam_header<W: Write>(
    writer: &mut W,
    ref_name: &str,
    ref_len: usize,
) -> std::io::Result<()> {
    writeln!(writer, "@HD\tVN:1.6\tSO:unsorted")?;
    writeln!(writer, "@SQ\tSN:{}\tLN:{}", ref_name, ref_len)?;
    writeln!(writer, "@PG\tID:freemap\tPN:freemap\tVN:0.0.1")?;
    Ok(())
}

// ============================================================================
// Main
// ============================================================================

fn print_usage() {
    eprintln!("freemap - Fast geometric sequence aligner");
    eprintln!();
    eprintln!("Usage: freemap [OPTIONS] <reference.fa> <reads.fq> <output.paf>");
    eprintln!(
        "       freemap [OPTIONS] -1 <R1.fq> -2 <R2.fq> <reference.fa> <output.sam>  # Paired-end"
    );
    eprintln!("       freemap -d <index.fmi> [OPTIONS] <reference.fa>   # Build index only");
    eprintln!();
    eprintln!("Presets (use -x, applied before other options):");
    eprintln!("  -x sr        short reads (Illumina): k=21, w=11, -R");
    eprintln!("  -x map-ont   ONT reads: k=15, w=10, f=20");
    eprintln!("  -x map-hifi  PacBio HiFi: k=19, w=19");
    eprintln!("  -x map-pb    PacBio CLR: k=17, w=10, f=20");
    eprintln!();
    eprintln!("Paired-end options:");
    eprintln!("  -1 FILE      first read file (R1)");
    eprintln!("  -2 FILE      second read file (R2)");
    eprintln!("  -I MIN:MAX   expected insert size range [0:1000]");
    eprintln!();
    eprintln!("Options:");
    eprintln!("  -k INT       k-mer size [19]");
    eprintln!("  -w INT       minimizer window size [19]");
    eprintln!("  -f INT       max k-mer frequency [200]");
    eprintln!("  -L INT       chaining lookback limit [16]");
    eprintln!("  -G INT       max gap difference for penalty [50]");
    eprintln!("  -S INT       gap penalty scaling factor [5]");
    eprintln!("  -t INT       threads [all]");
    eprintln!("  -c           generate detailed CIGAR (heuristic gap alignment)");
    eprintln!("  -g           TRAJECTORY mode: CIGAR from pure geometry (NO DP!)");
    eprintln!("  -a           output SAM format instead of PAF");
    eprintln!("  -H           HOMOPOLYMER compression mode (recommended for ONT!)");
    eprintln!("  -r           refine boundaries with micro-anchors (improves ONT accuracy)");
    eprintln!("  -R           SHORT READ mode: skip reverse strand if forward is strong (faster)");
    eprintln!("  -u           ULTRALONG mode: relaxed chaining for ONT ultralong reads");
    eprintln!("  -p           POLISH mode: refine gap CIGARs with O(n) greedy alignment");
    eprintln!("  --ransac-threshold FLOAT  RANSAC inlier threshold for trajectory regression [25.0]");
    eprintln!("  --polish-max-indel INT    Max indel size for CIGAR polishing [4]");
    eprintln!("  -d FILE      dump index to FILE (build index only, no alignment)");
    eprintln!("  -i FILE      load pre-built index from FILE");
    eprintln!("  -q           quiet mode");
    eprintln!("  -C           collision stats: count 32-bit hash collisions");
    eprintln!("  --multi      output secondary and supplementary alignments (for SV callers)");
    eprintln!("  --no-tiebreaker  disable primary-chromosome tiebreaker (ablation mode)");
    eprintln!("  --max-secondary N  max secondary alignments per read (default: 5)");
    eprintln!("  -h           show this help");
    eprintln!("  --version    show version");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  freemap -x sr -a ref.fa reads.fq out.sam             # Single-end Illumina");
    eprintln!("  freemap -x sr -a -1 R1.fq -2 R2.fq ref.fa out.sam    # Paired-end Illumina");
    eprintln!("  freemap -x map-ont -a ref.fa reads.fq out.sam        # ONT long reads");
    eprintln!("  freemap -x map-hifi -a ref.fa reads.fq out.sam       # PacBio HiFi");
    eprintln!("  freemap -x map-pb -a ref.fa reads.fq out.sam          # PacBio CLR");
    eprintln!();
    eprintln!("The -g flag enables trajectory-based alignment (CIGAR from geometry, no DP).");
    eprintln!("The -H flag enables homopolymer compression for ONT systematic errors.");
    eprintln!("The -p flag polishes trajectory CIGARs with base-level indel detection (DP-free).");
}

fn main() {
    let args: Vec<String> = env::args().collect();

    // Default values
    let mut k: usize = 19;
    let mut w: usize = 19;
    let mut max_freq: usize = 200;
    let mut max_lookback: usize = 16;
    let mut gap_max: i32 = 50;
    let mut gap_scale: i32 = 5;
    let mut threads: usize = 0;
    let mut quiet = false;
    let mut do_cigar = false;
    let mut sam_output = false;

    // Paired-end mode
    let mut r1_path: Option<String> = None;
    let mut r2_path: Option<String> = None;
    let mut insert_min: i64 = 0;
    let mut insert_max: i64 = 1000;
    let mut trajectory_mode = false;
    let mut refine_boundaries = false;
    let mut homopolymer_mode = false;
    let mut short_read_mode = false;
    let mut dump_index: Option<String> = None;
    let mut load_index_path: Option<String> = None;
    let mut collision_stats = false;
    let mut diagnostic_mode = false;
    let mut ultralong_mode = false;
    let mut custom_band_width: Option<i32> = None;
    let mut polish_mode = false; // Post-hoc DP refinement of large M blocks
    let mut ransac_threshold: f64 = 25.0; // RANSAC inlier threshold for trajectory regression
    let mut polish_max_indel: usize = 4; // Max indel size for CIGAR polishing
    let mut multi_mode = false; // Output secondary/supplementary alignments (--multi)
    let mut max_secondary: usize = 5; // Max secondary alignments (--max-secondary)
    let mut no_tiebreaker = false; // Disable primary-chromosome tiebreaker (--no-tiebreaker)
    let mut positional = Vec::new();

    // First pass: find preset (-x) and apply it
    let mut preset_name: Option<String> = None;
    for j in 1..args.len() {
        if args[j] == "-x" && j + 1 < args.len() {
            preset_name = Some(args[j + 1].clone());
            break;
        }
    }

    // Apply preset defaults FIRST (explicit options will override)
    if let Some(ref p) = preset_name {
        match p.as_str() {
            "sr" => {
                // Short reads (Illumina): optimized for 100-300bp reads
                k = 21;
                w = 11;
                max_freq = 500; // Higher threshold for repetitive regions in short reads
                short_read_mode = true;
                trajectory_mode = true;
                polish_mode = true; // Base-level gap comparison instead of geometric 'M'
                do_cigar = true;
            }
            "map-ont" => {
                // ONT reads: NO homopolymer compression (unlike PacBio)
                // minimap2 uses k=15, w=10 for ONT
                k = 15;
                w = 10;
                max_freq = 20; // Aggressive repetitive k-mer filtering for large genomes
                trajectory_mode = true;
                polish_mode = true; // Base-level gap comparison for better CIGAR
                do_cigar = true;
            }
            "map-hifi" => {
                // PacBio HiFi: high accuracy long reads (~1% error)
                // minimap2 uses k=19, w=19, NO -H
                k = 19;
                w = 19;
                max_freq = 20; // Aggressive repetitive k-mer filtering for large genomes
                trajectory_mode = true;
                polish_mode = true; // Base-level gap comparison for better CIGAR
                do_cigar = true;
            }
            "map-pb" | "map-clr" => {
                // PacBio CLR: higher error rate, longer reads
                // No homopolymer compression: HPC hurts accuracy on multi-chromosome
                // genomes due to coordinate decompression artifacts
                k = 17;
                w = 10;
                max_freq = 20; // Aggressive repetitive k-mer filtering for large genomes
                max_lookback = 64; // Larger lookback for CLR's higher error rate
                trajectory_mode = true;
                do_cigar = true;
            }
            _ => {}
        }
    }

    // Second pass: parse all options (overrides preset)
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-x" => {
                i += 1;
            } // Skip preset (already processed)
            "-k" => {
                i += 1;
                k = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(k);
            }
            "-w" => {
                i += 1;
                w = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(w);
            }
            "-f" => {
                i += 1;
                max_freq = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(max_freq);
            }
            "-L" => {
                i += 1;
                max_lookback = args
                    .get(i)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(max_lookback);
            }
            "-G" => {
                i += 1;
                gap_max = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(gap_max);
            }
            "-S" => {
                i += 1;
                gap_scale = args
                    .get(i)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(gap_scale);
            }
            "-t" => {
                i += 1;
                threads = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(0);
            }
            "-d" => {
                i += 1;
                dump_index = args.get(i).map(|s| s.to_string());
            }
            "-i" => {
                i += 1;
                load_index_path = args.get(i).map(|s| s.to_string());
            }
            "-1" => {
                i += 1;
                r1_path = args.get(i).map(|s| s.to_string());
            }
            "-2" => {
                i += 1;
                r2_path = args.get(i).map(|s| s.to_string());
            }
            "-I" => {
                i += 1;
                if let Some(s) = args.get(i) {
                    let parts: Vec<&str> = s.split(':').collect();
                    if parts.len() == 2 {
                        insert_min = parts[0].parse().unwrap_or(0);
                        insert_max = parts[1].parse().unwrap_or(1000);
                    }
                }
            }
            "-c" => {
                do_cigar = true;
            }
            "-g" => {
                trajectory_mode = true;
                do_cigar = true;
            }
            "-a" => {
                sam_output = true;
                do_cigar = true;
            }
            "-H" => {
                homopolymer_mode = true;
            }
            "--dp" => {
                // Optional DP polish: replace the greedy gap aligner with optimal affine-gap
                // DP on short inter-anchor gaps (better indel/mismatch placement for
                // base-level tasks like variant calling; DP is confined to short gaps so it
                // stays fast). Off by default (DP-free trajectory remains the default).
                DP_POLISH.store(true, std::sync::atomic::Ordering::Relaxed);
            }
            "-r" => {
                refine_boundaries = true;
            }
            "-R" => {
                short_read_mode = true;
            }
            "-q" => {
                quiet = true;
            }
            "-C" => {
                collision_stats = true;
            }
            "-D" => {
                diagnostic_mode = true;
            }
            "-u" => {
                ultralong_mode = true;
            }
            "--lookback" => {
                i += 1;
                max_lookback = args[i].parse().expect("--lookback requires integer");
            }
            "--band" => {
                i += 1;
                custom_band_width = Some(args[i].parse::<i32>().expect("--band requires integer"));
            }
            "-p" | "--polish" => {
                polish_mode = true;
                trajectory_mode = true;
                do_cigar = true;
            }
            "--ransac-threshold" => {
                i += 1;
                ransac_threshold = args[i].parse::<f64>().expect("--ransac-threshold requires a float");
            }
            "--polish-max-indel" => {
                i += 1;
                polish_max_indel = args[i].parse::<usize>().expect("--polish-max-indel requires an integer");
            }
            "--multi" => {
                multi_mode = true;
            }
            "--no-tiebreaker" => {
                no_tiebreaker = true;
            }
            "--max-secondary" => {
                i += 1;
                max_secondary = args[i].parse::<usize>().expect("--max-secondary requires an integer");
            }
            "-h" | "--help" => {
                print_usage();
                return;
            }
            "--version" => {
                println!("freemap {}", env!("CARGO_PKG_VERSION"));
                return;
            }
            s if !s.starts_with('-') => {
                positional.push(s.to_string());
            }
            _ => {
                eprintln!("Unknown option: {}", args[i]);
            }
        }
        i += 1;
    }

    // Detect paired-end mode
    let paired_end_mode = r1_path.is_some() && r2_path.is_some();

    // Print preset info
    if let Some(ref p) = preset_name {
        if !quiet {
            match p.as_str() {
                "sr" => eprintln!("[freemap] Preset 'sr': k={}, w={}, short-read mode", k, w),
                "map-ont" => eprintln!("[freemap] Preset 'map-ont': k={}, w={}", k, w),
                "map-hifi" => eprintln!("[freemap] Preset 'map-hifi': k={}, w={}", k, w),
                "map-pb" | "map-clr" => eprintln!(
                    "[freemap] Preset 'map-pb': k={}, w={}",
                    k, w
                ),
                _ => eprintln!("[freemap] Warning: unknown preset '{}', using defaults", p),
            }
        }
    }

    // Apply ultralong mode settings (for ONT ultralong reads)
    if ultralong_mode {
        max_lookback = 1024; // Much more lookback for sparse anchors with drift
        if !quiet {
            eprintln!("[freemap] ULTRALONG mode: max_lookback=1024, relaxed chaining");
        }
    }

    // For index-only mode (-d) or collision stats (-C), we only need reference
    let index_only_mode = (dump_index.is_some() || collision_stats) && positional.len() == 1;

    // Validate arguments based on mode
    if !index_only_mode {
        if paired_end_mode {
            // Paired mode: need reference and output
            if positional.len() < 2 {
                eprintln!("Error: Paired-end mode requires: -1 <R1.fq> -2 <R2.fq> <reference.fa> <output.sam>");
                print_usage();
                std::process::exit(1);
            }
        } else {
            // Single-end mode: need reference, reads, and output
            if positional.len() < 3 {
                print_usage();
                std::process::exit(1);
            }
        }
    }

    let ref_path = &positional[0];
    let reads_path = if !paired_end_mode && positional.len() > 1 {
        Some(&positional[1])
    } else {
        None
    };
    let output_path = if paired_end_mode {
        if positional.len() > 1 {
            Some(&positional[1])
        } else {
            None
        }
    } else {
        if positional.len() > 2 {
            Some(&positional[2])
        } else {
            None
        }
    };

    // Set thread count
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .ok();
    }
    let actual_threads = rayon::current_num_threads();

    // Load reference - either from embedded index (fast) or FASTA (slow)
    let start = Instant::now();
    let (global_ref, preloaded_index) = if let Some(ref index_path) = load_index_path {
        // Try loading index first to check for embedded reference
        match load_index_mmap(index_path, None) {
            Ok((mmap_idx, loaded_k, loaded_w, loaded_max_freq, loaded_homo))
                if mmap_idx.has_embedded_ref() =>
            {
                if !quiet {
                    eprintln!(
                        "[freemap] Loading embedded reference from index {}...",
                        index_path
                    );
                }
                let gref = mmap_idx.extract_reference().unwrap();
                (
                    gref,
                    Some((mmap_idx, loaded_k, loaded_w, loaded_max_freq, loaded_homo)),
                )
            }
            _ => {
                // No embedded ref or old format - need FASTA
                if !quiet {
                    eprintln!("[freemap] Loading reference...");
                }
                let ref_records = parse_fasta(ref_path);
                if ref_records.is_empty() {
                    eprintln!(
                        "[freemap] ERROR: No sequences found in reference file '{}'",
                        ref_path
                    );
                    eprintln!(
                        "[freemap] Make sure the file is in FASTA format (starts with '>')"
                    );
                    std::process::exit(1);
                }
                (GlobalReference::from_records(ref_records), None)
            }
        }
    } else {
        // Building index or no index specified - always need FASTA
        if !quiet {
            eprintln!("[freemap] Loading reference...");
        }
        let ref_records = parse_fasta(ref_path);
        if ref_records.is_empty() {
            eprintln!(
                "[freemap] ERROR: No sequences found in reference file '{}'",
                ref_path
            );
            eprintln!("[freemap] Make sure the file is in FASTA format (starts with '>')");
            std::process::exit(1);
        }
        (GlobalReference::from_records(ref_records), None)
    };
    let reference = &global_ref.sequence;

    if reference.is_empty() {
        eprintln!("[freemap] ERROR: Reference sequences are empty");
        std::process::exit(1);
    }
    if !quiet {
        eprintln!(
            "[freemap] Reference: {} chromosomes, {} bp total",
            global_ref.chroms.len(),
            reference.len()
        );
        if global_ref.chroms.len() <= 10 {
            for c in &global_ref.chroms {
                eprintln!("[freemap]   {} ({} bp)", c.name, c.len);
            }
        }
    }

    // Collision stats mode: just count collisions and exit
    if collision_stats {
        print_collision_stats(reference, k, w);
        return;
    }

    // Build or load index
    let index_start = Instant::now();

    // Build or load index - using IndexType enum to support both FxHashMap and MmapIndex
    let (index, compressed_ref, ref_pos_map): (IndexType, Vec<u8>, Vec<u32>) =
        if let Some((mmap_idx, loaded_k, loaded_w, _loaded_max_freq, loaded_homo)) =
            preloaded_index
        {
            // Index was already loaded with embedded reference (fast path)
            k = loaded_k;
            w = loaded_w;
            homopolymer_mode = loaded_homo;
            if !quiet {
                eprintln!(
                    "[freemap] Loaded index: k={}, w={}, homo={}",
                    k, w, homopolymer_mode
                );
            }
            (IndexType::Mmap(mmap_idx), Vec::new(), Vec::new())
        } else if let Some(ref index_path) = load_index_path {
            // Load pre-built index (v5, no embedded ref) - reference already loaded from FASTA
            if !quiet {
                eprintln!("[freemap] Loading pre-built index from {}...", index_path);
            }
            match load_index_mmap(index_path, Some(reference)) {
                Ok((mmap_idx, loaded_k, loaded_w, _loaded_max_freq, loaded_homo)) => {
                    k = loaded_k;
                    w = loaded_w;
                    homopolymer_mode = loaded_homo;
                    if !quiet {
                        eprintln!(
                            "[freemap] Loaded index: k={}, w={}, homo={}",
                            k, w, homopolymer_mode
                        );
                    }
                    (IndexType::Mmap(mmap_idx), Vec::new(), Vec::new())
                }
                Err(e) => {
                    eprintln!("[freemap] ERROR: {}", e);
                    std::process::exit(1);
                }
            }
        } else if homopolymer_mode {
        // Build homopolymer-compressed index
        if !quiet {
            eprintln!(
                "[freemap] Building HOMOPOLYMER-COMPRESSED index (k={}, w={})...",
                k, w
            );
        }

        // Compress reference first
        let (comp_ref, pos_map) = compress_homopolymers(reference);
        if !quiet {
            eprintln!(
                "[freemap] Compressed reference: {} bp → {} bp ({:.1}% reduction)",
                reference.len(),
                comp_ref.len(),
                100.0 * (1.0 - comp_ref.len() as f64 / reference.len() as f64)
            );
        }

        // Memory-efficient path: build SortedIndex, save, then reload as MmapIndex
        if let Some(ref index_path) = dump_index {
            if !quiet {
                eprintln!("[freemap] Saving index to {}...", index_path);
            }
            // Build SortedIndex (memory efficient)
            let idx_sorted = build_kmer_index(&comp_ref, k, w, max_freq);
            if let Err(e) = save_index(
                index_path,
                &idx_sorted,
                reference,
                &comp_ref,
                &pos_map,
                k,
                w,
                max_freq,
                true,
                &global_ref,
            ) {
                eprintln!("[freemap] WARNING: Failed to save index: {}", e);
            } else if !quiet {
                eprintln!("[freemap] Index saved successfully");
            }
            // Drop SortedIndex before loading MmapIndex to minimize peak memory
            drop(idx_sorted);

            // Reload as MmapIndex (memory-mapped, very efficient)
            match load_index_mmap(index_path, Some(reference)) {
                Ok((mmap_idx, _, _, _, _)) => (IndexType::Mmap(mmap_idx), Vec::new(), Vec::new()),
                Err(e) => {
                    eprintln!("[freemap] WARNING: Failed to reload index as mmap, falling back to HashMap: {}", e);
                    let idx_sorted = build_kmer_index(&comp_ref, k, w, max_freq);
                    (IndexType::Sorted(idx_sorted), comp_ref, pos_map)
                }
            }
        } else {
            // No saving needed - build SortedIndex directly (flat layout, low RAM)
            let idx_sorted = build_kmer_index(&comp_ref, k, w, max_freq);
            (IndexType::Sorted(idx_sorted), comp_ref, pos_map)
        }
    } else {
        // Build standard index
        if !quiet {
            eprintln!("[freemap] Building minimizer index (k={}, w={})...", k, w);
        }

        // Memory-efficient path: build SortedIndex, save, then reload as MmapIndex
        if let Some(ref index_path) = dump_index {
            if !quiet {
                eprintln!("[freemap] Saving index to {}...", index_path);
            }
            // Build SortedIndex (memory efficient)
            let idx_sorted = build_kmer_index(reference, k, w, max_freq);
            if let Err(e) = save_index(
                index_path,
                &idx_sorted,
                reference,
                &[],
                &[],
                k,
                w,
                max_freq,
                false,
                &global_ref,
            ) {
                eprintln!("[freemap] WARNING: Failed to save index: {}", e);
            } else if !quiet {
                eprintln!("[freemap] Index saved successfully");
            }
            // Drop SortedIndex before loading MmapIndex to minimize peak memory
            drop(idx_sorted);

            // Reload as MmapIndex (memory-mapped, very efficient)
            match load_index_mmap(index_path, Some(reference)) {
                Ok((mmap_idx, _, _, _, _)) => (IndexType::Mmap(mmap_idx), Vec::new(), Vec::new()),
                Err(e) => {
                    eprintln!("[freemap] WARNING: Failed to reload index as mmap, falling back to HashMap: {}", e);
                    let idx_sorted = build_kmer_index(reference, k, w, max_freq);
                    (IndexType::Sorted(idx_sorted), Vec::new(), Vec::new())
                }
            }
        } else {
            // No saving needed - build SortedIndex directly (flat layout, low RAM)
            let idx_sorted = build_kmer_index(reference, k, w, max_freq);
            (IndexType::Sorted(idx_sorted), Vec::new(), Vec::new())
        }
    };

    // Fix: when using MmapIndex in homopolymer mode, the compressed_ref and ref_pos_map Vecs
    // are empty (to avoid duplicating data already stored in the mmap). Extract slices from the
    // MmapIndex so alignment receives the actual compressed reference and position map.
    let (effective_compressed_ref, effective_pos_map): (Vec<u8>, Vec<u32>) =
        if compressed_ref.is_empty() && homopolymer_mode {
            match &index {
                IndexType::Mmap(ref mmap_idx) => {
                    (mmap_idx.compressed_ref().to_vec(), mmap_idx.pos_map().to_vec())
                }
                _ => (compressed_ref, ref_pos_map),
            }
        } else {
            (compressed_ref, ref_pos_map)
        };
    let compressed_ref = effective_compressed_ref;
    let ref_pos_map = effective_pos_map;

    let index_time = index_start.elapsed();
    if !quiet {
        eprintln!(
            "[freemap] Index: {} minimizers in {:.2}s",
            index.len(),
            index_time.as_secs_f64()
        );
    }

    // If index-only mode, exit here
    if index_only_mode {
        if !quiet {
            eprintln!("[freemap] Index-only mode complete.");
        }
        return;
    }

    // Load reads (single-end or paired-end)
    let load_reads_start = Instant::now();
    if !quiet {
        eprintln!("[freemap] Loading reads...");
    }

    let (reads, read_pairs): (Vec<FastqRecord>, Vec<ReadPair>) = if paired_end_mode {
        let pairs = load_paired_reads(r1_path.as_ref().unwrap(), r2_path.as_ref().unwrap());
        if !quiet {
            eprintln!(
                "[freemap] Loaded {} read pairs (paired-end mode), using {} threads",
                pairs.len(),
                actual_threads
            );
            eprintln!("[freemap] Insert size range: {}-{}", insert_min, insert_max);
        }
        (Vec::new(), pairs)
    } else {
        let reads = parse_fastq(reads_path.unwrap());
        if !quiet {
            eprintln!(
                "[freemap] Loaded {} reads, using {} threads",
                reads.len(),
                actual_threads
            );
        }
        (reads, Vec::new())
    };
    let load_reads_time = load_reads_start.elapsed();
    if !quiet {
        eprintln!(
            "[freemap] Read loading time: {:.2}s",
            load_reads_time.as_secs_f64()
        );
    }

    let total = if paired_end_mode {
        read_pairs.len() * 2
    } else {
        reads.len()
    };

    // Diagnostic mode: analyze anchor statistics
    if diagnostic_mode {
        eprintln!("[freemap] DIAGNOSTIC MODE: analyzing anchor statistics...");
        let sample_size = total.min(1000);
        let mut anchor_counts: Vec<usize> = Vec::with_capacity(sample_size);
        let mut read_lengths: Vec<usize> = Vec::with_capacity(sample_size);
        let mut chain_lengths: Vec<usize> = Vec::with_capacity(sample_size);

        for read in reads.iter().take(sample_size) {
            let seq = if homopolymer_mode {
                compress_homopolymers(&read.seq).0
            } else {
                read.seq.clone()
            };

            let max_occ = if ultralong_mode { 20 } else { 50 };
            let diag_limits =
                select_seeding_limits(k, w, max_occ, ultralong_mode, short_read_mode, false);
            let anchors = seed_read(&seq, &index, k, w, max_occ, &diag_limits);
            let n_anchors = anchors.len();
            anchor_counts.push(n_anchors);
            read_lengths.push(read.seq.len());

            if !anchors.is_empty() {
                let mut anchors_clone = anchors.clone();
                let diag_band = if ultralong_mode { 2048 } else { 64 };
                let chain_result = chain_anchors_topk(
                    &mut anchors_clone,
                    diag_band,
                    5,
                    max_lookback,
                    gap_max,
                    gap_scale,
                    ultralong_mode,
                );
                chain_lengths.push(chain_result.best_chain.len());
            }
        }

        anchor_counts.sort();
        read_lengths.sort();
        chain_lengths.sort();

        let median_idx = anchor_counts.len() / 2;
        let p90_idx = (anchor_counts.len() as f64 * 0.90) as usize;
        let p99_idx = (anchor_counts.len() as f64 * 0.99) as usize;

        eprintln!("=== ANCHOR DIAGNOSTICS (sample={}) ===", sample_size);
        eprintln!("Read lengths:");
        eprintln!("  Min:    {} bp", read_lengths.first().unwrap_or(&0));
        eprintln!(
            "  Median: {} bp",
            read_lengths.get(median_idx).unwrap_or(&0)
        );
        eprintln!("  P90:    {} bp", read_lengths.get(p90_idx).unwrap_or(&0));
        eprintln!("  Max:    {} bp", read_lengths.last().unwrap_or(&0));
        eprintln!("Anchors per read:");
        eprintln!("  Min:    {}", anchor_counts.first().unwrap_or(&0));
        eprintln!("  Median: {}", anchor_counts.get(median_idx).unwrap_or(&0));
        eprintln!("  P90:    {}", anchor_counts.get(p90_idx).unwrap_or(&0));
        eprintln!("  P99:    {}", anchor_counts.get(p99_idx).unwrap_or(&0));
        eprintln!("  Max:    {}", anchor_counts.last().unwrap_or(&0));
        if !chain_lengths.is_empty() {
            chain_lengths.sort();
            let ch_med = chain_lengths.len() / 2;
            let ch_p90 = (chain_lengths.len() as f64 * 0.90) as usize;
            eprintln!("Best chain length:");
            eprintln!("  Min:    {}", chain_lengths.first().unwrap_or(&0));
            eprintln!("  Median: {}", chain_lengths.get(ch_med).unwrap_or(&0));
            eprintln!("  P90:    {}", chain_lengths.get(ch_p90).unwrap_or(&0));
            eprintln!("  Max:    {}", chain_lengths.last().unwrap_or(&0));
        }
        eprintln!("=== END DIAGNOSTICS ===");
        return;
    }

    // Parallel alignment
    if !quiet {
        let mode_str = if homopolymer_mode && polish_mode {
            " with HOMOPOLYMER compression + POLISH"
        } else if homopolymer_mode {
            " with HOMOPOLYMER compression (no DP)"
        } else if polish_mode {
            " with TRAJECTORY mode + POLISH"
        } else if trajectory_mode && refine_boundaries {
            " with TRAJECTORY mode (no DP) + micro-anchor refinement"
        } else if trajectory_mode {
            " with TRAJECTORY mode (no DP)"
        } else if do_cigar {
            " with CIGAR"
        } else {
            ""
        };
        let pe_str = if paired_end_mode { " (paired-end)" } else { "" };
        eprintln!("[freemap] Aligning{}{}...", mode_str, pe_str);
        if no_tiebreaker {
            eprintln!("[freemap] Primary-chromosome tiebreaker DISABLED (--no-tiebreaker)");
        }
    }
    let align_start = Instant::now();
    let aligned_count = AtomicUsize::new(0);
    let band_width = custom_band_width.unwrap_or(if ultralong_mode {
        2048i32
    } else if short_read_mode {
        128i32 // Tighter diagonal band for 150bp reads (vs 512 for long reads)
    } else {
        512i32
    });

    // Results for paired-end or single-end mode
    type PairedResult = (
        String,
        Vec<u8>,
        Vec<u8>,
        Vec<u8>,
        Vec<u8>,
        Option<AlignmentResult>,
        Option<AlignmentResult>,
        bool,
        i64,
    );

    let mut multi_results: Vec<Option<(String, Vec<u8>, Vec<u8>, MultiAlignmentResult)>> = Vec::new();

    let (results, paired_results): (
        Vec<Option<(String, Vec<u8>, Vec<u8>, AlignmentResult)>>,
        Vec<PairedResult>,
    ) = if multi_mode && !paired_end_mode && !homopolymer_mode {
        // Multi-alignment mode: output secondary + supplementary alignments
        if !quiet {
            eprintln!("[freemap] Multi-alignment mode: max_secondary={}", max_secondary);
        }
        multi_results = reads
            .par_iter()
            .map(|read| {
                let result = align_read_multi(
                    &read.seq,
                    reference,
                    &index,
                    k,
                    w,
                    band_width,
                    do_cigar,
                    trajectory_mode,
                    refine_boundaries,
                    short_read_mode,
                    max_lookback,
                    gap_max,
                    gap_scale,
                    ultralong_mode,
                    max_freq,
                    Some(&global_ref),
                    false,
                    false,
                    polish_mode,
                    ransac_threshold,
                    polish_max_indel,
                    max_secondary,
                    no_tiebreaker,
                )?;
                aligned_count.fetch_add(1, Ordering::Relaxed);
                Some((read.name.clone(), read.seq.clone(), read.qual.clone(), result))
            })
            .collect();
        (Vec::new(), Vec::new())
    } else if paired_end_mode {
        // Paired-end alignment
        let paired_res: Vec<PairedResult> = read_pairs
            .par_iter()
            .map(|pair| {
                // Align R1
                let mut r1_result = align_read(
                    &pair.r1.seq,
                    reference,
                    &index,
                    k,
                    w,
                    band_width,
                    do_cigar,
                    trajectory_mode,
                    refine_boundaries,
                    short_read_mode,
                    max_lookback,
                    gap_max,
                    gap_scale,
                    ultralong_mode,
                    max_freq,
                    Some(&global_ref),
                    false,
                    false,
                    polish_mode,
                    ransac_threshold,
                    polish_max_indel,
                    no_tiebreaker,
                );
                if (!homopolymer_mode && (k <= 16 || short_read_mode))
                    && (r1_result.is_none()
                        || r1_result.as_ref().map(|r| r.mapq < 10).unwrap_or(true))
                {
                    r1_result = align_read(
                        &pair.r1.seq,
                        reference,
                        &index,
                        k,
                        w,
                        band_width,
                        do_cigar,
                        trajectory_mode,
                        refine_boundaries,
                        short_read_mode,
                        max_lookback,
                        gap_max,
                        gap_scale,
                        ultralong_mode,
                        max_freq,
                        Some(&global_ref),
                        true,
                        false,
                        polish_mode,
                        ransac_threshold,
                        polish_max_indel,
                        no_tiebreaker,
                    );
                }

                // Align R2
                let mut r2_result = align_read(
                    &pair.r2.seq,
                    reference,
                    &index,
                    k,
                    w,
                    band_width,
                    do_cigar,
                    trajectory_mode,
                    refine_boundaries,
                    short_read_mode,
                    max_lookback,
                    gap_max,
                    gap_scale,
                    ultralong_mode,
                    max_freq,
                    Some(&global_ref),
                    false,
                    false,
                    polish_mode,
                    ransac_threshold,
                    polish_max_indel,
                    no_tiebreaker,
                );
                if (!homopolymer_mode && (k <= 16 || short_read_mode))
                    && (r2_result.is_none()
                        || r2_result.as_ref().map(|r| r.mapq < 10).unwrap_or(true))
                {
                    r2_result = align_read(
                        &pair.r2.seq,
                        reference,
                        &index,
                        k,
                        w,
                        band_width,
                        do_cigar,
                        trajectory_mode,
                        refine_boundaries,
                        short_read_mode,
                        max_lookback,
                        gap_max,
                        gap_scale,
                        ultralong_mode,
                        max_freq,
                        Some(&global_ref),
                        true,
                        false,
                        polish_mode,
                        ransac_threshold,
                        polish_max_indel,
                        no_tiebreaker,
                    );
                }

                // Count aligned reads
                if r1_result.is_some() {
                    aligned_count.fetch_add(1, Ordering::Relaxed);
                }
                if r2_result.is_some() {
                    aligned_count.fetch_add(1, Ordering::Relaxed);
                }

                // Determine if proper pair and calculate insert size
                let (is_proper, insert_size) = match (&r1_result, &r2_result) {
                    (Some(r1), Some(r2)) => {
                        // Same chromosome check (both map to reference)
                        let (r1_chr, _) = global_ref.global_to_local(r1.ref_start);
                        let (r2_chr, _) = global_ref.global_to_local(r2.ref_start);
                        let same_chrom = r1_chr == r2_chr;

                        // Orientation check: R1 forward, R2 reverse (FR)
                        let correct_orientation = r1.strand == '+' && r2.strand == '-';

                        // Insert size: leftmost to rightmost
                        let isize = if r1.ref_start < r2.ref_start {
                            (r2.ref_end as i64) - (r1.ref_start as i64)
                        } else {
                            -((r1.ref_end as i64) - (r2.ref_start as i64))
                        };

                        let proper = same_chrom
                            && correct_orientation
                            && isize.abs() >= insert_min
                            && isize.abs() <= insert_max;
                        (proper, isize)
                    }
                    _ => (false, 0),
                };

                (
                    pair.name.clone(),
                    pair.r1.seq.clone(),
                    pair.r2.seq.clone(),
                    pair.r1.qual.clone(),
                    pair.r2.qual.clone(),
                    r1_result,
                    r2_result,
                    is_proper,
                    insert_size,
                )
            })
            .collect();
        (Vec::new(), paired_res)
    } else if homopolymer_mode {
        // Single-end homopolymer-compressed alignment
        let original_ref_len = reference.len();
        let single_res: Vec<Option<(String, Vec<u8>, Vec<u8>, AlignmentResult)>> = reads
            .par_iter()
            .map(|read| {
                let result = align_read_compressed(
                    &read.seq,
                    &compressed_ref,
                    &ref_pos_map,
                    reference,
                    &index,
                    k,
                    w,
                    band_width,
                    original_ref_len,
                    max_lookback,
                    gap_max,
                    gap_scale,
                    ultralong_mode,
                    max_freq,
                    Some(&global_ref),
                    false,
                    false,
                    polish_mode,
                    polish_max_indel,
                    no_tiebreaker,
                )?;
                aligned_count.fetch_add(1, Ordering::Relaxed);
                Some((read.name.clone(), read.seq.clone(), read.qual.clone(), result))
            })
            .collect();
        (single_res, Vec::new())
    } else {
        // Standard single-end alignment
        let single_res: Vec<Option<(String, Vec<u8>, Vec<u8>, AlignmentResult)>> = reads
            .par_iter()
            .map(|read| {
                let mut result = align_read(
                    &read.seq,
                    reference,
                    &index,
                    k,
                    w,
                    band_width,
                    do_cigar,
                    trajectory_mode,
                    refine_boundaries,
                    short_read_mode,
                    max_lookback,
                    gap_max,
                    gap_scale,
                    ultralong_mode,
                    max_freq,
                    Some(&global_ref),
                    false,
                    false,
                    polish_mode,
                    ransac_threshold,
                    polish_max_indel,
                    no_tiebreaker,
                )?;
                if (!homopolymer_mode && k <= 16) && result.mapq < 10 {
                    if let Some(recovered) = align_read(
                        &read.seq,
                        reference,
                        &index,
                        k,
                        w,
                        band_width,
                        do_cigar,
                        trajectory_mode,
                        refine_boundaries,
                        short_read_mode,
                        max_lookback,
                        gap_max,
                        gap_scale,
                        ultralong_mode,
                        max_freq,
                        Some(&global_ref),
                        true,
                        false,
                        polish_mode,
                        ransac_threshold,
                        polish_max_indel,
                        no_tiebreaker,
                    ) {
                        result = recovered;
                    }
                }
                aligned_count.fetch_add(1, Ordering::Relaxed);
                Some((read.name.clone(), read.seq.clone(), read.qual.clone(), result))
            })
            .collect();
        (single_res, Vec::new())
    };

    let align_time = align_start.elapsed();
    let aligned = aligned_count.load(Ordering::Relaxed);

    // Write output - convert global coordinates to local (per-chromosome)
    let output_start = Instant::now();
    let output_file = File::create(output_path.unwrap()).expect("Cannot create output");
    let mut writer = BufWriter::with_capacity(1 << 20, output_file);

    if sam_output {
        // Write SAM header with all chromosomes
        writeln!(writer, "@HD\tVN:1.6\tSO:unsorted").ok();
        for chrom in &global_ref.chroms {
            writeln!(writer, "@SQ\tSN:{}\tLN:{}", chrom.name, chrom.len).ok();
        }
        writeln!(writer, "@PG\tID:freemap\tPN:freemap\tVN:0.0.1").ok();

        if paired_end_mode {
            // Write paired-end SAM records
            for (name, r1_seq, r2_seq, r1_qual, r2_qual, r1_result, r2_result, is_proper, insert_size) in
                &paired_results
            {
                // Get chromosome info for both reads
                let r1_chrom = r1_result.as_ref().map(|r| {
                    let (idx, local_start) = global_ref.global_to_local(r.ref_start);
                    (
                        global_ref.chrom_name(idx),
                        local_start,
                        r.ref_end - r.ref_start,
                    )
                });
                let r2_chrom = r2_result.as_ref().map(|r| {
                    let (idx, local_start) = global_ref.global_to_local(r.ref_start);
                    (
                        global_ref.chrom_name(idx),
                        local_start,
                        r.ref_end - r.ref_start,
                    )
                });

                // Write R1
                if let Some(r1) = r1_result {
                    let (chrom_name, local_start, span) = r1_chrom.clone().unwrap();
                    let mut local_r1 = r1.clone();
                    local_r1.ref_start = local_start;
                    local_r1.ref_end = local_start + span;

                    let mate_chrom = r2_chrom.as_ref().map(|(c, _, _)| *c);
                    let mate_pos = r2_result.as_ref().map(|_r2| {
                        let (_, ls, _) = r2_chrom.clone().unwrap();
                        ls
                    });

                    let line = format_sam_line_paired(
                        &format!("{}/1", name),
                        r1_seq,
                        r1_qual,
                        chrom_name,
                        &local_r1,
                        mate_chrom,
                        mate_pos,
                        true,
                        *is_proper,
                        *insert_size,
                    );
                    writeln!(writer, "{}", line).ok();
                } else {
                    // R1 unmapped
                    let mate_chrom = r2_chrom.as_ref().map(|(c, _, _)| *c);
                    let mate_pos = r2_result.as_ref().map(|_| {
                        let (_, ls, _) = r2_chrom.clone().unwrap();
                        ls
                    });
                    let line = format_sam_line_unmapped_paired(
                        &format!("{}/1", name),
                        r1_seq,
                        r1_qual,
                        mate_chrom,
                        mate_pos,
                        true,
                    );
                    writeln!(writer, "{}", line).ok();
                }

                // Write R2
                if let Some(r2) = r2_result {
                    let (chrom_name, local_start, span) = r2_chrom.clone().unwrap();
                    let mut local_r2 = r2.clone();
                    local_r2.ref_start = local_start;
                    local_r2.ref_end = local_start + span;

                    let mate_chrom = r1_chrom.as_ref().map(|(c, _, _)| *c);
                    let mate_pos = r1_result.as_ref().map(|_| {
                        let (_, ls, _) = r1_chrom.clone().unwrap();
                        ls
                    });

                    let line = format_sam_line_paired(
                        &format!("{}/2", name),
                        r2_seq,
                        r2_qual,
                        chrom_name,
                        &local_r2,
                        mate_chrom,
                        mate_pos,
                        false,
                        *is_proper,
                        -*insert_size,
                    );
                    writeln!(writer, "{}", line).ok();
                } else {
                    // R2 unmapped
                    let mate_chrom = r1_chrom.as_ref().map(|(c, _, _)| *c);
                    let mate_pos = r1_result.as_ref().map(|_| {
                        let (_, ls, _) = r1_chrom.clone().unwrap();
                        ls
                    });
                    let line = format_sam_line_unmapped_paired(
                        &format!("{}/2", name),
                        r2_seq,
                        r2_qual,
                        mate_chrom,
                        mate_pos,
                        false,
                    );
                    writeln!(writer, "{}", line).ok();
                }
            }
        } else if multi_mode && !multi_results.is_empty() {
            // Write multi-alignment SAM records (primary + secondary + supplementary)
            let mut secondary_count = 0usize;
            let mut supplementary_count = 0usize;
            for opt in &multi_results {
                if let Some((name, seq, qual, multi)) = opt {
                    // Convert primary to local coordinates
                    let (pri_chrom_idx, pri_local_start) = global_ref.global_to_local(multi.primary.ref_start);
                    let pri_chrom_name = global_ref.chrom_name(pri_chrom_idx);
                    let mut local_primary = multi.primary.clone();
                    local_primary.ref_start = pri_local_start;
                    local_primary.ref_end = pri_local_start + (multi.primary.ref_end - multi.primary.ref_start);
                    if global_ref.is_blacklisted_contig(pri_chrom_idx) {
                        local_primary.mapq = 0;
                    }

                    // Convert supplementaries to local coordinates for SA tag
                    let local_supps: Vec<(String, AlignmentResult, usize, usize)> = multi.supplementaries.iter().map(|sup| {
                        let (idx, ls) = global_ref.global_to_local(sup.result.ref_start);
                        let cn = global_ref.chrom_name(idx).to_string();
                        let mut local_res = sup.result.clone();
                        local_res.ref_start = ls;
                        local_res.ref_end = ls + (sup.result.ref_end - sup.result.ref_start);
                        (cn, local_res, sup.clip_start, sup.clip_len)
                    }).collect();

                    // Build SA tag: list of all supplementary entries + primary (for supplementary records)
                    let primary_sa_entry = format_sa_entry(pri_chrom_name, &local_primary);
                    let supp_sa_entries: Vec<String> = local_supps.iter()
                        .map(|(cn, sup, _, _)| format_sa_entry(cn, sup))
                        .collect();

                    // SA tag on primary: lists all supplementaries
                    let sa_for_primary = if !supp_sa_entries.is_empty() {
                        let mut sa = supp_sa_entries.join(";");
                        sa.push(';');
                        Some(sa)
                    } else {
                        None
                    };

                    // Write primary
                    let line = format_sam_line_multi(
                        name, seq, qual, pri_chrom_name, &local_primary,
                        0, sa_for_primary.as_deref(),
                    );
                    writer.write_all(line.as_bytes()).ok();
                    writer.write_all(b"\n").ok();

                    // Write secondaries
                    for sec in &multi.secondaries {
                        let (sec_chrom_idx, sec_local_start) = global_ref.global_to_local(sec.ref_start);
                        let sec_chrom_name = global_ref.chrom_name(sec_chrom_idx);
                        let mut local_sec = sec.clone();
                        local_sec.ref_start = sec_local_start;
                        local_sec.ref_end = sec_local_start + (sec.ref_end - sec.ref_start);
                        let line = format_sam_line_multi(
                            name, seq, qual, sec_chrom_name, &local_sec,
                            0x100, None,
                        );
                        writer.write_all(line.as_bytes()).ok();
                        writer.write_all(b"\n").ok();
                        secondary_count += 1;
                    }

                    // Write supplementaries
                    // SA tag on supplementary: lists primary + other supplementaries
                    for (i, (sup_chrom_name, local_sup, clip_start, clip_len)) in local_supps.iter().enumerate() {
                        let mut sa_for_supp = primary_sa_entry.clone();
                        for (j, entry) in supp_sa_entries.iter().enumerate() {
                            if j != i {
                                sa_for_supp.push(';');
                                sa_for_supp.push_str(entry);
                            }
                        }
                        sa_for_supp.push(';');

                        let line = format_sam_line_supplementary(
                            name, seq, qual, sup_chrom_name, local_sup,
                            *clip_start, *clip_len, Some(&sa_for_supp),
                        );
                        writer.write_all(line.as_bytes()).ok();
                        writer.write_all(b"\n").ok();
                        supplementary_count += 1;
                    }
                }
            }
            if !quiet {
                eprintln!("[freemap] Multi-mode: {} secondary, {} supplementary alignments",
                    secondary_count, supplementary_count);
            }
        } else {
            // Write single-end SAM records - format in parallel, write sequentially
            let lines: Vec<Option<String>> = results.par_iter().map(|opt| {
                opt.as_ref().map(|(name, seq, qual, result)| {
                    let (chrom_idx, local_start) = global_ref.global_to_local(result.ref_start);
                    let chrom_name = global_ref.chrom_name(chrom_idx);
                    let mut local_result = result.clone();
                    local_result.ref_start = local_start;
                    local_result.ref_end = local_start + (result.ref_end - result.ref_start);
                    if global_ref.is_blacklisted_contig(chrom_idx) {
                        local_result.mapq = 0;
                    }
                    format_sam_line(name, seq, qual, chrom_name, &local_result)
                })
            }).collect();
            for line in &lines {
                if let Some(l) = line {
                    writer.write_all(l.as_bytes()).ok();
                    writer.write_all(b"\n").ok();
                }
            }
        }
    } else {
        // Write PAF records - format in parallel, write sequentially
        let lines: Vec<Option<String>> = results.par_iter().map(|opt| {
            opt.as_ref().map(|(name, seq, _qual, result)| {
                let (chrom_idx, local_start) = global_ref.global_to_local(result.ref_start);
                let chrom_name = global_ref.chrom_name(chrom_idx);
                let chrom_len = global_ref.chrom_len(chrom_idx);
                let mut local_result = result.clone();
                local_result.ref_start = local_start;
                local_result.ref_end = local_start + (result.ref_end - result.ref_start);
                if global_ref.is_blacklisted_contig(chrom_idx) {
                    local_result.mapq = 0;
                }
                format_paf_line(name, seq.len(), chrom_name, chrom_len, &local_result)
            })
        }).collect();
        for line in &lines {
            if let Some(l) = line {
                writer.write_all(l.as_bytes()).ok();
                writer.write_all(b"\n").ok();
            }
        }
    }
    writer.flush().ok();
    let output_time = output_start.elapsed();

    let total_time = start.elapsed();

    if !quiet {
        eprintln!();
        eprintln!("[freemap] === Results ===");
        eprintln!(
            "[freemap] Aligned: {}/{} ({:.1}%)",
            aligned,
            total,
            100.0 * aligned as f64 / total as f64
        );
        eprintln!("[freemap] Index time: {:.2}s", index_time.as_secs_f64());
        eprintln!(
            "[freemap] Read load time: {:.2}s",
            load_reads_time.as_secs_f64()
        );
        eprintln!("[freemap] Align time: {:.2}s", align_time.as_secs_f64());
        eprintln!("[freemap] Output time: {:.2}s", output_time.as_secs_f64());
        eprintln!("[freemap] Total time: {:.2}s", total_time.as_secs_f64());
        eprintln!(
            "[freemap] Throughput: {:.0} reads/sec",
            total as f64 / align_time.as_secs_f64()
        );

    }
}

// ============================================================================
// Unit Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write as IoWrite;

    // ============================================================================
    // Helper functions for tests
    // ============================================================================

    fn create_temp_file(content: &str, suffix: &str) -> String {
        let path = format!("/tmp/geomap_test_{}{}", std::process::id(), suffix);
        let mut file = File::create(&path).unwrap();
        file.write_all(content.as_bytes()).unwrap();
        path
    }

    fn cleanup_file(path: &str) {
        let _ = fs::remove_file(path);
    }

    // ============================================================================
    // Base conversion tests
    // ============================================================================

    #[test]
    fn test_complement() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'T'), b'A');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'G'), b'C');
        assert_eq!(complement(b'N'), b'N');
    }

    #[test]
    fn test_reverse_complement() {
        let seq = b"ACGT";
        let rc = reverse_complement(seq);
        assert_eq!(rc, b"ACGT"); // ACGT is self-complementary reversed

        let seq2 = b"AAACCC";
        let rc2 = reverse_complement(seq2);
        // AAACCC -> complement TTTGGG -> reversed GGGTTT
        assert_eq!(rc2, b"GGGTTT");
    }

    #[test]
    fn test_base_to_bits() {
        assert_eq!(base_to_bits(b'A'), 0);
        assert_eq!(base_to_bits(b'a'), 0);
        assert_eq!(base_to_bits(b'C'), 1);
        assert_eq!(base_to_bits(b'c'), 1);
        assert_eq!(base_to_bits(b'G'), 2);
        assert_eq!(base_to_bits(b'g'), 2);
        assert_eq!(base_to_bits(b'T'), 3);
        assert_eq!(base_to_bits(b't'), 3);
        assert_eq!(base_to_bits(b'N'), 4);
    }

    // ============================================================================
    // Homopolymer compression tests
    // ============================================================================

    #[test]
    fn test_compress_homopolymers_basic() {
        let seq = b"AAACCCGGGTTT";
        let (compressed, pos_map) = compress_homopolymers(seq);
        assert_eq!(compressed, b"ACGT");
        assert_eq!(pos_map.len(), 4);
        assert_eq!(pos_map[0], 0); // First A
        assert_eq!(pos_map[1], 3); // First C
        assert_eq!(pos_map[2], 6); // First G
        assert_eq!(pos_map[3], 9); // First T
    }

    #[test]
    fn test_compress_homopolymers_no_compression_needed() {
        let seq = b"ACGT";
        let (compressed, pos_map) = compress_homopolymers(seq);
        assert_eq!(compressed, b"ACGT");
        assert_eq!(pos_map, vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_compress_homopolymers_empty() {
        let seq = b"";
        let (compressed, pos_map) = compress_homopolymers(seq);
        assert!(compressed.is_empty());
        assert!(pos_map.is_empty());
    }

    #[test]
    fn test_decompress_pos() {
        let pos_map = vec![0u32, 3, 6, 9];
        assert_eq!(decompress_pos(0, &pos_map), 0);
        assert_eq!(decompress_pos(1, &pos_map), 3);
        assert_eq!(decompress_pos(2, &pos_map), 6);
        assert_eq!(decompress_pos(3, &pos_map), 9);
        // Beyond map - extrapolate
        assert_eq!(decompress_pos(4, &pos_map), 10);
    }

    // ============================================================================
    // K-mer and indexing tests
    // ============================================================================

    #[test]
    fn test_kmer_hash_deterministic() {
        let kmer1 = b"ACGT";
        let kmer2 = b"ACGT";
        assert_eq!(kmer_hash(kmer1), kmer_hash(kmer2));
    }

    #[test]
    fn test_kmer_hash_different() {
        let kmer1 = b"ACGT";
        let kmer2 = b"TGCA";
        assert_ne!(kmer_hash(kmer1), kmer_hash(kmer2));
    }

    #[test]
    fn test_build_kmer_index_basic() {
        let reference = b"ACGTACGTACGT";
        let index = build_kmer_index(reference, 4, 2, 100);
        // Should have entries for k-mers
        assert!(!index.is_empty());
    }

    #[test]
    fn test_build_kmer_index_with_n() {
        let reference = b"ACGTNNNNNACGT";
        let index = build_kmer_index(reference, 4, 2, 100);
        // N's should be skipped
        assert!(!index.is_empty());
    }

    #[test]
    fn test_build_kmer_index_empty() {
        let reference = b"AC"; // Too short for k=4
        let index = build_kmer_index(reference, 4, 2, 100);
        assert!(index.is_empty());
    }

    // ============================================================================
    // Seeding tests
    // ============================================================================

    #[test]
    fn test_seed_read_finds_matches() {
        let reference = b"ACGTACGTACGTACGTACGT";
        let index = build_kmer_index(reference, 4, 2, 100);
        let read = b"ACGTACGT";
        let limits = select_seeding_limits(4, 2, 50, false, false, false);
        let anchors = seed_read(read, &index, 4, 2, 50, &limits);
        // Should find some anchors
        assert!(!anchors.is_empty());
    }

    #[test]
    fn test_seed_read_no_matches() {
        let reference = b"AAAAAAAAAAAAAAAAAAAAAA";
        let index = build_kmer_index(reference, 4, 2, 100);
        let read = b"TTTTTTTT";
        let limits = select_seeding_limits(4, 2, 50, false, false, false);
        let _anchors = seed_read(read, &index, 4, 2, 50, &limits);
        // Should not find matches (different sequences)
        // Note: might find empty due to index filtering
    }

    // ============================================================================
    // Chaining tests
    // ============================================================================

    #[test]
    fn test_chain_anchors_empty() {
        let mut anchors: Vec<Anchor> = vec![];
        let result = chain_anchors_topk(&mut anchors, 64, 5, 16, 50, 5, false);
        assert!(result.best_chain.is_empty());
        assert_eq!(result.best_score, 0);
    }

    #[test]
    fn test_chain_anchors_single() {
        let mut anchors = vec![Anchor {
            read_start: 0,
            ref_start: 100,
            len: 15,
        }];
        let result = chain_anchors_topk(&mut anchors, 64, 5, 16, 50, 5, false);
        assert_eq!(result.best_chain.len(), 1);
        assert_eq!(result.best_score, 15);
    }

    #[test]
    fn test_chain_anchors_colinear() {
        let mut anchors = vec![
            Anchor {
                read_start: 0,
                ref_start: 100,
                len: 15,
            },
            Anchor {
                read_start: 20,
                ref_start: 120,
                len: 15,
            },
            Anchor {
                read_start: 40,
                ref_start: 140,
                len: 15,
            },
        ];
        let result = chain_anchors_topk(&mut anchors, 64, 5, 16, 50, 5, false);
        // All anchors should be chained together
        assert!(result.best_score > 15);
    }

    // ============================================================================
    // MAPQ calculation tests
    // ============================================================================

    #[test]
    fn test_topological_class_unique() {
        let chain_result = ChainResult {
            best_chain: vec![Anchor {
                read_start: 0,
                ref_start: 100,
                len: 15,
            }],
            best_score: 100,
            alternatives: vec![],
        };
        let class = classify_ambiguity(&chain_result, 1000);
        assert_eq!(class, TopologicalClass::Unique);
    }

    #[test]
    fn test_calculate_mapq_unique() {
        let chain_result = ChainResult {
            best_chain: vec![Anchor {
                read_start: 0,
                ref_start: 100,
                len: 15,
            }],
            best_score: 100,
            alternatives: vec![],
        };
        let mapq = calculate_topological_mapq(&chain_result, 150);
        // Should be high MAPQ for unique alignment
        assert!(mapq >= 30);
    }

    #[test]
    fn test_calculate_mapq_ambiguous() {
        let chain1 = vec![
            Anchor {
                read_start: 0,
                ref_start: 100,
                len: 15,
            },
            Anchor {
                read_start: 20,
                ref_start: 120,
                len: 15,
            },
            Anchor {
                read_start: 40,
                ref_start: 140,
                len: 15,
            },
        ];
        let chain2 = vec![
            Anchor {
                read_start: 0,
                ref_start: 200,
                len: 15,
            },
            Anchor {
                read_start: 20,
                ref_start: 220,
                len: 15,
            },
            Anchor {
                read_start: 40,
                ref_start: 240,
                len: 15,
            },
        ];
        let chain_result = ChainResult {
            best_chain: chain1.clone(),
            best_score: 45,
            alternatives: vec![(chain2, 45, 1)],
        };
        let mapq = calculate_topological_mapq(&chain_result, 150);
        // Should be lower MAPQ for ambiguous alignment
        assert!(mapq < 60);
    }

    // ============================================================================
    // CIGAR generation tests
    // ============================================================================

    #[test]
    fn test_cigar_to_string() {
        let ops = vec![('M', 10), ('I', 2), ('M', 5), ('D', 3), ('M', 8)];
        let cigar = cigar_to_string(&ops);
        assert_eq!(cigar, "10M2I5M3D8M");
    }

    #[test]
    fn test_cigar_to_string_empty() {
        let ops: Vec<(char, usize)> = vec![];
        let cigar = cigar_to_string(&ops);
        assert_eq!(cigar, "");
    }

    #[test]
    fn test_add_cigar_op_merge() {
        let mut ops: Vec<(char, usize)> = vec![('M', 5)];
        add_cigar_op(&mut ops, 'M', 3);
        assert_eq!(ops, vec![('M', 8)]);
    }

    #[test]
    fn test_add_cigar_op_no_merge() {
        let mut ops: Vec<(char, usize)> = vec![('M', 5)];
        add_cigar_op(&mut ops, 'I', 2);
        assert_eq!(ops, vec![('M', 5), ('I', 2)]);
    }

    // ============================================================================
    // Gap alignment tests
    // ============================================================================

    #[test]
    fn test_fast_gap_align_equal() {
        let query = b"ACGT";
        let target = b"ACGT";
        let ops = fast_gap_align(query, target);
        // Should be simple match
        assert!(!ops.is_empty());
    }

    #[test]
    fn test_fast_gap_align_insertion() {
        let query = b"ACGTACGT";
        let target = b"ACGT";
        let ops = fast_gap_align(query, target);
        // Query longer = insertion
        let total_len: usize = ops.iter().map(|(_, l)| *l).sum();
        assert!(total_len >= 4);
    }

    #[test]
    fn test_fast_gap_align_deletion() {
        let query = b"ACGT";
        let target = b"ACGTACGT";
        let ops = fast_gap_align(query, target);
        // Target longer = deletion
        let has_deletion = ops.iter().any(|(op, _)| *op == 'D');
        assert!(has_deletion);
    }

    #[test]
    fn test_fast_gap_align_empty() {
        let query: &[u8] = b"";
        let target: &[u8] = b"";
        let ops = fast_gap_align(query, target);
        assert!(ops.is_empty());
    }

    // ============================================================================
    // diagonal_align_refined() tests
    // ============================================================================

    #[test]
    fn test_diagonal_align_refined_perfect_match() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTACGTACGT";
        let ops = diagonal_align_refined(query, target, 4);
        // Perfect match → single M block
        assert_eq!(ops.len(), 1);
        assert_eq!(ops[0], ('M', 12));
    }

    #[test]
    fn test_diagonal_align_refined_mismatch() {
        let query = b"ACGTACGTACGT";
        let target = b"ACGTANGTACGT";
        let ops = diagonal_align_refined(query, target, 4);
        // Should contain only M ops (mismatches counted as M in this aligner)
        let total_read: usize = ops.iter().filter(|(op, _)| *op == 'M' || *op == 'I').map(|(_, l)| *l).sum();
        let total_ref: usize = ops.iter().filter(|(op, _)| *op == 'M' || *op == 'D').map(|(_, l)| *l).sum();
        assert_eq!(total_read, 12);
        assert_eq!(total_ref, 12);
    }

    #[test]
    fn test_diagonal_align_refined_small_insertion() {
        // Read has 2bp insertion in the middle
        let target = b"AAAAACCCCCGGGGGTTTTTT"; // 21bp
        let query = b"AAAAACCTTCCCGGGGGTTTTTT"; // 23bp (2bp inserted after pos 7)
        let ops = diagonal_align_refined(query, target, 4);
        // Should detect the insertion
        let has_ins = ops.iter().any(|(op, _)| *op == 'I');
        assert!(has_ins, "should detect insertion: {:?}", ops);
        let total_read: usize = ops.iter().filter(|(op, _)| *op == 'M' || *op == 'I').map(|(_, l)| *l).sum();
        assert_eq!(total_read, query.len(), "read consumption mismatch");
    }

    #[test]
    fn test_diagonal_align_refined_small_deletion() {
        // Read has 2bp deletion in the middle
        let target = b"AAAAACCCCCGGGGGTTTTTT"; // 21bp
        let query = b"AAAAACCCGGGGGTTTTTT"; // 19bp (2bp deleted)
        let ops = diagonal_align_refined(query, target, 4);
        // Should detect the deletion
        let has_del = ops.iter().any(|(op, _)| *op == 'D');
        assert!(has_del, "should detect deletion: {:?}", ops);
        let total_ref: usize = ops.iter().filter(|(op, _)| *op == 'M' || *op == 'D').map(|(_, l)| *l).sum();
        assert_eq!(total_ref, target.len(), "ref consumption mismatch");
    }

    #[test]
    fn test_diagonal_align_refined_empty_inputs() {
        assert!(diagonal_align_refined(b"", b"", 4).is_empty());
        assert_eq!(diagonal_align_refined(b"", b"ACGT", 4), vec![('D', 4)]);
        assert_eq!(diagonal_align_refined(b"ACGT", b"", 4), vec![('I', 4)]);
    }

    #[test]
    fn test_nw_affine_cigar() {
        // helpers
        let qlen = |ops: &[(char, usize)]| -> usize {
            ops.iter().filter(|(o, _)| matches!(o, 'M' | 'I')).map(|(_, l)| l).sum()
        };
        let tlen = |ops: &[(char, usize)]| -> usize {
            ops.iter().filter(|(o, _)| matches!(o, 'M' | 'D')).map(|(_, l)| l).sum()
        };
        // identical -> single M run
        assert_eq!(nw_affine_cigar(b"ACGTACGT", b"ACGTACGT"), vec![('M', 8)]);
        // single mismatch stays M (M does not distinguish match/mismatch)
        assert_eq!(nw_affine_cigar(b"ACGTACGT", b"ACGAACGT"), vec![('M', 8)]);
        // clean single-base deletion (target longer): query aligns with one D
        let d = nw_affine_cigar(b"ACGTACGT", b"ACGTTACGT");
        assert_eq!(qlen(&d), 8);
        assert_eq!(tlen(&d), 9);
        assert_eq!(d.iter().filter(|(o, _)| *o == 'D').map(|(_, l)| *l).sum::<usize>(), 1);
        // clean single-base insertion (query longer)
        let ins = nw_affine_cigar(b"ACGTTACGT", b"ACGTACGT");
        assert_eq!(qlen(&ins), 9);
        assert_eq!(tlen(&ins), 8);
        assert_eq!(ins.iter().filter(|(o, _)| *o == 'I').map(|(_, l)| *l).sum::<usize>(), 1);
        // empties
        assert!(nw_affine_cigar(b"", b"").is_empty());
        assert_eq!(nw_affine_cigar(b"", b"ACGT"), vec![('D', 4)]);
        assert_eq!(nw_affine_cigar(b"ACGT", b""), vec![('I', 4)]);
        // affine gaps: a contiguous 3bp deletion should be ONE D3, not split
        let d3 = nw_affine_cigar(b"ACGTACGT", b"ACGTGGGACGT");
        assert_eq!(qlen(&d3), 8);
        assert_eq!(tlen(&d3), 11);
        assert!(d3.iter().any(|(o, l)| *o == 'D' && *l == 3));
        // conservation on random-ish cases
        let cases: Vec<(&[u8], &[u8])> = vec![
            (b"ACGTACGTACGT", b"ACGTACGTACGT"),
            (b"ACGTTTACGT", b"ACGTACGT"),
            (b"ACGTACGT", b"ACGTTTACGT"),
            (b"AAACCCGGG", b"AAACCGGG"),
        ];
        for (q, t) in cases {
            let ops = nw_affine_cigar(q, t);
            assert_eq!(qlen(&ops), q.len());
            assert_eq!(tlen(&ops), t.len());
        }
    }

    #[test]
    fn test_diagonal_align_refined_conservation() {
        // For any input, total M+I must equal query len, total M+D must equal target len
        let cases: Vec<(&[u8], &[u8])> = vec![
            (b"ACGTACGT", b"ACGTACGT"),       // equal
            (b"ACGTTTACGT", b"ACGTACGT"),     // query longer (insertion)
            (b"ACGTACGT", b"ACGTTTACGT"),     // target longer (deletion)
            (b"AAAAAAAAAA", b"CCCCCCCCCC"),   // all mismatches
            (b"A", b"ACGT"),                   // very short query
        ];
        for (q, t) in cases {
            let ops = diagonal_align_refined(q, t, 4);
            let total_read: usize = ops.iter().filter(|(op, _)| *op == 'M' || *op == 'I').map(|(_, l)| *l).sum();
            let total_ref: usize = ops.iter().filter(|(op, _)| *op == 'M' || *op == 'D').map(|(_, l)| *l).sum();
            assert_eq!(total_read, q.len(), "read conservation failed for q={} t={}: {:?}", q.len(), t.len(), ops);
            assert_eq!(total_ref, t.len(), "ref conservation failed for q={} t={}: {:?}", q.len(), t.len(), ops);
        }
    }

    // ============================================================================
    // RANSAC regression tests
    // ============================================================================

    #[test]
    fn test_ransac_perfect_line() {
        let points = vec![(0.0, 100.0), (10.0, 110.0), (20.0, 120.0), (30.0, 130.0)];
        let (slope, intercept) = ransac_regression(&points, 50, 5.0);
        // Slope should be ~1.0, intercept ~100.0
        assert!((slope - 1.0).abs() < 0.1);
        assert!((intercept - 100.0).abs() < 5.0);
    }

    #[test]
    fn test_ransac_with_outliers() {
        let points = vec![
            (0.0, 100.0),
            (10.0, 110.0),
            (20.0, 120.0),
            (30.0, 130.0),
            (15.0, 500.0), // outlier
        ];
        let (slope, _intercept) = ransac_regression(&points, 100, 10.0);
        // Should still find ~1.0 slope despite outlier
        assert!((slope - 1.0).abs() < 0.2);
    }

    #[test]
    fn test_ransac_single_point() {
        let points = vec![(10.0, 110.0)];
        let (slope, intercept) = ransac_regression(&points, 50, 5.0);
        // With single point, should return slope=1.0
        assert_eq!(slope, 1.0);
        assert_eq!(intercept, 100.0);
    }

    // ============================================================================
    // Chains equivalence tests
    // ============================================================================

    #[test]
    fn test_chains_not_equivalent_different_scores() {
        let chain1 = vec![
            Anchor {
                read_start: 0,
                ref_start: 100,
                len: 15,
            },
            Anchor {
                read_start: 20,
                ref_start: 120,
                len: 15,
            },
            Anchor {
                read_start: 40,
                ref_start: 140,
                len: 15,
            },
        ];
        let chain2 = vec![
            Anchor {
                read_start: 0,
                ref_start: 200,
                len: 15,
            },
            Anchor {
                read_start: 20,
                ref_start: 220,
                len: 15,
            },
        ];
        // Different scores (3*15=45 vs 2*15=30)
        assert!(!chains_are_equivalent(&chain1, &chain2, 45, 30));
    }

    // ============================================================================
    // Index serialization tests
    // ============================================================================

    #[test]
    fn test_save_and_load_index_standard() {
        let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let k = 15;
        let w = 10;
        let max_freq = 100;
        let global_ref = GlobalReference {
            sequence: reference.to_vec(),
            chroms: vec![ChromInfo {
                name: "seq".to_string(),
                offset: 0,
                len: reference.len(),
            }],
        };

        // Build index
        let index = build_kmer_index(reference, k, w, max_freq);
        let compressed_ref: Vec<u8> = vec![];
        let pos_map: Vec<u32> = vec![];

        // Save index
        let index_path = format!("/tmp/geo_test_index_{}.fmi", std::process::id());
        save_index(
            &index_path,
            &index,
            reference,
            &compressed_ref,
            &pos_map,
            k,
            w,
            max_freq,
            false,
            &global_ref,
        )
        .unwrap();

        // Load index using mmap
        let (mmap_index, loaded_k, loaded_w, loaded_max_freq, loaded_homo) =
            load_index_mmap(&index_path, Some(reference)).unwrap();

        // Verify parameters
        assert_eq!(loaded_k, k);
        assert_eq!(loaded_w, w);
        assert_eq!(loaded_max_freq, max_freq);
        assert!(!loaded_homo);

        // Verify lookups work - check that we can find entries that exist in original
        for (hash, positions) in index.iter() {
            if let Some(mmap_positions) = mmap_index.get_positions(hash) {
                // Positions should match
                assert_eq!(positions, mmap_positions);
            }
        }

        // Cleanup
        cleanup_file(&index_path);
    }

    #[test]
    fn test_load_index_wrong_reference() {
        let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let k = 15;
        let w = 10;
        let max_freq = 100;
        let global_ref = GlobalReference {
            sequence: reference.to_vec(),
            chroms: vec![ChromInfo {
                name: "seq".to_string(),
                offset: 0,
                len: reference.len(),
            }],
        };

        // Build and save index
        let index = build_kmer_index(reference, k, w, max_freq);
        let index_path = format!("/tmp/geo_test_index_wrong_{}.fmi", std::process::id());
        save_index(
            &index_path,
            &index,
            reference,
            &[],
            &[],
            k,
            w,
            max_freq,
            false,
            &global_ref,
        )
        .unwrap();

        // Try to load with different reference
        let different_ref = b"TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA";
        let result = load_index_mmap(&index_path, Some(different_ref));
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("mismatch"));

        // Cleanup
        cleanup_file(&index_path);
    }

    #[test]
    fn test_load_index_invalid_file() {
        let reference = b"ACGT";
        let result = load_index_mmap("/nonexistent/path/to/index.fmi", Some(reference));
        assert!(result.is_err());
    }

    // ============================================================================
    // Full alignment pipeline tests
    // ============================================================================

    #[test]
    fn test_align_read_finds_position() {
        // Simple test: read is exact substring of reference
        // Use longer sequences to ensure k-mer matches
        let reference = b"NNNNNNNNNNACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNNNNNN";
        let read = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        let index = build_kmer_index(reference, 15, 10, 100);
        let result = align_read(
            read, reference, &index, 15, 10, 64, false, false, false, false, 16, 50, 5, false, 100,
            None, false, false, false, 25.0, 4, false,
        );

        // Should find an alignment
        assert!(result.is_some());
        let aln = result.unwrap();
        // Just verify it found something reasonable
        assert!(aln.ref_end > aln.ref_start);
    }

    #[test]
    fn test_align_read_reverse_complement() {
        // Read is reverse complement of reference substring
        // ACGT is self-complementary, so we use a non-palindromic sequence
        let reference = b"NNNNNNNNNNAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTNNNNNNNNNN";
        let fwd_seq = b"AAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTT";
        let read = reverse_complement(fwd_seq);

        let index = build_kmer_index(reference, 15, 10, 100);
        let result = align_read(
            &read, reference, &index, 15, 10, 64, false, false, false, false, 16, 50, 5, false,
            100, None, false, false, false, 25.0, 4, false,
        );

        // Should find an alignment (either strand)
        assert!(result.is_some());
    }

    // ============================================================================
    // FASTA/FASTQ parsing tests
    // ============================================================================

    #[test]
    fn test_parse_fasta() {
        let fasta_content = ">seq1 description\nACGT\nTGCA\n>seq2\nAAAA\n";
        let path = create_temp_file(fasta_content, ".fa");

        let records = parse_fasta(&path);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].0, "seq1");
        assert_eq!(records[0].1, b"ACGTTGCA");
        assert_eq!(records[1].0, "seq2");
        assert_eq!(records[1].1, b"AAAA");

        cleanup_file(&path);
    }

    #[test]
    fn test_parse_fastq() {
        let fastq_content = "@read1 comment\nACGT\n+\nIIII\n@read2\nTTTT\n+\nJJJJ\n";
        let path = create_temp_file(fastq_content, ".fq");

        let records = parse_fastq(&path);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "read1");
        assert_eq!(records[0].seq, b"ACGT");
        assert_eq!(records[1].name, "read2");
        assert_eq!(records[1].seq, b"TTTT");

        cleanup_file(&path);
    }

    // ============================================================================
    // PAF output format tests
    // ============================================================================

    #[test]
    fn test_format_paf_line() {
        let result = AlignmentResult {
            ref_start: 100,
            ref_end: 200,
            read_start: 0,
            read_end: 100,
            cigar: "100M".to_string(),
            chain_score: 95,
            mapq: 60,
            strand: '+',
            _chrom_idx: 0,
        };

        let paf = format_paf_line("read1", 100, "chr1", 1000, &result);
        let fields: Vec<&str> = paf.split('\t').collect();

        assert_eq!(fields[0], "read1"); // Query name
        assert_eq!(fields[1], "100"); // Query length
        assert_eq!(fields[2], "0"); // Query start
        assert_eq!(fields[3], "100"); // Query end
        assert_eq!(fields[4], "+"); // Strand
        assert_eq!(fields[5], "chr1"); // Target name
        assert_eq!(fields[6], "1000"); // Target length
        assert_eq!(fields[7], "100"); // Target start
        assert_eq!(fields[8], "200"); // Target end
        assert_eq!(fields[11], "60"); // MAPQ
        assert!(fields[12].starts_with("cg:Z:")); // CIGAR tag
    }

    // ============================================================================
    // Trajectory-based alignment tests
    // ============================================================================

    #[test]
    fn test_align_trajectory_basic() {
        let anchors = vec![
            Anchor {
                read_start: 0,
                ref_start: 100,
                len: 15,
            },
            Anchor {
                read_start: 20,
                ref_start: 120,
                len: 15,
            },
            Anchor {
                read_start: 40,
                ref_start: 140,
                len: 15,
            },
            Anchor {
                read_start: 60,
                ref_start: 160,
                len: 15,
            },
        ];
        let read = vec![b'A'; 100];
        let reference = vec![b'A'; 300];

        let result = align_trajectory_based(&anchors, &read, &reference, 60, 60, '+', false, false, 25.0, 4);
        assert!(result.is_some());

        let aln = result.unwrap();
        assert!(!aln.cigar.is_empty());
    }

    #[test]
    fn test_align_trajectory_too_few_anchors() {
        let anchors = vec![Anchor {
            read_start: 0,
            ref_start: 100,
            len: 15,
        }];
        let read = vec![b'A'; 100];
        let reference = vec![b'A'; 300];

        let result = align_trajectory_based(&anchors, &read, &reference, 15, 60, '+', false, false, 25.0, 4);
        // Should return None with < 2 anchors
        assert!(result.is_none());
    }

    // ============================================================================
    // SAM output tests
    // ============================================================================

    #[test]
    fn test_format_sam_line_forward() {
        let result = AlignmentResult {
            ref_start: 100,
            ref_end: 250,
            read_start: 0,
            read_end: 150,
            strand: '+',
            chain_score: 300,
            mapq: 60,
            cigar: "150M".to_string(),
            _chrom_idx: 0,
        };
        let seq = b"ACGTACGTACGT";
        let qual = b"IIIIIIIIIIII";
        let line = format_sam_line("read1", seq, qual, "chr1", &result);

        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[0], "read1"); // QNAME
        assert_eq!(fields[1], "0"); // FLAG (forward)
        assert_eq!(fields[2], "chr1"); // RNAME
        assert_eq!(fields[3], "101"); // POS (1-based)
        assert_eq!(fields[4], "60"); // MAPQ
        assert_eq!(fields[5], "150M"); // CIGAR
        assert_eq!(fields[9], "ACGTACGTACGT"); // SEQ
        assert_eq!(fields[10], "IIIIIIIIIIII"); // QUAL
    }

    #[test]
    fn test_format_sam_line_reverse() {
        let result = AlignmentResult {
            ref_start: 100,
            ref_end: 250,
            read_start: 0,
            read_end: 150,
            strand: '-',
            chain_score: 300,
            mapq: 55,
            cigar: "150M".to_string(),
            _chrom_idx: 0,
        };
        let seq = b"AACC";
        let qual = b"ABCD";
        let line = format_sam_line("read2", seq, qual, "chr1", &result);

        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[0], "read2");
        assert_eq!(fields[1], "16"); // FLAG (reverse)
        assert_eq!(fields[4], "55"); // MAPQ
        assert_eq!(fields[9], "GGTT"); // Reverse complement of AACC
        assert_eq!(fields[10], "DCBA"); // Reversed quality (not complemented)
    }

    #[test]
    fn test_parse_clip_sizes() {
        assert_eq!(parse_clip_sizes("100S50M"), (100, 0));
        assert_eq!(parse_clip_sizes("50M100S"), (0, 100));
        assert_eq!(parse_clip_sizes("50S100M200S"), (50, 200));
        assert_eq!(parse_clip_sizes("150M"), (0, 0));
        assert_eq!(parse_clip_sizes("100H50M"), (100, 0));
        assert_eq!(parse_clip_sizes("50M100H"), (0, 100));
    }

    #[test]
    fn test_format_sa_entry() {
        let result = AlignmentResult {
            ref_start: 999,
            ref_end: 1500,
            read_start: 0,
            read_end: 500,
            strand: '+',
            chain_score: 400,
            mapq: 60,
            cigar: "500M".to_string(),
            _chrom_idx: 0,
        };
        let entry = format_sa_entry("chr1", &result);
        assert_eq!(entry, "chr1,1000,+,500M,60,0");
    }

    #[test]
    fn test_format_sam_line_multi_primary() {
        let result = AlignmentResult {
            ref_start: 100,
            ref_end: 250,
            read_start: 0,
            read_end: 150,
            strand: '+',
            chain_score: 300,
            mapq: 60,
            cigar: "150M".to_string(),
            _chrom_idx: 0,
        };
        let line = format_sam_line_multi("read1", b"ACGT", b"IIII", "chr1", &result, 0, None);
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[0], "read1");
        assert_eq!(fields[1], "0"); // primary, forward
    }

    #[test]
    fn test_format_sam_line_multi_secondary() {
        let result = AlignmentResult {
            ref_start: 500,
            ref_end: 650,
            read_start: 0,
            read_end: 150,
            strand: '-',
            chain_score: 200,
            mapq: 0,
            cigar: "150M".to_string(),
            _chrom_idx: 0,
        };
        let line = format_sam_line_multi("read1", b"ACGT", b"IIII", "chr1", &result, 0x100, None);
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[1], "272"); // 0x100 | 0x10 = 256 + 16 = 272
    }

    #[test]
    fn test_format_sam_line_multi_supplementary() {
        let result = AlignmentResult {
            ref_start: 2000,
            ref_end: 2500,
            read_start: 0,
            read_end: 500,
            strand: '+',
            chain_score: 400,
            mapq: 30,
            cigar: "100S400M".to_string(),
            _chrom_idx: 0,
        };
        let sa = "chr1,101,+,150M,60,0;";
        let line = format_sam_line_multi("read1", b"ACGT", b"IIII", "chr1", &result, 0x800, Some(sa));
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[1], "2048"); // 0x800
        // CIGAR keeps soft clips (full SEQ is emitted)
        assert_eq!(fields[5], "100S400M");
        // SEQ should be present (SV callers need it)
        assert_eq!(fields[9], "ACGT");
        // Should have SA tag
        assert!(line.contains("SA:Z:chr1,101,+,150M,60,0;"));
    }

    #[test]
    fn test_sam_header() {
        let mut output = Vec::new();
        write_sam_header(&mut output, "chr1", 1000000).unwrap();
        let header = String::from_utf8(output).unwrap();

        assert!(header.contains("@HD\tVN:1.6"));
        assert!(header.contains("@SQ\tSN:chr1\tLN:1000000"));
        assert!(header.contains("@PG\tID:freemap"));
    }

    #[test]
    fn test_mmap_index_preserves_hpc_data() {
        // Build a multi-chromosome reference with homopolymer runs
        // chr1: 60bp with runs, chr2: 60bp with runs
        let chr1 = b"AAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTT";
        let chr2 = b"TTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAA";
        let mut reference = Vec::with_capacity(chr1.len() + chr2.len());
        reference.extend_from_slice(chr1);
        reference.extend_from_slice(chr2);

        let global_ref = GlobalReference {
            sequence: reference.clone(),
            chroms: vec![
                ChromInfo { name: "chr1".to_string(), offset: 0, len: chr1.len() },
                ChromInfo { name: "chr2".to_string(), offset: chr1.len(), len: chr2.len() },
            ],
        };

        // Compress reference
        let (comp_ref, pos_map) = compress_homopolymers(&reference);
        assert!(!comp_ref.is_empty());
        assert!(!pos_map.is_empty());

        // Build index on compressed reference
        let k = 15;
        let w = 10;
        let max_freq = 100;
        let idx = build_kmer_index(&comp_ref, k, w, max_freq);

        // Save index with HPC data
        let index_path = format!("/tmp/geo_test_hpc_mmap_{}.fmi", std::process::id());
        save_index(&index_path, &idx, &reference, &comp_ref, &pos_map, k, w, max_freq, true, &global_ref)
            .unwrap();

        // Load as MmapIndex
        let (mmap_idx, _, _, _, loaded_homo) = load_index_mmap(&index_path, Some(&reference)).unwrap();
        assert!(loaded_homo);

        // Verify compressed_ref and pos_map are preserved (not empty)
        let mmap_comp_ref = mmap_idx.compressed_ref();
        let mmap_pos_map = mmap_idx.pos_map();
        assert_eq!(mmap_comp_ref.len(), comp_ref.len(), "compressed_ref should be preserved in mmap");
        assert_eq!(mmap_pos_map.len(), pos_map.len(), "pos_map should be preserved in mmap");
        assert_eq!(mmap_comp_ref, comp_ref.as_slice());
        for i in 0..pos_map.len() {
            assert_eq!(mmap_pos_map[i], pos_map[i], "pos_map mismatch at index {}", i);
        }

        // Verify that decompressed positions map to correct chromosomes
        // A position in the second half of compressed space should map to chr2
        let mid_compressed = comp_ref.len() / 2;
        let decompressed = decompress_pos(mid_compressed, &pos_map);
        let (_chrom_idx, _) = global_ref.global_to_local(decompressed as usize);
        // The decompressed position should be valid (within reference bounds)
        assert!((decompressed as usize) < reference.len(),
            "decompressed position {} should be within reference length {}", decompressed, reference.len());

        // Verify that empty Vec (the bug) would give wrong results
        let empty_map: Vec<u32> = vec![];
        let wrong_pos = decompress_pos(mid_compressed, &empty_map);
        // With empty map, decompress_pos returns the compressed position unchanged
        assert_eq!(wrong_pos as usize, mid_compressed,
            "empty pos_map should return compressed position unchanged (the bug)");
        // The wrong position should differ from the correct one (proving the bug matters)
        assert_ne!(wrong_pos, decompressed,
            "empty pos_map gives wrong position: {} vs correct {}", wrong_pos, decompressed);

        // Cleanup
        cleanup_file(&index_path);
    }

    #[test]
    fn test_embedded_reference_roundtrip() {
        // Test that reference can be saved in the index and extracted without FASTA
        let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let k = 15;
        let w = 10;
        let max_freq = 100;
        let global_ref = GlobalReference {
            sequence: reference.to_vec(),
            chroms: vec![
                ChromInfo { name: "chr1".to_string(), offset: 0, len: 24 },
                ChromInfo { name: "chr2".to_string(), offset: 24, len: 23 },
            ],
        };

        let index = build_kmer_index(reference, k, w, max_freq);
        let index_path = format!("/tmp/geo_test_embedded_ref_{}.fmi", std::process::id());
        save_index(&index_path, &index, reference, &[], &[], k, w, max_freq, false, &global_ref)
            .unwrap();

        // Load without external reference (standalone mode)
        let (mmap_idx, loaded_k, _, _, _) = load_index_mmap(&index_path, None).unwrap();
        assert_eq!(loaded_k, k);
        assert!(mmap_idx.has_embedded_ref());

        // Extract reference and verify
        let extracted = mmap_idx.extract_reference().unwrap();
        assert_eq!(extracted.sequence, reference.to_vec());
        assert_eq!(extracted.chroms.len(), 2);
        assert_eq!(extracted.chroms[0].name, "chr1");
        assert_eq!(extracted.chroms[0].offset, 0);
        assert_eq!(extracted.chroms[0].len, 24);
        assert_eq!(extracted.chroms[1].name, "chr2");
        assert_eq!(extracted.chroms[1].offset, 24);
        assert_eq!(extracted.chroms[1].len, 23);

        // Verify lookups still work
        for (hash, positions) in index.iter() {
            if let Some(mmap_positions) = mmap_idx.get_positions(hash) {
                assert_eq!(positions, mmap_positions);
            }
        }

        cleanup_file(&index_path);
    }
}
