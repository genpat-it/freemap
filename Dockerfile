FROM rust:1.83-slim AS builder

WORKDIR /build
COPY Cargo.toml Cargo.lock ./
COPY src/ src/

RUN cargo build --release

FROM debian:bookworm-slim

COPY --from=builder /build/target/release/freemap /usr/local/bin/freemap

ENTRYPOINT ["freemap"]
CMD ["-h"]
