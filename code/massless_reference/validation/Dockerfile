# pySecDec numerical validator for the massive masters (no Mathematica license).
FROM python:3.11-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
      g++ make libgsl-dev zlib1g-dev ca-certificates \
      m4 bzip2 lbzip2 xz-utils file perl autoconf libtool \
    && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir numpy sympy pySecDec
WORKDIR /work
