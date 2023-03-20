# Docker Image for running MfiX-Exa ctests

FROM ubuntu:20.04

LABEL maintainer="Justin Weber <justin.weber@netl.doe.gov>"

ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE 1
ENV TZ=America/New_York
ENV DEBIAN_FRONTEND=noninteractive

SHELL ["/bin/bash", "-o", "pipefail", "-c"]


RUN apt-get -qq update\
    && apt-get -qq upgrade -y\
    && apt-get -qq -y install --no-install-recommends \
    ca-certificates \
    curl=7.* \
    gnupg=2.* \
    wget \
    dpkg \
    && wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-keyring_1.0-1_all.deb \
    && dpkg -i cuda-keyring_1.0-1_all.deb \
    && apt-get -qq update \
    && apt-get -qq -y install --no-install-recommends \
    build-essential=12.* \
    ccache=3.* \
    clang-12=1:12.* \
    cuda-command-line-tools-11-3=11.3.* \
    cuda-compiler-11-3=11.3.* \
    cuda-cupti-dev-11-3=11.3.* \
    cuda-minimal-build-11-3=11.3.* \
    cuda-nvml-dev-11-3=11.3.* \
    cuda-nvtx-11-3=11.3.* \
    g++=4:9.* gfortran-9=9.* \
    git=1:2.* \
    libcurand-dev-11-3=10.* \
    libopenmpi-dev=4.* \
    openmpi-bin=4.* \
    pkg-config=0.29.* \
    python-is-python3=3.* \
    python3-pip=20.* \
    python3-setuptools=45.* \
    python3-venv=3.* \
    shellcheck=0.* \
    && apt-get clean && rm -rf /var/lib/apt/lists/* \
    && pip3 install --no-cache-dir \
    cmake==3.* \
    codespell==2.2.1 \
    ninja==1.* \
    && curl -L https://github.com/mozilla/sccache/releases/download/v0.2.15/sccache-v0.2.15-x86_64-unknown-linux-musl.tar.gz -- \
    | tar xz \
    && cp sccache-v0.2.15-x86_64-unknown-linux-musl/sccache /usr/local/bin/sccache \
    && rm -rf sccache-v0.2.15-x86_64-unknown-linux-musl/ \
    && chmod +x /usr/local/bin/sccache

# install csg dependencies
RUN apt-get -qq update\
    && apt-get -qq -y install --no-install-recommends \
    libgmp-dev \
    libmpfr-dev \
    libboost-all-dev

RUN git clone --depth 1 --branch v2.13.7 https://github.com/catchorg/Catch2
WORKDIR /Catch2
RUN cmake -S . -B build -DCATCH_BUILD_TESTING=OFF
RUN cmake --build build --target install
WORKDIR /

RUN git clone --depth 1 --branch v5.3 https://github.com/CGAL/cgal
WORKDIR /cgal
RUN cmake -S . -B build
RUN cmake --build build --target install
WORKDIR /

RUN git clone --branch 3.2.2 https://github.com/taocpp/PEGTL
WORKDIR /PEGTL
RUN cmake -S . -B build -DPEGTL_BUILD_TESTS=OFF -DPEGTL_BUILD_EXAMPLES=OFF
RUN cmake --build build --target install
WORKDIR /
    

ENV PATH="/usr/local/cuda-11.3/bin:${PATH}"

RUN printf '#!/bin/bash\nenv $@' > /usr/local/bin/srun \
    && chmod +x /usr/local/bin/srun \
    && useradd --create-home -s /bin/bash user

WORKDIR /home/user
USER user

CMD [ "/bin/bash" ]
