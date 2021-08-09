# Docker Image for running MfiX-Exa ctests

FROM ubuntu:20.04

LABEL maintainer="Mark Meredith <mark.meredith@netl.doe.gov>"

ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE 1
ENV DEBIAN_FRONTEND=noninteractive

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get -qq update \
    && apt-get -qq upgrade \
    && apt-get -qq -y install --no-install-recommends \
    ca-certificates=20210119~20.04.1 \
    curl=7.68.0-1ubuntu2.6 \
    gnupg=2.2.19-3ubuntu2.1 \
    && curl -L https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/7fa2af80.pub > 7fa2af80.pub \
    && apt-key add 7fa2af80.pub \
    && rm 7fa2af80.pub \
    && echo "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64 /" > /etc/apt/sources.list.d/cuda.list \
    && apt-get -qq update \
    && apt-get -qq -y install --no-install-recommends \
    build-essential=12.8ubuntu1.1 \
    ccache=3.7.7-1 \
    cuda-command-line-tools-11-3=11.3.0-1 \
    cuda-compiler-11-3=11.3.0-1 \
    cuda-cupti-dev-11-3=11.3.58-1 \
    cuda-minimal-build-11-3=11.3.0-1 \
    cuda-nvml-dev-11-3=11.3.58-1 \
    cuda-nvtx-11-3=11.3.58-1 \
    g++=4:9.3.0-1ubuntu2 gfortran-9=9.3.0-17ubuntu1~20.04 \
    git=1:2.25.1-1ubuntu3.1 \
    libcurand-dev-11-3=10.2.4.58-1 \
    libopenmpi-dev=4.0.3-0ubuntu1 \
    openmpi-bin=4.0.3-0ubuntu1 \
    pkg-config=0.29.1-0ubuntu4 \
    python3-pip=20.0.2-5ubuntu1.6 \
    python3-setuptools=45.2.0-1 \
    python3-venv=3.8.2-0ubuntu2 \
    shellcheck=0.7.0-2build2 \
    && apt-get clean && rm -rf /var/lib/apt/lists/* \
    && pip3 install --no-cache-dir \
    cmake==3.21.0 \
    codespell==2.1.0 \
    conan==1.38.0 \
    ninja==1.10.0.post3 \
    && curl -L https://github.com/mozilla/sccache/releases/download/v0.2.15/sccache-v0.2.15-x86_64-unknown-linux-musl.tar.gz -- \
    | tar xz \
    && cp sccache-v0.2.15-x86_64-unknown-linux-musl/sccache /usr/local/bin/sccache \
    && rm -rf sccache-v0.2.15-x86_64-unknown-linux-musl/ \
    && chmod +x /usr/local/bin/sccache

ENV PATH="/usr/local/cuda-11.3/bin:${PATH}"

RUN printf '#!/bin/bash\nenv $@' > /usr/local/bin/srun \
    && chmod +x /usr/local/bin/srun \
    && useradd --create-home -s /bin/bash user

WORKDIR /home/user
USER user

CMD [ "/bin/bash" ]
