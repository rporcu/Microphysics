# Docker Image for running MfiX-Exa ctests

FROM ubuntu:20.04

LABEL maintainer="Justin Weber <justin.weber@netl.doe.gov>"

ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE 1
ENV DEBIAN_FRONTEND=noninteractive

SHELL ["/bin/bash", "-o", "pipefail", "-c"]


RUN apt-get -qq update\
    && apt-get -qq upgrade -y\
    && apt-get -qq -y install --no-install-recommends \
    ca-certificates=20210119~20.* \
    curl=7.* \
    gnupg=2.* \
    && curl -L https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/7fa2af80.pub > 7fa2af80.pub \
    && apt-key add 7fa2af80.pub \
    && rm 7fa2af80.pub \
    && echo "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64 /" > /etc/apt/sources.list.d/cuda.list \
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
    codespell==2.* \
    conan==1.* \
    ninja==1.* \
    && curl -L https://github.com/mozilla/sccache/releases/download/v0.2.15/sccache-v0.2.15-x86_64-unknown-linux-musl.tar.gz -- \
    | tar xz \
    && cp sccache-v0.2.15-x86_64-unknown-linux-musl/sccache /usr/local/bin/sccache \
    && rm -rf sccache-v0.2.15-x86_64-unknown-linux-musl/ \
    && chmod +x /usr/local/bin/sccache

ENV PATH="/usr/local/cuda-11.3/bin:${PATH}"

#RUN printf '#!/bin/bash\nenv $@' > /usr/local/bin/srun \
#    && chmod +x /usr/local/bin/srun \

RUN useradd --create-home -s /bin/bash user

WORKDIR /home/user
USER user

CMD [ "/bin/bash" ]
