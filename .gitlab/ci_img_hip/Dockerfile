# Docker Image for running MfiX-Exa ctests

FROM ubuntu:20.04

LABEL maintainer="Justin Weber <justin.weber@netl.doe.gov>"

ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE 1
ENV TZ=America/New_York
ENV DEBIAN_FRONTEND=noninteractive

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# update and add rocm repo
RUN apt-get -qq update\
    && apt-get -qq -y install --no-install-recommends \
    apt-utils \
    ca-certificates \
    curl \
    gnupg \
    wget \
    && wget -q -O - https://repo.radeon.com/rocm/rocm.gpg.key | apt-key add - \
    && echo 'deb [arch=amd64] https://repo.radeon.com/rocm/apt/debian/ ubuntu main' | tee /etc/apt/sources.list.d/rocm.list \
    && apt-get -qq update \
    && apt-get -qq upgrade -y

# install python  first, rocm-llvm requires it?
RUN apt-get install -y --no-install-recommends \
    python-is-python3 \
    python3 \
    python3-pip \
    python3-setuptools \
    python3-venv

# install rocm
RUN apt-get install -y --no-install-recommends \
    build-essential \
    ccache \
    clang-12 \
    gfortran \
    git \
    libnuma-dev \
    libopenmpi-dev \
    llvm-amdgpu \
    openmpi-bin \
    pkg-config \
    rocm-dev \
    roctracer-dev \
    rocprofiler-dev \
    rocrand-dev \
    rocprim-dev

# python build tools
RUN apt-get clean && rm -rf /var/lib/apt/lists/* \
    && pip3 install --no-cache-dir \
    cmake==3.* \
    conan==1.* \
    ninja==1.*

# add rocm to path
ENV PATH="/opt/rocm/llvm/bin:/opt/rocm/bin:/opt/rocm/profiler/bin:/opt/rocm/opencl/bin:${PATH}"

RUN printf '#!/bin/bash\nenv $@' > /usr/local/bin/srun \
    && chmod +x /usr/local/bin/srun \
    && useradd --create-home -s /bin/bash user

WORKDIR /home/user
USER user

CMD [ "/bin/bash" ]
