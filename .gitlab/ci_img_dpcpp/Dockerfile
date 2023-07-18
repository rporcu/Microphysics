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
    && wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null \
    && echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list \
    && apt-get -qq update \
    && apt-get -qq upgrade -y

# install python  first
RUN apt-get install -y --no-install-recommends \
    python-is-python3 \
    python3 \
    python3-pip \
    python3-setuptools \
    python3-venv

# install dpcpp
RUN apt-get install -y --no-install-recommends \
    build-essential \
    ccache \
    clang-12=1:12.* \
    g++ \
    gfortran \
    git \
    intel-oneapi-compiler-dpcpp-cpp \
    intel-oneapi-mkl-devel \
    libopenmpi-dev \
    openmpi-bin \
    pkg-config

# python build tools
RUN apt-get clean && rm -rf /var/lib/apt/lists/* \
    && pip3 install --no-cache-dir \
    cmake==3.* \
    ninja==1.*

RUN printf '#!/bin/bash\nenv $@' > /usr/local/bin/srun \
    && chmod +x /usr/local/bin/srun \
    && useradd --create-home -s /bin/bash user

WORKDIR /home/user
USER user

CMD [ "/bin/bash" ]
