# Docker Image for Intel DPC++ compiler

FROM ubuntu:20.04
LABEL maintainer="Mark Meredith <mark.meredith@netl.doe.gov>"

ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE 1
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update \
  && apt-get -qq -y install --no-install-recommends \
  ca-certificates=20210119~20.04.1 \
  git=1:2.25.1-1ubuntu3 \
  gnupg=2.2.19-3ubuntu2.1 \
  wget=1.20.3-1ubuntu1 \
  && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN wget --no-check-certificate https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB \
    && apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB \
    && rm GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB
RUN echo "deb https://apt.repos.intel.com/oneapi all main" \
  > /etc/apt/sources.list.d/oneAPI.list

RUN apt-get -qq update \
  && apt-get -qq -y install --no-install-recommends \
  build-essential=12.8ubuntu1.1 \
  cmake=3.16.3-1ubuntu1 \
  g++=4:9.3.0-1ubuntu2 \
  gfortran=4:9.3.0-1ubuntu2 \
  intel-oneapi-dpcpp-cpp-compiler \
  intel-oneapi-mkl-devel=2021.1.1-52 \
  libopenmpi-dev=4.0.3-0ubuntu1 \
  ninja-build=1.10.0-1build1 \
  openmpi-bin=4.0.3-0ubuntu1 \
  python3-venv=3.8.2-0ubuntu2 \
  && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN useradd --create-home -s /bin/bash user
WORKDIR /home/user
USER user

CMD [ "/bin/bash" ]
