# Docker Image for AMD ROCm

FROM ubuntu:20.04

LABEL maintainer="Mark Meredith <mark.meredith@netl.doe.gov>"

ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE 1
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update \
  && apt-get -qq -y install --no-install-recommends \
  ca-certificates=20210119~20.04.1 \
  gnupg=2.2.19-3ubuntu2.1 \
  wget=1.20.3-1ubuntu1 \
  && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN wget --no-check-certificate http://repo.radeon.com/rocm/rocm.gpg.key \
    && apt-key add rocm.gpg.key \
    && rm rocm.gpg.key

RUN echo 'deb [arch=amd64] http://repo.radeon.com/rocm/apt/debian/ xenial main' \
  > /etc/apt/sources.list.d/rocm.list

RUN echo "export PATH=$PATH:/opt/rocm/bin:/opt/rocm/profiler/bin:/opt/rocm/opencl/bin" \
  > /etc/profile.d/rocm.sh

RUN apt-get -qq update \
  && apt-get -qq -y install --no-install-recommends \
  build-essential=12.8ubuntu1.1 \
  cmake=3.16.3-1ubuntu1 \
  gfortran=4:9.3.0-1ubuntu2 \
  git=1:2.25.1-1ubuntu3 \
  libnuma-dev=2.0.12-1 \
  libopenmpi-dev=4.0.3-0ubuntu1 \
  ninja-build=1.10.0-1build1 \
  openmpi-bin=4.0.3-0ubuntu1 \
  python3-pip=20.0.2-5ubuntu1.1 \
  python3-venv=3.8.2-0ubuntu2 \
  rocm-dev=4.1.0.40100-26 \
  rocrand=2.10.7.40100-26 \
  && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN pip3 install conan==1.32.1

RUN useradd --create-home -s /bin/bash user
WORKDIR /home/user
USER user

CMD [ "/bin/bash" ]
