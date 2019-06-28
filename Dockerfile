# Docker Image for runnning MfiX-Exa ctests

FROM ubuntu:19.04

MAINTAINER Mark Meredith <mark.meredith@netl.doe.gov>

RUN apt-get -qq update \
  && apt-get -qq -y install \
  autoconf \
  automake \
  ccache \
  cmake \
  g++ \
  gcc \
  gfortran \
  git \
  libopenmpi-dev \
  make \
  openmpi-bin \
  python \
  wget

RUN wget -qO- "https://cmake.org/files/v3.14/cmake-3.14.1-Linux-x86_64.tar.gz" | tar --strip-components=1 -xz -C /usr/local
RUN printf '#!/bin/bash\nenv $@' > /usr/local/bin/srun
RUN chmod +x /usr/local/bin/srun

CMD [ "/bin/bash" ]
