# Docker Image for runnning MfiX-Exa ctests

FROM nvidia/cuda:10.2-devel-ubuntu18.04

MAINTAINER Mark Meredith <mark.meredith@netl.doe.gov>

RUN apt-get -qq update \
  && apt-get -qq -y install \
  ccache \
  gfortran \
  git \
  libopenmpi-dev \
  openmpi-bin \
  python \
  wget

RUN wget https://bootstrap.pypa.io/get-pip.py
RUN python get-pip.py
RUN python -m pip install cmake
RUN python -m pip install ninja

RUN printf '#!/bin/bash\nenv $@' > /usr/local/bin/srun
RUN chmod +x /usr/local/bin/srun

RUN useradd --create-home -s /bin/bash user
WORKDIR /home/user
USER user

CMD [ "/bin/bash" ]
