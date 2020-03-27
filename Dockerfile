# Docker Image for running MfiX-Exa ctests

FROM nvidia/cuda:10.2-devel-ubuntu18.04

MAINTAINER Mark Meredith <mark.meredith@netl.doe.gov>

RUN apt-get -qq update \
  && apt-get -qq -y install \
  ccache \
  gfortran \
  git \
  libopenmpi-dev \
  openmpi-bin \
  python3 \
  python3-setuptools \
  wget

RUN wget https://bootstrap.pypa.io/get-pip.py
RUN python3 get-pip.py
RUN python3 -m pip install cmake
RUN python3 -m pip install ninja

RUN printf '#!/bin/bash\nenv $@' > /usr/local/bin/srun
RUN chmod +x /usr/local/bin/srun
RUN ln -s /usr/bin/python3 /usr/local/bin/python

RUN useradd --create-home -s /bin/bash user
WORKDIR /home/user
USER user

CMD [ "/bin/bash" ]
