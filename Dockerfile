# Docker Image for running MfiX-Exa ctests

FROM ubuntu:20.04

LABEL maintainer="Mark Meredith <mark.meredith@netl.doe.gov>"

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -qq update \
    && apt-get -qq -y install --no-install-recommends \
    build-essential=12.8ubuntu1.1 \
    ccache=3.7.7-1 \
    git=1:2.25.1-1ubuntu3 \
    libopenmpi-dev=4.0.3-0ubuntu1 \
    openmpi-bin=4.0.3-0ubuntu1 \
    python3-pip=20.0.2-5ubuntu1.1 \
    python3-setuptools=45.2.0-1 \
    python3-venv=3.8.2-0ubuntu2 \
    shellcheck=0.7.0-2build2 \
    && rm -rf /var/lib/apt

RUN pip3 install --no-cache-dir \
    cmake==3.18.4 \
    codespell==2.0.0 \
    conan==1.34.0 \
    ninja==1.10.0.post2

RUN printf '#!/bin/bash\nenv $@' > /usr/local/bin/srun
RUN chmod +x /usr/local/bin/srun

RUN useradd --create-home -s /bin/bash user
WORKDIR /home/user
USER user

CMD [ "/bin/bash" ]
