# modified from https://github.com/ContinuumIO/docker-images/blob/master/anaconda/Dockerfile

FROM ubuntu:16.04

MAINTAINER Mark Meredith <mark.meredith@netl.doe.gov>

RUN apt-get update --fix-missing && apt-get install -y wget ca-certificates \
    git gfortran cmake make numdiff openmpi-bin libopenmpi-dev bzip2

RUN apt-get install -y curl grep sed tcsh dpkg && \
    TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

ENV PATH /opt/conda/bin:$PATH

# http://bugs.python.org/issue19846
# > At the moment, setting "LANG=C" on a Linux system *fundamentally breaks Python 3*, and that's not OK.
ENV LANG C.UTF-8

ENTRYPOINT [ "/usr/bin/tini", "--" ]
RUN groupadd -r gitlab && useradd -r -g gitlab gitlab
USER gitlab
CMD [ "/bin/bash" ]
