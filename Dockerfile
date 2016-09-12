# modified from https://github.com/ContinuumIO/docker-images/blob/master/anaconda/Dockerfile

FROM ubuntu:16.04

MAINTAINER Mark Meredith <mark.meredith@netl.doe.gov>

RUN apt-get update --fix-missing && apt-get install -y wget ca-certificates \
    git gfortran cmake make openmpi-bin libopenmpi-dev bzip2

CMD [ "/bin/bash" ]
