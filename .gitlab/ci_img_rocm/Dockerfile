# Docker Image for running MfiX-Exa ctests

FROM rocm/dev-ubuntu-20.04

LABEL maintainer="Justin Weber <justin.weber@netl.doe.gov>"

ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE 1
ENV DEBIAN_FRONTEND=noninteractive

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get -qq update\
    && apt-get -qq upgrade -y\
    && apt-get -qq -y install --no-install-recommends \
    build-essential=12.* \
    ccache=3.* \
    git=1:2.* \
    libopenmpi-dev=4.* \
    openmpi-bin=4.* \
    pkg-config=0.29.* \
    python-is-python3=3.* \
    python3-pip=20.* \
    python3-setuptools=45.* \
    python3-venv=3.* \
    && apt-get clean && rm -rf /var/lib/apt/lists/* \
    && pip3 install --no-cache-dir \
    cmake==3.* \
    conan==1.* \
    ninja==1.* \

RUN printf '#!/bin/bash\nenv $@' > /usr/local/bin/srun \
    && chmod +x /usr/local/bin/srun \
    && useradd --create-home -s /bin/bash user

WORKDIR /home/user
USER user

CMD [ "/bin/bash" ]
