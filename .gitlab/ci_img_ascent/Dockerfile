#Docker Image for ascent

FROM ubuntu:20.04

LABEL maintainer="Deepak Rangarajan <deepak.rangarajan@netl.doe.gov>"

ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE 1
ENV TZ=America/New_York
ENV DEBIAN_FRONTEND=noninteractive

SHELL ["/bin/bash", "-o", "pipefail", "-c"]


RUN apt-get -qq update\
    && apt-get -qq upgrade -y\
    && apt-get -qq -y install --no-install-recommends \
    ca-certificates \
    curl=7.* \
    gnupg=2.* \
    && apt-get -qq update \
    && apt-get -qq -y install --no-install-recommends \
    build-essential=12.* \
    ccache=3.* \
    g++=4:9.* gfortran-9=9.* \
    git=1:2.* \
    libopenmpi-dev=4.* \
    openmpi-bin=4.* \
    pkg-config=0.29.* \
    python-is-python3=3.* \
    python3-pip=20.* \
    python3-setuptools=45.* \
    python3-venv=3.* \
    shellcheck=0.* \
    && apt-get clean && rm -rf /var/lib/apt/lists/* \
    && pip3 install --no-cache-dir \
    cmake==3.* \
    ninja==1.* \
    && curl -L https://github.com/mozilla/sccache/releases/download/v0.2.15/sccache-v0.2.15-x86_64-unknown-linux-musl.tar.gz -- \
    | tar xz \
    && cp sccache-v0.2.15-x86_64-unknown-linux-musl/sccache /usr/local/bin/sccache \
    && rm -rf sccache-v0.2.15-x86_64-unknown-linux-musl/ \
    && chmod +x /usr/local/bin/sccache

# Install conduit
RUN git clone --recursive --depth 1 --branch v0.8.0 https://github.com/LLNL/conduit.git
WORKDIR /conduit
RUN cmake -S src -B build \
	       -DCMAKE_BUILD_TYPE=Release \
	       -DENABLE_MPI=ON \
	       -DENABLE_OPENMP=OFF \
          -DENABLE_CUDA=OFF \
          -G Ninja
RUN cmake --build build --target install
WORKDIR /

# Install vtk-m
RUN git clone --branch v1.7.0 https://gitlab.kitware.com/vtk/vtk-m.git
WORKDIR /vtk-m
RUN git checkout d5143ad7ed2e55025cec2f49fb3f0ab46045a8e2
RUN cmake -S . -B build \
          -DCMAKE_BUILD_TYPE=Release \
          -DVTKm_ENABLE_MPI=ON \
          -DVTKm_ENABLE_OPENMP=OFF \
          -DVTKm_ENABLE_CUDA=OFF \
          -DVTKm_USE_64BIT_IDS=OFF \
          -DVTKm_USE_DOUBLE_PRECISION=ON \
          -DVTKm_USE_DEFAULT_TYPES_FOR_ASCENT=ON \
          -DVTKm_NO_DEPRECATED_VIRTUAL=ON \
          -G Ninja
RUN cmake --build build --target install
WORKDIR /

# Install vtk-h
RUN git clone --recursive --branch v0.8.0 https://github.com/Alpine-DAV/vtk-h.git
WORKDIR /vtk-h
RUN git checkout be86a73f6dc53176e076002f9520bc37b16b5119
RUN cmake -S src -B build \
          -DVTKM_DIR=/usr/local \
          -DENABLE_MPI=ON \
          -DENABLE_OPENMP=OFF \
          -DENABLE_CUDA=OFF \
          -G Ninja
RUN cmake --build build --target install
WORKDIR /

# Install ascent
RUN git clone --recursive --depth 1 --branch v0.8.0 https://github.com/Alpine-DAV/ascent.git
WORKDIR /ascent
RUN cmake -S src -B build \
          -DVTKM_DIR=/usr/local \
          -DVTKH_DIR=/usr/local \
          -DCONDUIT_DIR=/usr/local \
	  -DENABLE_FORTRAN=OFF \
	  -DENABLE_PYTHON=OFF \
	  -DENABLE_DOCS=OFF \
          -DBUILD_SHARED_LIBS=ON \
          -G Ninja
RUN cmake --build build --target install
WORKDIR /

RUN printf '#!/bin/bash\nenv $@' > /usr/local/bin/srun \
    && chmod +x /usr/local/bin/srun \
    && useradd --create-home -s /bin/bash user

WORKDIR /home/user
USER user

CMD [ "/bin/bash" ]
