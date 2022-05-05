# Docker Image for running ascent test

FROM alpinedav/ascent:0.8.0

LABEL maintainer="Deepak Rangarajan <deepak.rangarajan@netl.doe.gov>"

ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE 1
ENV TZ=America/New_York
ENV DEBIAN_FRONTEND=noninteractive

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get -qq update\
    && apt-get -qq upgrade -y\
    && apt-get -qq -y install --no-install-recommends \
    ccache=3.* \
    ninja-build=1.* \
    libopenmpi-dev=2.* \
    openmpi-bin=2.* \
    && apt-get clean && rm -rf /var/lib/apt/lists/* \
    && curl -L https://github.com/mozilla/sccache/releases/download/v0.2.15/sccache-v0.2.15-x86_64-unknown-linux-musl.tar.gz -- \
    | tar xz \
    && cp sccache-v0.2.15-x86_64-unknown-linux-musl/sccache /usr/local/bin/sccache \
    && rm -rf sccache-v0.2.15-x86_64-unknown-linux-musl/ \
    && chmod +x /usr/local/bin/sccache

RUN printf '#!/bin/bash\nenv $@' > /usr/local/bin/srun \
    && chmod +x /usr/local/bin/srun \
    && useradd --create-home -s /bin/bash user

WORKDIR /home/user
USER user

CMD [ "/bin/bash" ]
