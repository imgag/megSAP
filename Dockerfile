ARG UBUNTU_VERSION=16

FROM ubuntu:16.04 AS base-16
RUN DEBIAN_FRONTEND=noninteractive apt-get update && apt-get -y install \
    bzip2 \
    default-jre \
    git \ 
    perl-base \
    php7.0-cli \ 
    php7.0-xml \ 
    php7.0-mysql \
    python-matplotlib \ 
    tabix \
    unzip \
    wget

FROM ubuntu:18.04 AS base-18
RUN DEBIAN_FRONTEND=noninteractive apt-get update && apt-get -y install \
    bzip2 \
    default-jre \
    git \ 
    perl-base \ 
    php7.2-cli \ 
    php7.2-xml \ 
    php7.2-mysql \ 
    python-matplotlib \ 
    tabix \ 
    unzip \ 
    wget

FROM base-${UBUNTU_VERSION} AS tools-ubuntu-16
RUN DEBIAN_FRONTEND=noninteractive apt-get update && apt-get -y install \
    build-essential \ 
    cmake \ 
    cpanminus \
    libbz2-dev \ 
    liblzma-dev \ 
    libncurses5-dev \ 
    libpng-dev \ 
    libqt5sql5-mysql \ 
    libqt5xmlpatterns5-dev \ 
    libssl-dev \ 
    mysql-client \
    qt5-default \ 
    qt5-qmake \ 
    qtbase5-dev 

FROM base-${UBUNTU_VERSION} AS tools-ubuntu-18
RUN DEBIAN_FRONTEND=noninteractive apt-get update && apt-get -y install \
    build-essential \ 
    cmake \ 
    cpanminus \
    libbz2-dev \ 
    liblzma-dev \ 
    libmysqlclient-dev \
    libncurses5-dev \ 
    libqt5sql5-mysql \ 
    libpng-dev \
    libqt5xmlpatterns5-dev \ 
    libssl-dev \
    qt5-default \ 
    qt5-qmake \ 
    qtbase5-dev

FROM tools-ubuntu-${UBUNTU_VERSION} AS build
ADD . /megSAP
WORKDIR /megSAP/data
RUN chmod 755 download_*.sh && ./download_tools.sh && ./download_tools_vep.sh

FROM base-${UBUNTU_VERSION}
RUN useradd -d /home/ubuntu -ms /bin/bash -g root -p ubuntu ubuntu
WORKDIR /home/ubuntu
COPY --from=build --chown=ubuntu:nogroup /megSAP/ /home/ubuntu/megSAP/
COPY --from=build --chown=ubuntu:nogroup /megSAP/data/dbs /home/ubuntu/megSAP/data/dbs_static
COPY --from=build --chown=ubuntu:nogroup /root/.cpanm/ /home/ubuntu/.cpanm/
COPY --from=build /usr/local/share/perl/ /usr/local/share/perl/
COPY --from=build /usr/local/lib/x86_64-linux-gnu/perl/ /usr/local/lib/x86_64-linux-gnu/perl/
COPY --from=build /usr/local/bin/ /usr/local/bin/

WORKDIR /home/ubuntu/megSAP
ENTRYPOINT ["/bin/sh"]
