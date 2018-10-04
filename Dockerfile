ARG UBUNTU_VERSION=16

FROM ubuntu:16.04 AS base-16
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -y install \
    bzip2 \
    default-jre \
    git \ 
    perl-base \
    php7.0-cli \ 
    php7.0-xml \ 
    php7.0-mysql \
    python-matplotlib \ 
    python-numpy \
    python-pysam \
    r-base-core \ 
    r-cran-optparse \ 
    r-cran-robustbase \ 
    r-cran-foreach \ 
    r-cran-doparallel \ 
    r-cran-mass \
    tabix \
    unzip \
    wget

FROM ubuntu:18.04 AS base-18
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -y install \
    bzip2 \
    default-jre \
    git \ 
    perl-base \ 
    php7.2-cli \ 
    php7.2-xml \ 
    php7.2-mysql \ 
    python-matplotlib \
    python-numpy \
    python-pysam \
    r-base-core \ 
    r-cran-optparse \ 
    r-cran-robustbase \ 
    r-cran-foreach \ 
    r-cran-doparallel \ 
    r-cran-mass \
    tabix \ 
    unzip \ 
    wget

FROM base-${UBUNTU_VERSION} AS tools-ubuntu-16
RUN apt-get update &&  DEBIAN_FRONTEND=noninteractive apt-get -y install \
    build-essential \ 
    cmake \ 
    cpanminus \
    libbz2-dev \ 
    liblzma-dev \ 
    libncurses5-dev \ 
    libpng-dev \ 
    libmysqlclient-dev \
    libqt5sql5-mysql \ 
    libqt5xmlpatterns5-dev \ 
    libssl-dev \ 
    qt5-default \ 
    qt5-qmake \ 
    qtbase5-dev 

FROM base-${UBUNTU_VERSION} AS tools-ubuntu-18
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -y install \
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
RUN chmod 755 download_*.sh && ./download_tools.sh

FROM base-${UBUNTU_VERSION}
COPY --from=build /megSAP/ /megSAP/
COPY --from=build /megSAP/data/dbs/ /megSAP/data/dbs_static/
COPY --from=build /usr/local/share/perl/ /usr/local/share/perl/
COPY --from=build /usr/local/lib/x86_64-linux-gnu/perl/ /usr/local/lib/x86_64-linux-gnu/perl/
COPY --from=build /usr/local/bin/ /usr/local/bin/

WORKDIR /megSAP
ENTRYPOINT ["/bin/bash"]
