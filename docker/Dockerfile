FROM ubuntu:latest
MAINTAINER Téo Lemane

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y dist-upgrade \
    && apt-get install -y --no-install-recommends wget

RUN apt-get install -y git cmake gcc g++ zlib1g zlib1g-dev

RUN cd /opt \
    && mkdir kmtricks && cd kmtricks \
    && wget https://github.com/tlemane/kmtricks/releases/download/v1.3.0/kmtricks-v1.3.0-sources.tar.gz \
    && tar -zxvf kmtricks-v1.3.0-sources.tar.gz

RUN cd /opt \
    && cd kmtricks \
    && mkdir build \
    && cd build \
    && cmake .. -DKMER_LIST="32 64 96 128" -DWITH_MODULES=ON -DWITH_HOWDE=ON -DWITH_SOCKS=ON -DNATIVE=OFF \
    && make -j8 \
    && chmod +x ../bin/kmtricks \
    && cp ../bin/kmtricks /bin

RUN cd /opt \
    && cd kmtricks/build \
    && cmake .. -DWITH_PLUGIN=ON \
    && make -j8 \
    && chmod +x ../bin/kmtricks \
    && cp ../bin/kmtricks /bin/kmtricksp

RUN rm -rf /opt/kmtricks

WORKDIR /tmp

ENTRYPOINT ["kmtricks"]
