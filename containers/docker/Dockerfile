FROM ubuntu:20.04

# INSTALL
COPY . /gramtools
COPY ./containers/docker/install.sh /gramtools/
WORKDIR /gramtools
# to avoid tzdata user prompt
ENV DEBIAN_FRONTEND=noninteractive
RUN chmod +x ./install.sh && ./install.sh

# CLEANUP
WORKDIR /
RUN rm -rf ./gramtools

ENTRYPOINT ["gramtools"]
CMD []
