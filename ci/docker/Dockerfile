FROM ubuntu:20.04

# INSTALL
COPY . /gramtools
COPY ./ci/docker/purge.sh /gramtools/
COPY ./ci/install_deps.sh /gramtools/
WORKDIR /gramtools
# to avoid tzdata user prompt
ENV DEBIAN_FRONTEND=noninteractive
RUN chmod +x ./install_deps.sh && ./install_deps.sh
RUN pip3 install .
RUN chmod +x ./purge.sh && ./purge.sh

# CLEANUP
WORKDIR /
RUN rm -rf ./gramtools

ENTRYPOINT ["gramtools"]
CMD []
