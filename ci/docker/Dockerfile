FROM ubuntu:20.04

COPY . /gramtools
COPY ./ci/docker/purge.sh /gramtools/
COPY ./ci/install_deps.sh /gramtools/

# INSTALL
WORKDIR /gramtools
# to avoid tzdata user prompt
ENV DEBIAN_FRONTEND=noninteractive
RUN chmod +x ./install_deps.sh && ./install_deps.sh
RUN pip3 install .

# CLEANUP
WORKDIR /
RUN chmod +x ./gramtools/purge.sh && ./gramtools/purge.sh
RUN rm -rf ./gramtools

ENTRYPOINT ["gramtools"]
CMD []
