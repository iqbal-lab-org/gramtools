#!/usr/bin/env bash
set -e

# Notes
# * Dependencies:
#       * cortex: used in variant discovery. depends on r-base

apt-get update
apt-get install -y \
    build-essential \
    cmake \
	automake \
	git \
	liblzma-dev \
	libbz2-dev \
	zlib1g-dev \
	pkg-config \
	python3 \
	python3-pip \
	python3-setuptools \
	python \
	python-dev \
	r-base \
	wget

pip3 install pip==20.0.2 # upgrade pip
