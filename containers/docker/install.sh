#!/usr/bin/env bash

apt-get update
apt-get install -y \
    build-essential \
    libbz2-dev \
    liblzma-dev \
    cmake \
    git \
    python3-pip \
    autoconf \
    sudo \
    wget

sudo pip3 install git+https://github.com/iqbal-lab-org/gramtools

apt-get autoremove --purge -y \
    cmake \
    git \
    autoconf \
    wget

apt-get clean

export SUDO_FORCE_REMOVE=yes
apt-get autoremove --purge -y sudo
