#!/usr/bin/env bash
set -e


apt-get autoremove --purge -y \
    build-essential \
    cmake \
    git \
    autoconf \
    wget

apt-get clean

export SUDO_FORCE_REMOVE=yes
apt-get autoremove --purge -y sudo

conan remove -f "*"
pip3 uninstall -y conan
rm -rf ./root/.cache
