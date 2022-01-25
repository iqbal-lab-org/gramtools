ExternalProject_Add(htslib
        DOWNLOAD_COMMAND  wget https://github.com/samtools/htslib/releases/download/1.10/htslib-1.10.tar.bz2 --timestamping
        DOWNLOAD_DIR      "${CMAKE_CURRENT_BINARY_DIR}/download"
        CONFIGURE_COMMAND bash -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src/htslib-1.10 && \
            ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR} \
            --disable-bz2 --disable-lzma --enable-libcurl=no"
        BUILD_COMMAND     bash -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src/htslib-1.10 && make && make install"
        INSTALL_COMMAND   ""
        TEST_COMMAND      "")

ExternalProject_Add_Step(htslib extract_tar
        COMMAND tar -xjf ${CMAKE_CURRENT_BINARY_DIR}/download/htslib-1.10.tar.bz2 -C ${CMAKE_CURRENT_BINARY_DIR}/src
        DEPENDEES download
        DEPENDERS configure)
