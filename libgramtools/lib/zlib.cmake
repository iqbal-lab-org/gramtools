execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/src/zlib-1.2.11)

ExternalProject_Add(zlib
        DOWNLOAD_COMMAND  wget https://github.com/madler/zlib/archive/v1.2.11.tar.gz --timestamping -O zlib.tar.gz
        DOWNLOAD_DIR      "${CMAKE_CURRENT_BINARY_DIR}/download"
        CONFIGURE_COMMAND bash -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src/zlib-1.2.11 && ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR}"
        BUILD_COMMAND     bash -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src/zlib-1.2.11 && make install"
        INSTALL_COMMAND   ""
        TEST_COMMAND      "")

ExternalProject_Add_Step(zlib extract_tar
        COMMAND tar -xzf ${CMAKE_CURRENT_BINARY_DIR}/download/zlib.tar.gz -C ${CMAKE_CURRENT_BINARY_DIR}/src
        DEPENDEES download
        DEPENDERS configure)