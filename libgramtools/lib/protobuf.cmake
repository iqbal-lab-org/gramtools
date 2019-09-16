execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/src/protobuf-3.9.1)

ExternalProject_Add(protobuf
        DOWNLOAD_COMMAND  wget https://github.com/protocolbuffers/protobuf/releases/download/v3.9.1/protobuf-cpp-3.9.1.tar.gz --timestamping
        DOWNLOAD_DIR      "${CMAKE_CURRENT_BINARY_DIR}/download"
        CONFIGURE_COMMAND bash -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src/protobuf-3.9.1 && ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR}"
        BUILD_COMMAND     bash -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src/protobuf-3.9.1 && make && make check && make install"
        INSTALL_COMMAND   ""
        TEST_COMMAND      "")

ExternalProject_Add_Step(protobuf extract_tar
        COMMAND tar -xzf ${CMAKE_CURRENT_BINARY_DIR}/download/protobuf-cpp-3.9.1.tar.gz -C ${CMAKE_CURRENT_BINARY_DIR}/src
        DEPENDEES download
        DEPENDERS configure)