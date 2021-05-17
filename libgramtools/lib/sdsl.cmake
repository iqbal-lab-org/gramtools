execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/src/sdsl-lite-2.1.1)

ExternalProject_Add(sdsl
        DOWNLOAD_COMMAND  wget https://github.com/simongog/sdsl-lite/releases/download/v2.1.1/sdsl-lite-2.1.1.tar.gz.offline.install.gz --timestamping
        DOWNLOAD_DIR      "${CMAKE_CURRENT_BINARY_DIR}/download"
        CONFIGURE_COMMAND ""
        BUILD_COMMAND     bash -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src/sdsl-lite-2.1.1 && ./install.sh ${CMAKE_CURRENT_BINARY_DIR} && rm -rf ${CMAKE_CURRENT_BINARY_DIR}/include/{gtest,gmock}"
        INSTALL_COMMAND   ""
        TEST_COMMAND      "")

ExternalProject_Add_Step(sdsl extract_tar
        COMMAND tar -xzf ${CMAKE_CURRENT_BINARY_DIR}/download/sdsl-lite-2.1.1.tar.gz.offline.install.gz -C ${CMAKE_CURRENT_BINARY_DIR}/src
        DEPENDEES download
        DEPENDERS configure)
