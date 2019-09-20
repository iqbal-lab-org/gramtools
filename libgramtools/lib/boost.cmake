execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/src/boost_1_65_1)

ExternalProject_Add(boost
        DOWNLOAD_COMMAND  wget https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.gz --timestamping
        DOWNLOAD_DIR      "${CMAKE_CURRENT_BINARY_DIR}/download"
        CONFIGURE_COMMAND bash -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src/boost_1_65_1 && ./bootstrap.sh --with-libraries=random,program_options,timer,system,filesystem,serialization --prefix=${CMAKE_CURRENT_BINARY_DIR}"
        BUILD_COMMAND     bash -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src/boost_1_65_1 && ./bjam install link=static"
        INSTALL_COMMAND   ""
        TEST_COMMAND      "")

ExternalProject_Add_Step(boost extract_tar
        COMMAND tar -xvzf ${CMAKE_CURRENT_BINARY_DIR}/download/boost_1_65_1.tar.gz -C ${CMAKE_CURRENT_BINARY_DIR}/src
        DEPENDEES download
        DEPENDERS configure)