set(BOOST_SRC_DIR ${CMAKE_CURRENT_BINARY_DIR}/src/boost_1_65_1)

execute_process(COMMAND mkdir -p ${BOOST_SRC_DIR})

ExternalProject_Add(boost
        DOWNLOAD_COMMAND  wget https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.gz --timestamping
        DOWNLOAD_DIR      "${CMAKE_CURRENT_BINARY_DIR}/download"
        CONFIGURE_COMMAND bash -c "cd ${BOOST_SRC_DIR} && ./bootstrap.sh \
                        --with-libraries=program_options,timer,system,serialization \
                        --prefix=${BOOST_SRC_DIR} 1>/dev/null"
        BUILD_COMMAND     bash -c "cd ${BOOST_SRC_DIR} && ./bjam install link=static 1>/dev/null"
        INSTALL_COMMAND   ""
        TEST_COMMAND      "")

ExternalProject_Add_Step(boost extract_tar
        COMMAND tar -xzf ${CMAKE_CURRENT_BINARY_DIR}/download/boost_1_65_1.tar.gz -C ${CMAKE_CURRENT_BINARY_DIR}/src
        DEPENDEES download
        DEPENDERS configure)
