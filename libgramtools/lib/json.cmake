ExternalProject_Add(nlohmann_json
        DOWNLOAD_COMMAND  wget https://github.com/nlohmann/json/releases/download/v3.7.3/json.hpp --timestamping
        DOWNLOAD_DIR      "${CMAKE_CURRENT_BINARY_DIR}/include/nlohmann"
        CONFIGURE_COMMAND ""
        BUILD_COMMAND     ""
        INSTALL_COMMAND   ""
        TEST_COMMAND      "")

