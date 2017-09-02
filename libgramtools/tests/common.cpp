#include <vector>
#include <string>


std::vector<uint8_t> encode_read(const std::string &read) {
    std::vector<uint8_t> encoded_read;
    for (uint16_t i = 0; i < read.length(); i++) {
        if (read[i] == 'A' or read[i] == 'a') encoded_read.push_back(1);
        if (read[i] == 'C' or read[i] == 'c') encoded_read.push_back(2);
        if (read[i] == 'G' or read[i] == 'g') encoded_read.push_back(3);
        if (read[i] == 'T' or read[i] == 't') encoded_read.push_back(4);
    }
    return encoded_read;
}
