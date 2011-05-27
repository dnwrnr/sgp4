#include "Globals.h"

void TrimLeft(std::string& str) {

    std::string whitespaces(" \t\f\n\r");

    if (!str.empty()) {
        std::string::size_type pos = str.find_first_not_of(whitespaces);

        if (pos != std::string::npos)
            str.erase(0, pos);
        else
            str.clear();
    }
}

void TrimRight(std::string& str) {

    std::string whitespaces(" \t\f\n\r");

    if (!str.empty()) {
        std::string::size_type pos = str.find_last_not_of(whitespaces);

        if (pos != std::string::npos)
            str.erase(pos + 1);
        else
            str.clear();
    }
}

void Trim(std::string& str) {

    TrimLeft(str);
    TrimRight(str);
}

