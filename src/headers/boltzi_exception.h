#pragma once

#include<exception>

namespace boltzi {
    class BoltziException : public std::exception {
        private:
        const char* str; 
        public:

        // this might refer to non existant memory
        const char* what() const noexcept override {return str;}
        BoltziException(const char* str) : str(str) {}
    }; 
}
