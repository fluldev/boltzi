#pragma once
#include <memory>
#include <sstream>
#include <iostream>

namespace boltzi {
namespace config {
class Silencer {  // Singleton
    private:
    static std::unique_ptr<Silencer> instance;    
    bool no_stdout = false;
    bool less_stdout = false;
    public:
    static bool silent() {
        if(!instance)
            return false;
        return instance->no_stdout;
    }
    static bool less_output() {
        if(!instance)
            return false;
        return instance->less_stdout;
    }

    static void set(bool setting) {
        if(!instance)
            instance = std::unique_ptr<Silencer>(new Silencer);
        instance->no_stdout = setting;
    }

    static void set_less(bool setting) {
        if(!instance)
            instance = std::unique_ptr<Silencer>(new Silencer);
        instance->less_stdout = setting;
    }

    static void print(const std::stringstream& output) {
        if(!silent())
            std::cout<<output.str();
    }
    static void print(const std::string& output) {
        if(!silent())
            std::cout<<output;
    }
};
}
}
