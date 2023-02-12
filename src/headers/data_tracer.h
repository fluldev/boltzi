#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <mutex>

namespace boltzi{
template<typename... T>
class DataTracer {
    public:
        using row_type = std::tuple<T...>;
    private:
        std::vector<row_type> data;
        std::string filename;
        std::mutex lock;

        template<size_t i>
        void append_to_stream(std::ofstream& s, const row_type& t) {
            s << std::get<i>(t);
            if constexpr(i+1 != std::tuple_size<row_type>::value)
                s << ",";
            if constexpr(i+2 < std::tuple_size<row_type>::value)
                append_to_stream<i+1>(s, t);  
            else if constexpr(i+1 < std::tuple_size<row_type>::value)
                s << std::get<i+1>(t);
        }


    public:
        DataTracer(std::string filename) {this->filename = filename;}

        void save() {
            std::lock_guard<std::mutex> guard(lock);
            std::ofstream file(filename);
            for(auto row : data) {
                append_to_stream<0>(file, row);
                file << "\n";
            }
            file.close();
        }

        void trace(row_type&& d) {
            std::lock_guard<std::mutex> guard(lock);
            data.push_back(std::move(d));
        }
        
        void trace(const row_type& d) {
            std::lock_guard<std::mutex> guard(lock);
            data.push_back(d);
        }
};
}
