#pragma once

#include <initializer_list>
#include <iterator>
#include <utility>
#include <vector>
#include <array>
#include <numeric> 
#include <cassert>
#include <armadillo>
#include <cmath>

#include <iostream>

#include "../boltzi_types.h"
#include "../boltzi_config.h"


namespace boltzi {
namespace util{
    template<unsigned int n_dims, typename T>
    class Grid {
        private:
        std::vector<T> data;
        std::array<unsigned int, n_dims> index_functional;

        inline unsigned int lin_index(std::initializer_list<unsigned int> nd_index) const {
            return std::inner_product(
                    nd_index.begin(), nd_index.end(), 
                    index_functional.begin(), 0u); 
        }
        
        public:
        class iterator {
            private:
            typename std::vector<T>::iterator v_it;
            const Grid* owner;
            public:
            iterator() = delete;
            iterator(const Grid* owner, typename std::vector<T>::iterator beg) :  v_it(beg), owner(owner) {}

            iterator operator++() {
                ++v_it; 
                return *this;
            }
            
            iterator operator--() {
                --v_it; 
                return *this;
            }

            iterator operator+(const unsigned b) const {
                return iterator{owner, v_it+b};
            }

            int operator-(const iterator& b) const {
                return v_it-b.v_it;
            }

            bool operator!=(const iterator& b) const {
                return v_it != b.v_it;
            }
            bool operator==(const iterator& b) const {
                return v_it == b.v_it;
            }

            std::pair<std::array<unsigned,n_dims>, T&> operator*() {
                std::array<unsigned,n_dims> ret_idx;
                unsigned cnt = (v_it - owner->data.begin());
                for(int i = n_dims-1; i>=0; --i) {
                    ret_idx[i] =  cnt / owner->index_functional[i];
                    cnt %= owner->index_functional[i];
                }
                return std::pair<std::array<unsigned,n_dims>, T&>{ret_idx, *v_it};
            }

            std::pair<std::array<unsigned,n_dims>, const T&> operator*() const {
                std::array<unsigned,n_dims> ret_idx;
                unsigned cnt = (v_it - owner->data.begin());
                for(int i = n_dims-1; i>=0; --i) {
                    ret_idx[i] =  cnt / owner->index_functional[i];
                    cnt %= owner->index_functional[i];
                }
                return std::pair<std::array<unsigned,n_dims>, const T&>{ret_idx, *v_it};
            }
        };

        std::array<unsigned int, n_dims> shape;
        const unsigned int dims = n_dims;

        const std::vector<T>& get_data() const {return data;}

        Grid() = delete;
        Grid(std::initializer_list<unsigned int> shape) {
            assert(n_dims == shape.size());
            std::copy(shape.begin(), shape.end(), this->shape.begin());

            // could be deduced from the loop later in this constructor
            unsigned int length = std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<>{});
            data.reserve(length);
            for(unsigned i=0; i<length;++i)  // can't use fill_n since it doesn't work with atomics
                data.push_back(T());

            index_functional[0] = 1;
            for(auto i = 1u; i<n_dims; ++i)
                index_functional[i] = index_functional[i-1] * this->shape[i-1];
        }

        T& at(std::initializer_list<unsigned int> nd_index) {
            return data.at(
                lin_index(nd_index)
            ); 
        }
        
        const T& at(std::initializer_list<unsigned int> nd_index) const {
            return data.at(
                lin_index(nd_index)
            ); 
        }

        iterator begin() {return iterator{this, data.begin()};}
        iterator end() {return iterator{this, data.end()};}
        size_t size() const {return data.size();}
    };
    //
    //TODO: center positions of Triangle and Tetrahedron might be wrong
    class Tetrahedron {
        public:
        vec3 points[4];
        double volume;
        vec3 center;

        Tetrahedron() = delete;

        Tetrahedron(vec3 p1, vec3 p2, vec3 p3, vec3 p4) {
            points[0] = p1;
            points[1] = p2;
            points[2] = p3;
            points[3] = p4;
            volume = std::abs(arma::dot(arma::cross(points[1]-points[0], points[2]-points[0]), points[3]-points[0]))/6;
            center = (points[1]-points[0])/2 + (points[2]-points[0])/4 + (points[3]-points[0])/8 + points[0];
        }

        inline bool is_inside(const vec3& r) const {
            double volume_1 = std::abs(arma::dot(arma::cross(points[1]-points[0], points[2]-points[0]), r-points[0]))/6;
            double volume_2 = std::abs(arma::dot(arma::cross(points[1]-r, points[2]-r), points[3]-r))/6;
            double volume_3 = std::abs(arma::dot(arma::cross(r-points[0], points[2]-points[0]), points[3]-points[0]))/6;
            double volume_4 = std::abs(arma::dot(arma::cross(points[1]-points[0], r-points[0]), points[3]-points[0]))/6;
            return (volume_1 + volume_2 + volume_3 + volume_4) <= volume*config::INSIDE_CHECK_VOLUME_MARGIN;
        }
    };

    class Triangle {
        public:
        vec3 points[3];
        double area;
        vec3 center;

        Triangle() = delete;

        Triangle(vec3 p1, vec3 p2, vec3 p3) {
            points[0] = p1;
            points[1] = p2;
            points[2] = p3;
            area = arma::norm(arma::cross(points[1]-points[0], points[2]-points[0]))/2;
            center = (points[1]-points[0])/2 + (points[2]-points[0])/4 + points[0];
        }

        inline bool is_inside(const vec3& r) const {
            double area_1 = arma::norm(arma::cross(r-points[0], r-points[1]))/2;
            double area_2 = arma::norm(arma::cross(r-points[1], r-points[2]))/2;
            double area_3 = arma::norm(arma::cross(r-points[0], r-points[2]))/2;
            return (area_1 + area_2 + area_3) <= area*config::INSIDE_CHECK_AREA_MARGIN;
        }
    };
}}
