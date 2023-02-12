#pragma once

#include <cmath>
#include<random>
#include<functional>

#include <iostream>
#include <sstream>

#include "../boltzi_types.h"
#include "../boltzi_exception.h"
#include "constants.h"
#include "basics.h"
#include "containers.h"
#include "../silencer.h"



namespace boltzi{
using config::Silencer;
namespace util{
    // random utilities should be thread_local!
    class RandomUtil {
        private:  
        std::default_random_engine eng;
        std::uniform_real_distribution<double> rand_real;

        public:
        RandomUtil() { 
            std::random_device r(std::string("/dev/urandom"));
            eng = std::default_random_engine(r());
            rand_real = std::uniform_real_distribution<double>(0, 1);
        }

        double draw_real() {
            return rand_real(eng);
        }


        unsigned int draw_uint(unsigned max_val) {
            return std::uniform_int_distribution<unsigned>(0, max_val)(eng);
        }

        template<typename T>
        double draw_element(const T& container) {
            std::uniform_int_distribution<unsigned int> rand_int(0, container.size()-1);
            return container[rand_int(eng)];
        }
        
        template<typename T>
        unsigned int draw_index_weighted(const T& weights) {
            double weight = std::accumulate(weights.cbegin(), weights.cend(), 0.) 
                                * rand_real(eng);
            double prev_sum = 0;
            for(auto it = weights.cbegin(); it!= weights.cend(); ++it) {
                prev_sum+=*it;
                if(weight<=prev_sum) return it-weights.begin();
            }
            throw BoltziException("Did not find matching bin in RandomUtil::draw_index_weighted.");
        }

        inline vec3 random_dir() {
            return sphere_to_cart({
                2*PI*draw_real(),
                std::acos(1-2*draw_real())
            });
        }

        inline vec3 random_dir_area(vec3 normal, vec3 x, vec3 y) {
            vec3 cart = sphere_to_cart({
                2*PI*draw_real(),
                std::asin(std::sqrt(draw_real()))
            });
            return cart[2] * normal + cart[0] * x + cart[1] * y;
        }

        inline double accept_reject(std::function<double(double)> fn, double a, double b, double max) {
            for(unsigned i=0;i<config::MAX_ACCEPT_REJECT_ATTEMPTS;++i) {
                double x = draw_real()*(b-a) + a;
                double y = draw_real()*max;
                if(y<=fn(x)) return x;
            } 
            Silencer::print(std::stringstream{}<<"Accept Reject Problem: "<<max<<", "<<fn(b)<<", "<<fn((a+b)/2)<<"\n");
            return draw_real()*(b-a) + a;
        }

        inline vec3 random_pos_cell(vec3 r0, vec3 cell_size) {
            return r0+vec3{
                draw_real()*cell_size[0], 
                draw_real()*cell_size[1], 
                draw_real()*cell_size[2], 
            };
        }
        
        inline vec3 random_pos_triangle(const Triangle& tri) {
            double z_1 = 1-std::sqrt(draw_real());
            double z_2 = draw_real();
            return  z_1*(tri.points[1]-tri.points[0])
                    + z_2*(1-z_1)*(tri.points[2]-tri.points[0])
                    + tri.points[0];
        }


        inline vec3 random_pos_tetrahedron(const Tetrahedron& tetra) {
            double z_1 = 1-std::pow(draw_real(),1./3);
            double z_2 = 1-std::sqrt(draw_real());
            double z_3 = draw_real();
            return  z_1*(tetra.points[1]-tetra.points[0])
                    + z_2*(1-z_1)*(tetra.points[2]-tetra.points[0])
                    + z_3*(1-z_2)*(1-z_1)*(tetra.points[3]-tetra.points[0])
                    + tetra.points[0];
        }
    };

    // random number generator utility class shared on thread level
    extern thread_local RandomUtil random;


    class BinDist {
        private:
        double a,b;
        unsigned int n_bins;
        std::vector<double> bin_edges;

        public:
        BinDist(std::function<double(double)> fn, double a, double b, unsigned int n_bins) {
            this->n_bins = n_bins;
            this->a = a;
            this->b = b;
            bin_edges.reserve(n_bins+1);
            bin_edges.push_back(0);
            for(unsigned i = 1; i<=n_bins; ++i)
                bin_edges.push_back(
                        fn((i-.5)/n_bins*(b-a) + a)
                        + bin_edges[i-1]
                );
        }

        double operator() () {
            double weight = random.draw_real() * bin_edges.back();
            for(auto it = bin_edges.begin(); it!=bin_edges.end(); ++it)
                if(*it >= weight) return (it-bin_edges.begin() - 1 + random.draw_real())*(b-a)/n_bins+a;
            throw BoltziException("No matching bin in BinDist.");
        }
    };
}}
