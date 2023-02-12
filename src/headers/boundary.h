#pragma once

#include <cstddef>
#include <vector>
#include <memory>
#include <optional>
#include <tuple>
#include <armadillo>
#include "boltzi_types.h" 
#include "boltzi_exception.h"
#include "boltzi_util.h"
#include "phonon.h"
#include "util/containers.h"

namespace boltzi {
    class Boundary : public std::enable_shared_from_this<Boundary> {
        // TODO: this should be only in an inherited Boundary VertexBoundary to leave open the
        // possiblity of raymarching boundaries
        private:
        static double constant_weight(vec3) {return 1;}
        public:
        struct Primitive {
            util::Triangle tri;
            vec3 normal;
            vec3 e_1;
            vec3 e_2;
        };
        using primitive_iterator = std::vector<Primitive>::iterator;
        // This class depends on the type of Boundary. Classes which implement Boundary must implement there own CollisionEvent
        // since it should be used by each boundary to query further information about a collision.
        struct CollisionEvent{
            size_t boundary_idx;
            vec3 position;
            std::shared_ptr<const Boundary> boundary;

            // this is a special property only belongig to vertex boundaries.
            unsigned primitive;
            bool negative;
        };

        std::vector<Primitive> primitives;

        Boundary() = delete;
        Boundary(const std::vector<Primitive>& primitives) : primitives(primitives) {}

        vec3 lowest_corner() const {
            vec3 res = primitives.front().tri.points[0];
            for(auto& pr : primitives) { 
                res[0] = std::min({res[0], pr.tri.points[0][0], pr.tri.points[1][0], pr.tri.points[2][0]});
                res[1] = std::min({res[1], pr.tri.points[0][1], pr.tri.points[1][1], pr.tri.points[2][1]});
                res[2] = std::min({res[2], pr.tri.points[0][2], pr.tri.points[1][2], pr.tri.points[2][2]});
            }
            return res;
        }

        vec3 highest_corner() const {
            vec3 res = primitives.front().tri.points[0];
            for(auto& pr : primitives) { 
                res[0] = std::max({res[0], pr.tri.points[0][0], pr.tri.points[1][0], pr.tri.points[2][0]});
                res[1] = std::max({res[1], pr.tri.points[0][1], pr.tri.points[1][1], pr.tri.points[2][1]});
                res[2] = std::max({res[2], pr.tri.points[0][2], pr.tri.points[1][2], pr.tri.points[2][2]});
            }
            return res;
        }
        vec3 reflect(const PhononData& phonon, std::shared_ptr<Boundary::CollisionEvent> coll) const {
            // TODO: reflects k directly which could change the phonon energy in an non isotropic material!
            return phonon.k - 2*arma::dot(phonon.k, primitives[coll->primitive].normal) * primitives[coll->primitive].normal;
        }

        std::shared_ptr<CollisionEvent> collision_point(vec3 orig, vec3 dir_vec, std::shared_ptr<CollisionEvent> last) const;
        std::vector<std::tuple<vec3,vec3,vec3,vec3>> random_positions(unsigned n, std::function<double(vec3)> weights=constant_weight) const;  // needed since the phonon propertys are decided 
        
        virtual double area() const {
            double sum = 0;
            for(const auto& tri : primitives)
                sum+=tri.tri.area;
            return sum;
        }
        
        virtual double surface_integral(std::function<double(vec3)> fn) {
            double res = 0;
            for(const auto& prim : primitives)
                res+=fn(prim.tri.center)*prim.tri.area;
            return res;

        }
                                                              // within sample since only sample knows the material at the spawn point
        virtual bool is_active() const = 0;
        virtual bool is_generator() const = 0;
        virtual double p_absorb(const PhononData& phonon, std::shared_ptr<Boundary::CollisionEvent> coll) const = 0;
        virtual double p_spec(const PhononData& phonon, std::shared_ptr<Boundary::CollisionEvent> coll) const = 0;
        virtual vec3 diffusive_draw(const PhononData& phonon, std::shared_ptr<Boundary::CollisionEvent> coll) const = 0;
        virtual double surface_temperature(vec3 pos) const = 0;
    };

    class SampleGeometry {
        private:
        std::vector<util::Tetrahedron> tetrahedra;

        static double constant_weight(vec3) {return 1;}

        public:
        double scale;
        std::vector<std::shared_ptr<Boundary>> boundaries;
        std::pair<vec3, vec3> box;  // first is min second is max

        bool is_inside_box(const vec3& pos) const {
            return (box.first[0] <= pos[0]) && (pos[0] < box.second[0])
                    && (box.first[1] <= pos[1]) && (pos[1] < box.second[1])
                    && (box.first[2] <= pos[2]) && (pos[2] < box.second[2]);
        }


        SampleGeometry() = delete;
        SampleGeometry(const std::vector<std::shared_ptr<Boundary>>& boundaries, double scale=1.0);

        bool is_inside(const vec3& r) const;  // finding the correct tetrahedra would be easier with a mapping rectangular cell -> tetrahedras that overlap with this cell
        std::shared_ptr<Boundary::CollisionEvent> collision_point(vec3 orig, vec3 direction, std::shared_ptr<Boundary::CollisionEvent> last) const;

        std::vector<vec3> random_positions(unsigned n, std::function<double(vec3)> weights = constant_weight) const;

        double volume_integral(std::function<double(vec3)> integrand) const {
            double res = 0;
            for(auto tetra : tetrahedra)
                res+=tetra.volume * integrand(tetra.center);
            return res;
        }

        // vec3 numeric_gradient(std::function<double(vec3)> fn) const;  // always better to use an analytical gradient since this method is expensive
                                                                         // calculation would be done be determining the comprising tetrahedra first and afterwards 
                                                                         // summing the delta/edge_length * edge_direction over all tetrahedra edges divided by the
                                                                         // amount of edges (6)

        double volume() const {
            double sum = 0;
            for(const auto& tetra : tetrahedra)
                sum+=tetra.volume;
            return sum;
        }
        
        double area() const {
            double sum = 0;
            for(const auto& bound : boundaries)
                sum+=bound->area();
            return sum;
        }

        vec3 center() const {
            return (box.first + box.second) / 2;
        }
    };
}
