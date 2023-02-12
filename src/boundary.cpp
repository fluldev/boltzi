#include <algorithm>
#include <iterator>
#include <tuple>
#include <vector>
#include <memory>
#include <optional>
#include <armadillo>
#include <tetgen.h>
#include <string>

#include <iostream>


#include "headers/boltzi_config.h"
#include "headers/boltzi_types.h"
#include "headers/boundary.h"
#include "headers/util/random.h"

namespace boltzi {
    std::shared_ptr<Boundary::CollisionEvent> Boundary::collision_point(vec3 orig, vec3 dir_vec, std::shared_ptr<Boundary::CollisionEvent> last) const {
        std::vector<std::tuple<double, vec3, unsigned, bool>> candidates;
        for(auto it = primitives.cbegin(); it!=primitives.cend();++it) {
            // TODO: should be within primitives for the possibility of other primitives
            double eq_left = arma::dot(dir_vec, it->normal);
            if(eq_left>=0 && is_active()) continue; // moves away from the surface 
            vec3 pos_diff = it->tri.points[0] - orig;
            double eq_right = arma::dot(pos_diff, it->normal);
            double lambd = eq_right / eq_left;
            if(lambd<0) continue;  // collision_point is in the wrong direction (in time)
            vec3 surf_coll_pt = orig + dir_vec*lambd;
            if(it->tri.is_inside(surf_coll_pt)) {
                // TODO: should depend on scale (makes sure that a collison with a measurement boundary doesn't happen twice)
                if(!is_active() && last && arma::norm(last->position - surf_coll_pt) < config::DOUBLE_COLL_IGNORE_RADIUS)
                    continue;
                candidates.push_back({arma::norm(surf_coll_pt - orig), surf_coll_pt, it-primitives.cbegin(), eq_left<0});
            }
        }
        if(candidates.empty()) return nullptr;
        auto min_ele = std::min_element(candidates.begin(), candidates.end(), 
            [] (const std::tuple<double, vec3, unsigned, bool> &a, const std::tuple<double, vec3, unsigned, bool>  &b) -> bool {return std::get<0>(a) < std::get<0>(b);}
        );
        return std::shared_ptr<Boundary::CollisionEvent>{new Boundary::CollisionEvent{
            999999,
            std::get<1>(*min_ele), 
            std::shared_ptr<const Boundary>(this->shared_from_this()), 
            std::get<2>(*min_ele),
            std::get<3>(*min_ele)
        }};
    }
    
    SampleGeometry::SampleGeometry(const std::vector<std::shared_ptr<Boundary>>& boundaries, double scale) : boundaries(boundaries), scale(scale) {
        struct Face{
            unsigned p1,p2,p3; 
        };
        
        if(boundaries.empty()) throw BoltziException("SampleGeometry must have at least one boundary.");

        //accumulate vertices and faces
        std::vector<vec3> unique_vertices; 
        std::vector<Face> faces;
        for(auto boundary : boundaries) {
            // passive boundaries do not physically exist
            if(!boundary->is_active())
                continue;
            for(auto prim : boundary->primitives) {
                unsigned i=0;
                unsigned face_assoc[3];
                for(vec3 vertex : prim.tri.points) {
                    auto index = std::find_if(unique_vertices.begin(), unique_vertices.end(),
                        [this,&vertex] (const vec3& ele) -> bool {return arma::norm(ele-vertex) < config::WELD_RADIUS * this->scale;}
                    );
                    if(index == unique_vertices.end()) {
                        unique_vertices.push_back(vertex);
                        face_assoc[i] = unique_vertices.size()-1;
                    }
                    else {
                        face_assoc[i] = index-unique_vertices.begin();
                    }
                    ++i;
                }
                faces.push_back(Face{face_assoc[0], face_assoc[1], face_assoc[2]});
            }
        }

        // find box that fits all boundaries
        for(auto it = unique_vertices.begin(); it!=unique_vertices.end(); ++it) {
            box.first[0] = std::min(box.first[0], (*it)[0]);
            box.second[0] = std::max(box.second[0], (*it)[0]);
            box.first[1] = std::min(box.first[1], (*it)[1]);
            box.second[1] = std::max(box.second[1], (*it)[1]);
            box.first[2] = std::min(box.first[2], (*it)[2]);
            box.second[2] = std::max(box.second[2], (*it)[2]);
        }
        
        // translate into tetgen format
        // TODO: unsure if those have destructors -> memory leak?
        tetgenio in, out;
        in.numberofpoints = unique_vertices.size();
        in.pointlist = new REAL [in.numberofpoints*3];
        for(unsigned i=0; i<unique_vertices.size(); ++i) {
            in.pointlist[i*3] = unique_vertices[i][0]; 
            in.pointlist[i*3+1] = unique_vertices[i][1]; 
            in.pointlist[i*3+2] = unique_vertices[i][2]; 
        }

        in.numberoffacets = faces.size();
        in.facetlist = new tetgenio::facet[in.numberoffacets];
        for(unsigned i=0; i<faces.size(); ++i) {
            in.facetlist[i].numberofpolygons = 1;
            in.facetlist[i].numberofholes = 0;
            in.facetlist[i].holelist = NULL;
            in.facetlist[i].polygonlist = new tetgenio::polygon[1];
            auto polygon = in.facetlist[i].polygonlist;
            polygon->numberofvertices = 3;
            polygon->vertexlist = new int[3];
            polygon->vertexlist[0] = faces[i].p1;
            polygon->vertexlist[1] = faces[i].p2;
            polygon->vertexlist[2] = faces[i].p3;
        }


        // TODO: facetmarks should be used for isothermal boundaries
        in.facetmarkerlist = new int[in.numberoffacets];
        for(int i = 0; i<in.numberoffacets; ++i)
            in.facetmarkerlist[i] = 0;

        // TODO: optimize settings for tetrahedralize
        char opt[] = "pq1.414a0.1Q";
        tetrahedralize(opt, &in, &out);

        //tetrahedra.resize(out.numberoftetrahedra);
        for(int i=0; i<out.numberoftetrahedra; ++i)
            tetrahedra.push_back(util::Tetrahedron{
                // vertex index starts with 1
                vec3{out.pointlist[out.tetrahedronlist[i*4+0]*3], out.pointlist[out.tetrahedronlist[i*4+0]*3+1], out.pointlist[out.tetrahedronlist[i*4+0]*3+2]},
                vec3{out.pointlist[out.tetrahedronlist[i*4+1]*3], out.pointlist[out.tetrahedronlist[i*4+1]*3+1], out.pointlist[out.tetrahedronlist[i*4+1]*3+2]},
                vec3{out.pointlist[out.tetrahedronlist[i*4+2]*3], out.pointlist[out.tetrahedronlist[i*4+2]*3+1], out.pointlist[out.tetrahedronlist[i*4+2]*3+2]},
                vec3{out.pointlist[out.tetrahedronlist[i*4+3]*3], out.pointlist[out.tetrahedronlist[i*4+3]*3+1], out.pointlist[out.tetrahedronlist[i*4+3]*3+2]},
            });
    }

    bool SampleGeometry::is_inside(const vec3& r) const {
        for(const auto& tetra : tetrahedra)
            if(tetra.is_inside(r))
                return true;
        return false;
    }
    std::shared_ptr<Boundary::CollisionEvent> SampleGeometry::collision_point(vec3 orig, vec3 direction, std::shared_ptr<Boundary::CollisionEvent> last) const {
        std::vector<std::shared_ptr<Boundary::CollisionEvent>> candidates;
        for(auto it = boundaries.begin(); it!=boundaries.end(); ++it) {
            auto coll = (*it)->collision_point(orig, direction, last);
            if(coll == nullptr) continue;  // no collision occured
            coll->boundary_idx = it-boundaries.begin();
            candidates.push_back(coll);
        }

        if(candidates.empty()) return nullptr;

        auto it = std::min_element(candidates.begin(), candidates.end(), 
            [orig] (const auto& a, const auto& b) -> bool {return arma::norm(a->position-orig) < arma::norm(b->position-orig);}
        );

        return *it;
    }

    std::vector<vec3> SampleGeometry::random_positions(unsigned n, std::function<double(vec3)> weights) const {
        std::vector<double> tetra_weights;
        std::vector<vec3> res;
        res.reserve(n);
        tetra_weights.reserve(tetrahedra.size());
        for(auto& tetra : tetrahedra)
            tetra_weights.push_back(tetra.volume * weights(tetra.center));

        for(unsigned i = 0; i<n; ++i) {
            res.push_back(
                util::random.random_pos_tetrahedron(
                    tetrahedra[util::random.draw_index_weighted(tetra_weights)]
            ));
        }
        return res;
    }
    

    //TODO: inefficient way of handling this (returned vector is way to large)
    std::vector<std::tuple<vec3,vec3,vec3,vec3>> Boundary::random_positions(unsigned n, std::function<double(vec3)> weights) const {
        std::vector<double> tri_weights;
        std::vector<std::tuple<vec3,vec3,vec3,vec3>> res;
        res.reserve(n);
        tri_weights.reserve(primitives.size());
        for(auto prim : primitives)
            tri_weights.push_back(prim.tri.area * weights(prim.tri.center));

        for(unsigned i = 0; i<n; ++i) {
            const auto& prim = primitives[util::random.draw_index_weighted(tri_weights)];
            res.push_back(std::tuple<vec3,vec3,vec3,vec3>{
                util::random.random_pos_triangle(
                    prim.tri
                ),
                prim.normal,
                prim.e_1,
                prim.e_2
            });
        }
        return res;
    }
}
