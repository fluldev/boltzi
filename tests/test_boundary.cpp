#include <utility>
#include <memory>
#include <vector>
#include <iostream>

#include "boltzi_types.h"
#include "phonon.h"
#include "boundary.h"
#include "boundaries/specular_vacuum_boundary.h"

using namespace std;

int main() {
    using namespace boltzi;
    boundaries::SpecularVacuumBoundary bound({
        Boundary::Primitive{{vec3{0,0,0}, vec3{0,1,0}, vec3{0,0,1}}, vec3{1,0,0}, vec3{0,1,0}, vec3{0,0,1}}
    });

    PhononData ph(vec3{-1,1,0}, vec3{0,0,0}, 0, 0, 1);
    bound.reflect(ph, std::shared_ptr<Boundary::CollisionEvent>(new Boundary::CollisionEvent{vec3{0,.5,.25}, std::shared_ptr<Boundary>(&bound), 0})).print();
    bound.lowest_corner().print();
    bound.highest_corner().print();

    // Format: Corner1, Corner2, Corner3, Normal, e_1, e_2
    // e_1, e_2 don't matter here
    SampleGeometry sg(std::vector<std::shared_ptr<Boundary>>{
        std::shared_ptr<Boundary>(new boundaries::SpecularVacuumBoundary(std::vector<Boundary::Primitive>{
            Boundary::Primitive{{vec3{0,0,0}, vec3{0,1,0}, vec3{0,0,1}}, vec3{1,0,0}, vec3{0,1,0}, vec3{0,0,1}},
            Boundary::Primitive{{vec3{0,1,1}, vec3{0,1,0}, vec3{0,0,1}}, vec3{1,0,0}, vec3{0,1,0}, vec3{0,0,1}},
        })),
        std::shared_ptr<Boundary>(new boundaries::SpecularVacuumBoundary(std::vector<Boundary::Primitive>{
            Boundary::Primitive{{vec3{1,0,0}, vec3{1,1,0}, vec3{1,0,1}}, vec3{-1,0,0}, vec3{0,1,0}, vec3{0,0,1}},
            Boundary::Primitive{{vec3{1,1,1}, vec3{1,1,0}, vec3{1,0,1}}, vec3{-1,0,0}, vec3{0,1,0}, vec3{0,0,1}},
        })),

        std::shared_ptr<Boundary>(new boundaries::SpecularVacuumBoundary(std::vector<Boundary::Primitive>{
            Boundary::Primitive{{vec3{0,0,0}, vec3{0,1,0}, vec3{1,0,0}}, vec3{0,0,1}, vec3{0,1,0}, vec3{0,0,1}},
            Boundary::Primitive{{vec3{1,1,0}, vec3{0,1,0}, vec3{1,0,0}}, vec3{0,0,1}, vec3{0,1,0}, vec3{0,0,1}},
        })),
        std::shared_ptr<Boundary>(new boundaries::SpecularVacuumBoundary(std::vector<Boundary::Primitive>{
            Boundary::Primitive{{vec3{0,0,1}, vec3{0,1,1}, vec3{1,0,1}}, vec3{1,0,-1}, vec3{0,1,0}, vec3{0,0,1}},
            Boundary::Primitive{{vec3{1,1,1}, vec3{0,1,1}, vec3{1,0,1}}, vec3{1,0,-1}, vec3{0,1,0}, vec3{0,0,1}},
        })),

        std::shared_ptr<Boundary>(new boundaries::SpecularVacuumBoundary(std::vector<Boundary::Primitive>{
            Boundary::Primitive{{vec3{0,0,0}, vec3{1,0,0}, vec3{0,0,1}}, vec3{0,1,0}, vec3{0,1,0}, vec3{0,0,1}},
            Boundary::Primitive{{vec3{1,0,1}, vec3{1,0,0}, vec3{0,0,1}}, vec3{0,1,0}, vec3{0,1,0}, vec3{0,0,1}},
        })),
        std::shared_ptr<Boundary>(new boundaries::SpecularVacuumBoundary(std::vector<Boundary::Primitive>{
            Boundary::Primitive{{vec3{0,1,0}, vec3{1,1,0}, vec3{0,1,1}}, vec3{0,-1,0}, vec3{0,1,0}, vec3{0,0,1}},
            Boundary::Primitive{{vec3{1,1,1}, vec3{1,1,0}, vec3{0,1,1}}, vec3{0,-1,0}, vec3{0,1,0}, vec3{0,0,1}},
        }))
    });

    cout<<sg.is_inside(vec3{.5,.5,.5})<<endl;
    cout<<sg.is_inside(vec3{1.5,.5,.5})<<endl;

    //std::optional<std::pair<const Boundary&, std::pair<vec3, unsigned>>> collision_check(vec3 orig, vec3 direction, double path_length) const;
    auto coll = sg.collision_point(vec3{.5,.5,.5}, vec3{1,1,0});
    if(coll!=nullptr) 
        coll->position.print();
    else cout<<"No collision."<<endl;
    return 0;
}
