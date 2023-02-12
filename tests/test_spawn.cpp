#include <fstream>
#include <utility>
#include <memory>
#include <vector>
#include <iostream>
#include <iostream>

#include "boltzi_types.h"
#include "phonon.h"
#include "boundary.h"
#include "boundaries/specular_vacuum_boundary.h"

using namespace std;

int main() {
    using namespace boltzi;
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

    auto res = sg.random_positions(50000);

    ofstream output("test_spawn.out.csv");

    for(auto pos : res) {
        output << pos[0] << "," << pos[1] << "," << pos[2] << endl;
    }
    output.close();

    return 0;
}
