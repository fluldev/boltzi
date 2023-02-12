#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>

#include <armadillo>

#include "scatter_processes/simple_parametric.h"
#include "branches/linear_dispersion_branch.h"
#include "samples/simple_sample.h"
#include "boundary.h"
#include "boundaries/simple_isothermal_boundary.h"
#include "boltzi_util.h"

using namespace boltzi;

void create_primitives(std::vector<Boundary::Primitive>& out, vec3 p1, vec3 p2, vec3 p3, vec3 p4) {
    // points in clockwise order when looking from the outside
    vec3 normal = -arma::normalise(arma::cross(p2-p1, p4-p1));
    vec3 e1 = arma::normalise(p2-p1);
    vec3 e2 = arma::cross(normal, e1);
    out.push_back(Boundary::Primitive{util::Triangle{p1, p2, p4}, normal, e1, e2});
    out.push_back(Boundary::Primitive{util::Triangle{p2, p3, p4}, normal, e1, e2});
}

int main() {
    double omega_debye = 300e11;
    double n_eff = 1e31;
    double v_g = 1;
    double init_temp = 40;
    double deviational_temp = 10;
    std::vector<Boundary::Primitive> primitives;
    vec3 c[] = {
        vec3{0,0,0},
        vec3{0,1,0},
        vec3{1,1,0},
        vec3{1,0,0},

        vec3{0,0,1},
        vec3{0,1,1},
        vec3{1,1,1},
        vec3{1,0,1}
    };
    create_primitives(primitives, c[0], c[1], c[2], c[3]);
    create_primitives(primitives, c[7], c[6], c[5], c[4]);

    create_primitives(primitives, c[3], c[7], c[4], c[0]);
    create_primitives(primitives, c[1], c[5], c[6], c[2]);

    create_primitives(primitives, c[0], c[4], c[5], c[1]);
    create_primitives(primitives, c[2], c[6], c[7], c[3]);

    std::cout<<"Created Primitive."<<std::endl;

    SampleGeometry geometry(std::vector<std::shared_ptr<Boundary>>{
        std::shared_ptr<Boundary>(new boundaries::SimpleIsothermalBoundary(primitives, 1)),
    });


    std::cout<<"Created SampleGeometry and Boundary"<<std::endl;
    
    auto scatter_prc = {
        std::shared_ptr<ThreePhononProcess>(new scatter::SimpleParametricThreePhononProcess<1,1>(1, v_g, omega_debye)),
        std::shared_ptr<ThreePhononProcess>(new scatter::SimpleParametricThreePhononProcess<2,2>(1, v_g, omega_debye))
    };
    std::vector<std::shared_ptr<Branch>> branches{
        std::shared_ptr<Branch>(new lindisp::LinearDispBranch(
            std::vector<std::shared_ptr<TwoPhononProcess>>{},
            std::vector<std::shared_ptr<ThreePhononProcess>>(scatter_prc),
            v_g, omega_debye
        ))  
    };

    Material material(branches);

    std::cout<<"Created Material."<<std::endl;

    samples::SimpleSample sample(
            n_eff,  // n_eff
            std::make_tuple(1,1,1),  // n_bins
            geometry,
            material,
            [deviational_temp] (vec3)->double {return deviational_temp;}
    );


    std::cout<<"Created SimpleSample."<<std::endl;
    
    std::cout<<"Spawning volume phonons."<<std::endl;
    sample.spawn_phonons([init_temp] (vec3) -> double {return init_temp;});  // constant temperature of one Kelvin
    std::cout<<"Volume phonons spawned."<<std::endl;

    std::ofstream file("test_lin_material.out.csv");
    for(auto pos : sample.phonons.phonons) {
        for(short i = 0; i<3; ++i)
            file<<pos.r[i]<<",";
        file<<arma::norm(pos.k)<<"\n";
    }
    file.close();

    std::cout<<"Calculate Temperature"<<std::endl;
    sample.calc_bin_temps();
    std::cout<<"Average temperature should be "<<init_temp<<": "<<sample.avg_temp()<<std::endl;

    std::cout<<"Average pseudo temperatures should be around "<<init_temp<<": "<<std::endl;
    auto pts = sample.avg_pseudo_temps();
    for(auto& [material_idx, material] : pts)
        for(auto& [branch_idx, branch] : material) {
            std::cout<<"Material "<<material_idx<<", Branch "<<branch_idx<<": "<<std::endl;
            for(auto pt : branch)
                std::cout<<pt<<std::endl;
        }

    sample.phonons.clear();

    std::cout<<"Spawning surface phonons."<<std::endl;
    sample.spawn_surface_phonons(1);
    std::cout<<"Surface phonons spawned."<<std::endl;

    std::ofstream file2("test_lin_material_surf.out.csv");
    for(auto pos : sample.phonons.phonons) {
        for(short i = 0; i<3; ++i)
            file2<<pos.r[i]<<",";
        for(short i = 0; i<3; ++i)
            file2<<pos.k[i]<<",";
        file2<<"\n";
    }
    file2.close();

    return 0;
}
