#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>

#include <armadillo>

#include "simulation_controller.h"
#include "scatter_processes/simple_parametric.h"
#include "branches/linear_dispersion_branch.h"
#include "samples/simple_sample.h"
#include "boundaries/specular_vacuum_boundary.h"
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

    std::vector<std::shared_ptr<Boundary>> bounds{
        std::shared_ptr<Boundary>(new boundaries::SimpleIsothermalBoundary(primitives, 1)),
    };
    SampleGeometry geometry(bounds);


    std::cout<<"Created SampleGeometry and Boundary"<<std::endl;
    
    auto scatter_prc = {
        std::shared_ptr<ThreePhononProcess>(new scatter::SimpleParametricThreePhononProcess<1,1>(1e-6, v_g, omega_debye)),
        std::shared_ptr<ThreePhononProcess>(new scatter::SimpleParametricThreePhononProcess<2,2>(1e-6, v_g, omega_debye))
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

    NonLinearDSMC simulation(sample, 0.1);
    std::cout<<"Created SimulationController"<<std::endl;

    // add test phonon
    sample.spawn_phonons([] (vec3)->double {return 1;});
    sample.overwrite_temps([] (vec3) -> double {return 1;});

    std::ofstream file("test_advection.out.csv");

    auto phonon = sample.phonons.phonons.begin();
    for(unsigned i = 0; i<2000; ++i) {
        simulation.advection_step(phonon, .1);
        simulation.collision_step(phonon, .1);
        for(short i = 0; i<3; ++i) file<<phonon->r[i]<<",";
        for(short i = 0; i<3; ++i) file<<phonon->k[i]<<",";
        file<<"\n";
    }
    file.close();
}
