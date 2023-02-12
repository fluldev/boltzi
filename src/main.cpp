#include <exception>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <vector>
#include <memory>
#include <string>
#include <sstream>

#include <yaml-cpp/node/node.h>
#include <yaml-cpp/yaml.h>
#include <exprtk.hpp>

// boltzi framework
#include "headers/boltzi.h"

// Specialisations
// Boundaries
#include "headers/boltzi_exception.h"
#include "headers/boundaries/vacuum_boundary.h"
#include "headers/boundaries/simple_isothermal_boundary.h"

// Branches
#include "headers/boundary.h"
#include "headers/branches/linear_dispersion_branch.h"

// Samples
#include "headers/material.h"
#include "headers/sample.h"
#include "headers/sample_geometry_factory.h"
#include "headers/samples/simple_sample.h"

// Scatter processes
#include "headers/scatter_processes/simple_parametric.h"
#include "headers/simulation_controller.h"

// Artifacts
#include "headers/artifacts.h"
#include "headers/util/constants.h"
#include "headers/util/pool.h"

// Silent mode
#include "headers/silencer.h"

using namespace boltzi;
using config::Silencer;

void help_msg() {
    std::cout<<"Usage: boltzi CONFIG [FLAG|OPTION...]\n";
    std::cout<<"Flags:\n";
    std::cout<<"--help, -h:\tShow this help.\n";
    std::cout<<"--silent:\tDo not produce any output to stdout (stderr is still active).\n";
    std::cout<<"--less_output:\tProduce less output to stdout (stderr is still active).\n";
    std::cout<<"Options:\n";
    std::cout<<"-a, --artifact NAME TIME [FILENAME]:\tMeasure NAME at time TIME and store the results in file with "
                                                      "optional name FILENAME (default name is \"NAME_TIME.boltzi_artifact\").\n";
    std::cout<<"-r, --resolution NX NY NZ:\tSet the measurement resolution to NX,NY,NZ along the x-,y- and z-axis for"
                                            "all following artifacts until another such option is encountered.\n";
    std::cout<<"-s, --save TIME [FILENAME]:\tSave simulation state at time TIME to file with optional name FILENAME (default name is \"state_TIME.boltzi_state\").\n";
    std::cout<<"-l, --load FILE:\tLoad simulation state from file FILE.\n";
    std::cout<<"-t, --time TIME:\tStarting time of the simulation (default is 0). Should match the save-time of the loaded simulation state if -l/--load is used.\n";
    std::cout<<"-C, --directory PATH:\tSets the directory in which the following artifacts get stored. Config artifacts get saved in the directory of the last of those options.\n";
    std::cout<<"-b:\tFilename of boundary current logs (default is boundarycurrent).";
    exit(0);
}
void eval_scatter_processes(
        std::vector<std::shared_ptr<TwoPhononProcess>>& two_phonon_processes, 
        std::vector<std::shared_ptr<ThreePhononProcess>>& three_phonon_processes, 
        const YAML::Node& this_branch
) {
    for(const auto& prc_arr : this_branch["two_phonon_processes"]) {
        const std::pair<YAML::Node, YAML::Node>& prc = *(prc_arr.begin());
        if(prc.first.as<std::string>() == std::string("parametric")) {
            double fit_constant = prc.second["fit_constant"].as<double>();
            unsigned omega_power = prc.second["wavevector_power"].as<unsigned>();
            unsigned temp_power = prc.second["temperature_power"].as<unsigned>();

            double debye_frequency;
            if(this_branch["debye_frequency"])
                debye_frequency = this_branch["debye_frequency"].as<double>();
            else 
                debye_frequency = this_branch["debye_temperature"].as<double>() * util::BMC/util::HBAR;

            two_phonon_processes.push_back(std::shared_ptr<TwoPhononProcess>(new scatter::SimpleParametricTwoPhononProcess{
                temp_power,
                omega_power,
                fit_constant,
                this_branch["velocity"].as<double>(),
                debye_frequency 
            }));
        }
        else if(prc.first.as<std::string>() == std::string("parametric_exp")) {
            double fit_constant = prc.second["fit_constant"].as<double>();
            unsigned omega_power = prc.second["wavevector_power"].as<unsigned>();
            unsigned temp_power = prc.second["temperature_power"].as<unsigned>();
            double expfac = util::HBAR * this_branch["debye_frequency"].as<double>() / util::BMC;
            if(prc.second["exp_constant"])
                expfac = prc.second["exp_constant"].as<double>();

            double debye_frequency;
            if(this_branch["debye_frequency"])
                debye_frequency = this_branch["debye_frequency"].as<double>();
            else 
                debye_frequency = this_branch["debye_temperature"].as<double>() * util::BMC/util::HBAR;

            two_phonon_processes.push_back(std::shared_ptr<TwoPhononProcess>(new scatter::SimpleParametricExpTwoPhononProcess{
                temp_power,
                omega_power,
                fit_constant,
                this_branch["velocity"].as<double>(),
                debye_frequency,
                expfac
            }));
        }
        else {
            std::cerr<<"Unkown two phonon process type \""<<prc.first.as<std::string>()<<"\". Quitting."<<std::endl;
            exit(-1);
        }
    }
    for(const auto& prc_arr : this_branch["three_phonon_processes"]) {
        const std::pair<YAML::Node, YAML::Node>& prc = *(prc_arr.begin());
        if(prc.first.as<std::string>() == std::string("parametric")) {
            double fit_constant = prc.second["fit_constant"].as<double>();
            unsigned omega_power = prc.second["wavevector_power"].as<unsigned>();
            unsigned temp_power = prc.second["temperature_power"].as<unsigned>();

            double debye_frequency;
            if(this_branch["debye_frequency"])
                debye_frequency = this_branch["debye_frequency"].as<double>();
            else 
                debye_frequency = this_branch["debye_temperature"].as<double>() * util::BMC/util::HBAR;

            bool is_umklapp=false;
            if(prc.second["is_umklapp"])
                is_umklapp = prc.second["is_umklapp"].as<bool>();

            three_phonon_processes.push_back(std::shared_ptr<ThreePhononProcess>(new scatter::SimpleParametricThreePhononProcess{
                temp_power,
                omega_power,
                fit_constant,
                this_branch["velocity"].as<double>(),
                debye_frequency,
                is_umklapp
            }));
        }
        else if(prc.first.as<std::string>() == std::string("parametric_exp")) {
            double fit_constant = prc.second["fit_constant"].as<double>();
            unsigned omega_power = prc.second["wavevector_power"].as<unsigned>();
            unsigned temp_power = prc.second["temperature_power"].as<unsigned>();
            double expfac = util::HBAR * this_branch["debye_frequency"].as<double>() / util::BMC;
            if(prc.second["exp_constant"])
                expfac = prc.second["exp_constant"].as<double>();
            
            double debye_frequency;
            if(this_branch["debye_frequency"])
                debye_frequency = this_branch["debye_frequency"].as<double>();
            else 
                debye_frequency = this_branch["debye_temperature"].as<double>() * util::BMC/util::HBAR;
            
            bool is_umklapp=false;
            if(prc.second["is_umklapp"])
                is_umklapp = prc.second["is_umklapp"].as<bool>();

            three_phonon_processes.push_back(std::shared_ptr<ThreePhononProcess>(new scatter::SimpleParametricExpThreePhononProcess{
                temp_power,
                omega_power,
                fit_constant,
                this_branch["velocity"].as<double>(),
                debye_frequency,
                is_umklapp,
                expfac
            }));
        }
        else {
            std::cerr<<"Unkown three phonon process type \""<<prc.first.as<std::string>()<<"\". Quitting."<<std::endl;
            exit(-1);
        }
    } 

}

// THIS CAN'T BE USED IN MULTIPLE THREADS!
// Use MultithreadExprWrapper instead in a multithreaded environment.
struct ExprWrapper {
    double x,y,z;
    double scale;
    std::string expr_str{};
    exprtk::expression<double> expr;
    exprtk::symbol_table<double> sym_tbl;
    ExprWrapper() = default;
    ExprWrapper(std::string expr_str, double scale) : expr_str{expr_str}, scale(scale) { 
        exprtk::parser<double> expr_parser;
        sym_tbl.add_variable("x", x);
        sym_tbl.add_variable("y", y);
        sym_tbl.add_variable("z", z);
        expr.register_symbol_table(sym_tbl);
        if(!expr_parser.compile(expr_str, expr)) {
            throw BoltziException((std::stringstream("Can't parse expression \"")<<expr_str<<"\".").str().c_str());
        }
    }
    ExprWrapper(const ExprWrapper& other) : ExprWrapper(other.expr_str, other.scale) {}
    ExprWrapper& operator=(ExprWrapper other) {
        expr_str = other.expr_str;
        scale = other.scale;
        exprtk::parser<double> expr_parser;
        sym_tbl.add_variable("x", x);
        sym_tbl.add_variable("y", y);
        sym_tbl.add_variable("z", z);
        expr.register_symbol_table(sym_tbl);
        if(!expr_parser.compile(expr_str, expr)) {
            throw BoltziException((std::stringstream("Can't parse expression \"")<<expr_str<<"\".").str().c_str());
        }
        return *this;
    }
    double operator()(const vec3& r) {
        x = r[0]/scale;
        y = r[1]/scale;
        z = r[2]/scale;
        return expr.value();
    }
    vec3 deriv(const vec3& r) {
        x = r[0]/scale;
        y = r[1]/scale;
        z = r[2]/scale;
        return vec3{
            exprtk::derivative(expr, "x"),
            exprtk::derivative(expr, "y"),
            exprtk::derivative(expr, "z")
        } / scale;
    }
};

struct MultithreadExprWrapper {
    std::unordered_map<std::thread::id, ExprWrapper> tbl;
    ExprWrapper to_copy;
    MultithreadExprWrapper(const ExprWrapper& obj) : to_copy(obj) {}

    MultithreadExprWrapper() {}
    MultithreadExprWrapper& operator=(const ExprWrapper& other) {
        to_copy = other;
        tbl.clear();
        return *this;
    }

    inline void create_on_first() {
        if(!tbl.contains(std::this_thread::get_id()))
            tbl[std::this_thread::get_id()] = to_copy;
    }

    double operator()(const vec3& r) {
        create_on_first(); 
        return tbl[std::this_thread::get_id()](r);
    }

    vec3 deriv(const vec3& r) {
        create_on_first();
        return tbl[std::this_thread::get_id()].deriv(r);
    }
};

int main(int argc, char* argv[]) {
    if(argc < 2) {
        std::cerr<<"Missing config file. Quitting."<<std::endl;
        return -1;
    }

    std::vector<std::shared_ptr<SimulationController::Artifact>> artfcs;
    
    if(strcmp(argv[1], "-h")==0 || strcmp(argv[1],  "--help")==0) help_msg();

    const char* config_file_str = argv[1];
    YAML::Node config;
    try {
         config = YAML::LoadFile(config_file_str);
    }
    catch(...) {
        std::cerr<<"Can not find or parse config file \""<<config_file_str<<"\". Quitting."<<std::endl;
        return -1;
    }

    std::tuple<unsigned, unsigned, unsigned> current_resolution; 
    // The following must be thread local in order to function
    // in parallel.
    ExprWrapper init_temperature;
    MultithreadExprWrapper deviational_temperature_expr;
    double time;
    double time_step;
    std::shared_ptr<SampleGeometry> geometry;
    std::vector<std::shared_ptr<Material>> materials;
    std::shared_ptr<SimulationController> simulation;
    std::shared_ptr<Sample> sample;

    std::string out_path;
    std::string boundary_current_log_fname = "boundarycurrent";
    
    try {
        current_resolution = std::make_tuple(
            config["resolution"]["x"].as<unsigned>(),
            config["resolution"]["y"].as<unsigned>(),
            config["resolution"]["z"].as<unsigned>()
        );


        time = config["time"].as<double>();
        if(time<0) {
            std::cerr<<"Config entry \"time\" can't be negative. Quitting."<<std::endl;
            return -1;
        }
        
        time_step = config["time_step"].as<double>();
        if(time_step<0) {
            std::cerr<<"Config entry \"time\" can't be negative. Quitting."<<std::endl;
            return -1;
        }
        if(time_step>time) {
            std::cerr<<"Config entry \"time\" can't be bigger than config entry \"time\". Quitting."<<std::endl;
            return -1;
        }
        for(const auto& mat_arr : config["materials"]) {
            // TODO: material name should be used in samples with multiple materials.
            const std::pair<YAML::Node,YAML::Node>& mat = *(mat_arr.begin());
            std::vector<std::shared_ptr<Branch>> branches;

            for(const auto& branch_arr : mat.second["branches"]) {
                const std::pair<YAML::Node,YAML::Node> branch = *(branch_arr.begin());
                if(branch.first.as<std::string>() == std::string("linear_dispersion")) {
                    double velocity = branch.second["velocity"].as<double>();
                    double omega_debye = branch.second["debye_frequency"].as<double>();
                    std::vector<std::shared_ptr<TwoPhononProcess>> two_phonon_processes;
                    std::vector<std::shared_ptr<ThreePhononProcess>> three_phonon_processes;

                    eval_scatter_processes(two_phonon_processes, three_phonon_processes, branch.second);
                    branches.push_back(std::shared_ptr<Branch>(new lindisp::LinearDispBranch{
                        std::move(two_phonon_processes),
                        std::move(three_phonon_processes),
                        velocity,
                        omega_debye
                    }));
                }
                else {
                    std::cerr<<"Unkown branch type \""<<branch.first.as<std::string>()<<"\". Quitting."<<std::endl;
                    return -1;
                }
            }
            materials.push_back(std::shared_ptr<Material>(new Material{
                branches
            }));
        }

        double scale = 1;
        if(config["sample"]["scale"])
            scale = config["sample"]["scale"].as<double>();
        geometry = std::shared_ptr<SampleGeometry>(new SampleGeometry(load_geometry(config["sample"]["file"].as<std::string>(), scale)));
        
        bool zero_gradient = false;
        if(config["zero_gradient"])
            zero_gradient=config["zero_gradient"].as<bool>();

        deviational_temperature_expr = ExprWrapper(config["deviational_temperature"].as<std::string>(), scale);
        init_temperature = ExprWrapper(config["init_temperature"].as<std::string>(), scale);


        // TODO: Type of sample could be automatically detected as simple if there is only one material.
        if(config["sample"]["type"].as<std::string>() == std::string("simple")) {
            sample = std::shared_ptr<Sample>(new samples::SimpleSample{
                config["effective_phonons"].as<double>(),
                current_resolution,
                *geometry,
                *(materials.front()),
                [&deviational_temperature_expr](vec3 r) ->double {return deviational_temperature_expr(r);},
                [&deviational_temperature_expr](vec3 r) ->vec3 {return deviational_temperature_expr.deriv(r);}
            });
        }
        else {
            std::cerr<<"Unkown sample type \""<<config["sample"]["type"].as<std::string>()<<"\". Quitting."<<std::endl;
            return -1;
        }

        if(config["simulation"].as<std::string>() == std::string("nonlinear"))
            simulation = std::shared_ptr<SimulationController>(new NonLinearDSMC{
                *sample,
                time_step,
                config["boundary_current_binwidth"].as<double>(),
                zero_gradient
            });
        else if(config["simulation"].as<std::string>() == std::string("linear")) {
            simulation = std::shared_ptr<SimulationController>(new LinearDSMC{
                *sample,
                time_step,
                config["boundary_current_binwidth"].as<double>(),
            });
        }
        else {
            std::cerr<<"Unkown simulation type \""<<config["simulation"].as<std::string>()<<"\". Quitting."<<std::endl;
            return -1;
        }

        // TODO: Check if values are valid.
        for(const auto& artifact_arr : config["artifacts"]) {
            const std::pair<YAML::Node,YAML::Node>& artifact = *(artifact_arr.begin());
            if(artifact.first.as<std::string>() == std::string("temperature")) {
                std::string name = "temperature";
                if(artifact.second["name"])
                    name = artifact.second["name"].as<std::string>();
                artfcs.push_back(std::shared_ptr<SimulationController::Artifact>(new artifacts::Temperature{
                    name,
                    artifact.second["time"].as<double>(),
                    std::make_tuple(
                        artifact.second["resolution"]["x"].as<unsigned>(), 
                        artifact.second["resolution"]["y"].as<unsigned>(), 
                        artifact.second["resolution"]["z"].as<unsigned>()
                    ),    
                }));
            }
            else if(artifact.first.as<std::string>() == std::string("heatcurrent")) {
                std::string name = "heatcurrent";
                if(artifact.second["name"])
                    name = artifact.second["name"].as<std::string>();
                artfcs.push_back(std::shared_ptr<SimulationController::Artifact>(new artifacts::HeatCurrent{
                    out_path+name,
                    artifact.second["time"].as<double>(),
                    std::make_tuple(
                        artifact.second["resolution"]["x"].as<unsigned>(), 
                        artifact.second["resolution"]["y"].as<unsigned>(), 
                        artifact.second["resolution"]["z"].as<unsigned>()
                    )    
                }));
            }
            else if(artifact.first.as<std::string>() == std::string("meanfreepath")) {
                std::string name = "meanfreepath";
                if(artifact.second["name"])
                    out_path+name = artifact.second["name"].as<std::string>();
                artfcs.push_back(std::shared_ptr<SimulationController::Artifact>(new artifacts::MeanFreePath{
                    name,
                    artifact.second["time"].as<double>(),
                    std::make_tuple(
                        artifact.second["resolution"]["x"].as<unsigned>(), 
                        artifact.second["resolution"]["y"].as<unsigned>(), 
                        artifact.second["resolution"]["z"].as<unsigned>()
                    )    
                }));
            }
            else if(artifact.first.as<std::string>() == std::string("meanfreetime")) {
                std::string name = "meanfreetime";
                if(artifact.second["name"])
                    out_path+name = artifact.second["name"].as<std::string>();
                artfcs.push_back(std::shared_ptr<SimulationController::Artifact>(new artifacts::MeanFreeTime{
                    name,
                    artifact.second["time"].as<double>(),
                    std::make_tuple(
                        artifact.second["resolution"]["x"].as<unsigned>(), 
                        artifact.second["resolution"]["y"].as<unsigned>(), 
                        artifact.second["resolution"]["z"].as<unsigned>()
                    )    
                }));
            }
            else if(artifact.first.as<std::string>() == std::string("simulatedphonons")) {
                std::string name = "simulatedphonons";
                if(artifact.second["name"])
                    out_path+name = artifact.second["name"].as<std::string>();
                artfcs.push_back(std::shared_ptr<SimulationController::Artifact>(new artifacts::SimulatedPhononCount{
                    name,
                    artifact.second["time"].as<double>(),
                    std::make_tuple(
                        artifact.second["resolution"]["x"].as<unsigned>(), 
                        artifact.second["resolution"]["y"].as<unsigned>(), 
                        artifact.second["resolution"]["z"].as<unsigned>()
                    )    
                }));
            }
            else if(artifact.first.as<std::string>() == std::string("internalenergydensity")) {
                std::string name = "internalenergydensity";
                if(artifact.second["name"])
                    out_path+name = artifact.second["name"].as<std::string>();
                artfcs.push_back(std::shared_ptr<SimulationController::Artifact>(new artifacts::InternalEnergyDensity{
                    name,
                    artifact.second["time"].as<double>(),
                    std::make_tuple(
                        artifact.second["resolution"]["x"].as<unsigned>(), 
                        artifact.second["resolution"]["y"].as<unsigned>(), 
                        artifact.second["resolution"]["z"].as<unsigned>()
                    )    
                }));
            }
            else {
                std::cerr<<"Unkown artifact \""<<artifact.second.as<std::string>()<<"\". Quitting."<<std::endl;
                return -1;
            }
        }
    }
    catch(BoltziException& exc) {
        std::cerr<<exc.what()<<" Quitting."<<std::endl;
        return -1;
    }
    catch(std::exception& exc) {
        // TODO: Print out what name was missing in the config on an exception.
        std::cerr<<exc.what()<<std::endl;
        std::cerr<<"Missing entry in config. Quitting."<<std::endl;
        return -1;
    } 

    bool ignore_init = false;
    double t0 = 0;
    bool time_set = false;

    for(int arg_idx = 2; arg_idx<argc; ++arg_idx) {
        if(strcmp(argv[arg_idx], "-h")==0 || strcmp(argv[arg_idx], "--help")==0) help_msg();
        else if(strcmp(argv[arg_idx], "--silent")==0) Silencer::set(true);
        else if(strcmp(argv[arg_idx], "--less_output")==0) Silencer::set_less(true);
        else if(strcmp(argv[arg_idx], "-a")==0 || strcmp(argv[arg_idx], "--artifact")==0) {
            if(!((arg_idx+2)<argc)) {
                std::cerr<<"Option \""<<argv[arg_idx-1]<<"\" without value(s). Quitting."<<std::endl;
                return -1;
            }
            ++arg_idx;
            char* name_str = argv[arg_idx];
            ++arg_idx;
            double time; 
            try {
                if(std::strlen(argv[arg_idx])==0) throw std::invalid_argument("");
                time = std::stod(argv[arg_idx]);
            }
            catch(std::invalid_argument) {
                std::cerr<<"Can't interpret time \""<<argv[arg_idx]<<"\". Quitting."<<std::endl;
                return -1;
            }
            std::string filename;
            if(arg_idx+1<argc && argv[arg_idx+1][0] != '-') {
                ++arg_idx;
                filename = argv[arg_idx]; 
            }
            if(strcmp(name_str, "temperature")==0) {
                if(filename.empty()) filename = "temperature";
                artfcs.push_back(std::shared_ptr<SimulationController::Artifact>(new artifacts::Temperature{
                    out_path+filename,
                    time,
                    current_resolution
                }));
            }
            else if(strcmp(name_str, "heatcurrent")==0) {
                if(filename.empty()) filename = "heatcurrent";
                artfcs.push_back(std::shared_ptr<SimulationController::Artifact>(new artifacts::HeatCurrent{
                    out_path+filename,
                    time,
                    current_resolution
                }));
            }
            else if(strcmp(name_str, "meanfreepath")==0) {
                if(filename.empty()) filename = "meanfreepath";
                artfcs.push_back(std::shared_ptr<SimulationController::Artifact>(new artifacts::MeanFreePath{
                    out_path+filename,
                    time,
                    current_resolution
                }));
            }
            else if(strcmp(name_str, "meanfreetime")==0) {
                if(filename.empty()) filename = "meanfreetime";
                artfcs.push_back(std::shared_ptr<SimulationController::Artifact>(new artifacts::MeanFreeTime{
                    out_path+filename,
                    time,
                    current_resolution
                }));
            }
            else if(strcmp(name_str, "simulatedphonons")==0) {
                if(filename.empty()) filename = "simulatedphonons";
                artfcs.push_back(std::shared_ptr<SimulationController::Artifact>(new artifacts::SimulatedPhononCount{
                    out_path+filename,
                    time,
                    current_resolution
                }));
            }
            else if(strcmp(name_str, "internalenergydensity")==0) {
                if(filename.empty()) filename = "internalenergydensity";
                artfcs.push_back(std::shared_ptr<SimulationController::Artifact>(new artifacts::InternalEnergyDensity{
                    out_path+filename,
                    time,
                    current_resolution
                }));
            }
            else {
                std::cerr<<"Unkown artifact \""<<name_str<<"\". Quitting."<<std::endl;
                return -1;
            }
        }
        else if(strcmp(argv[arg_idx], "-s")==0 || strcmp(argv[arg_idx], "--save")==0) {
            if(!((arg_idx+1)<argc)) {
                std::cerr<<"Option \""<<argv[arg_idx-1]<<"\" without value(s). Quitting."<<std::endl;
                return -1;
            }
            ++arg_idx;
            double time; 
            try {
                if(std::strlen(argv[arg_idx])==0) throw std::invalid_argument("");
                time = std::stod(argv[arg_idx]);
            }
            catch(std::invalid_argument&) {
                std::cerr<<"Can't interpret time \""<<argv[arg_idx]<<"\". Quitting."<<std::endl;
                return -1;
            }
            std::string filename;
            if(arg_idx+1<argc && argv[arg_idx+1][0] != '-') {
                ++arg_idx;
                filename = argv[arg_idx]; 
            }
            else 
                filename = "state";
            simulation->add_savepoint(time, out_path+filename);
        }
        else if(strcmp(argv[arg_idx], "-l")==0 || strcmp(argv[arg_idx], "--load")==0) {
            if(ignore_init)
                throw BoltziException((std::stringstream{}<<"Multiple uses of \""<<argv[arg_idx]<<"\" make no sense.").str().c_str());
            if(!((arg_idx+1)<argc)) {
                std::cerr<<"Option \""<<argv[arg_idx-1]<<"\" without value(s). Quitting."<<std::endl;
                return -1;
            }
            ++arg_idx;
            std::string fname = argv[arg_idx]; 
            if(fname.length()==0) throw std::invalid_argument("");
            Silencer::print(std::stringstream{}<<"Loading initial phonons from file "<<fname<<".\n");
            sample->phonons.load(fname);
            ignore_init = true;
        }
        else if(strcmp(argv[arg_idx], "-C")==0 || strcmp(argv[arg_idx], "--directory")==0) {
            if(!((arg_idx+1)<argc)) {
                std::cerr<<"Option \""<<argv[arg_idx-1]<<"\" without value(s). Quitting."<<std::endl;
                return -1;
            }
            ++arg_idx;
            out_path = argv[arg_idx]; 
            if(out_path.back() != '/')
                out_path+='/';
            if(out_path.length()==0) throw std::invalid_argument("");
        }
        else if(strcmp(argv[arg_idx], "-b")==0) {
            if(!((arg_idx+1)<argc)) {
                std::cerr<<"Option \""<<argv[arg_idx-1]<<"\" without value(s). Quitting."<<std::endl;
                return -1;
            }
            ++arg_idx;
            boundary_current_log_fname = argv[arg_idx]; 
            if(boundary_current_log_fname.length()==0) throw std::invalid_argument("");
        }
        else if(strcmp(argv[arg_idx], "-t")==0 || strcmp(argv[arg_idx], "--time")==0) {
            if(time_set)
                throw BoltziException((std::stringstream{}<<"Multiple uses of \""<<argv[arg_idx]<<"\" make no sense.").str().c_str());
            time_set = true;
            if(!((arg_idx+1)<argc)) {
                std::cerr<<"Option \""<<argv[arg_idx-1]<<"\" without value(s). Quitting."<<std::endl;
                return -1;
            }
            ++arg_idx;
            double time; 
            try {
                if(std::strlen(argv[arg_idx])==0) throw std::invalid_argument("");
                time = std::stod(argv[arg_idx]);
            }
            catch(std::invalid_argument&) {
                std::cerr<<"Can't interpret time \""<<argv[arg_idx]<<"\". Quitting."<<std::endl;
                return -1;
            }
            t0 = time;
        }
        else if(strcmp(argv[arg_idx], "-r")==0 || strcmp(argv[arg_idx], "--resolution")==0) {
            if(!(arg_idx+3<argc)){
                std::cerr<<"Insufficient amount of parameters to \""<<argv[arg_idx]<<"\". Quitting."<<std::endl;
                return(-1);
            }
            try{
                ++arg_idx;
                if(std::strlen(argv[arg_idx])==0) throw std::invalid_argument("");
                int res = std::stoi(argv[arg_idx]);
                if(res<=0) 
                    throw std::exception();   
                std::get<0>(current_resolution) = static_cast<unsigned>(res); 
                ++arg_idx;
                if(std::strlen(argv[arg_idx])==0) throw std::invalid_argument("");
                res = std::stoi(argv[arg_idx]);
                if(res<=0) 
                    throw std::exception();   
                std::get<1>(current_resolution) = static_cast<unsigned>(res); 
                ++arg_idx;
                if(std::strlen(argv[arg_idx])==0) throw std::invalid_argument("");
                res = std::stoi(argv[arg_idx]);
                if(res<=0) 
                    throw std::exception();   
                std::get<2>(current_resolution) = static_cast<unsigned>(res); 
            }
            catch(...) {
                std::cerr<<"Invalid resolution \""<<argv[arg_idx]<<"\". Quitting."<<std::endl;
                exit(-1);
            }
        }
        else {
            std::cerr<<"Unknown flag/option \""<<argv[arg_idx]<<"\". Quitting."<<std::endl;
            return -1;
        }
    }
    
    simulation->boundary_current_log_fname = out_path+boundary_current_log_fname;

    for(auto artifact : artfcs) {
        simulation->add_artifact(artifact);
    }

    if(!ignore_init) {
        Silencer::print("Spawning initial phonons.\n");
        simulation->initialize(
            [&init_temperature](vec3 r) mutable ->double {return init_temperature(r);}, t0
        );
    }
    Silencer::print("Starting simulation.\n");
    try {
    simulation->run(time, t0);
    }
    catch(BoltziException& exc){
        std::cerr<<exc.what()<<" Quitting."<<std::endl;
        exit(-1);
    }
    Silencer::print("Finished simulation successfully.\n");
    
    util::pool.terminate();

    return 0;
}
