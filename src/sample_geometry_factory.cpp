#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <regex>

#include <OBJ_Loader.h>
#include <armadillo>
#include <algorithm>

#include "headers/boltzi_exception.h"
#include "headers/util/containers.h"
#include "headers/sample_geometry_factory.h"
#include "headers/boundary.h"
#include "headers/boundaries/simple_isothermal_boundary.h"
#include "headers/boundaries/specular_vacuum_boundary.h"
#include "headers/boundaries/bk_vacuum_boundary.h"
#include "headers/boundaries/lambert_vacuum_boundary.h"
#include "headers/boundaries/even_vacuum_boundary.h"
#include "headers/boundaries/measure_boundary.h"

namespace {
std::string parentheses{
    R"(\(((?:(?:\d+.\d*|\d*.\d+|\d)(?:e-?\d+)?\s*,\s*)*(?:(?:\d+.\d*|\d*.\d+|\d)(?:e-?\d+)?\s*)?)\))"
};
std::regex vacuum_boundary{std::string{"vacuum_specular"}+parentheses};
std::regex rough_bk{std::string{"rough_bk"}+parentheses};
std::regex rough_lambert{std::string{"rough_lambert"}+parentheses};
std::regex rough_even{std::string{"rough_even"}+parentheses};
std::regex isothermal_boundary{std::string{"isothermal_simple"}+parentheses};
std::regex measure_boundary{std::string{"measure"}+parentheses};

std::string next_param(std::string& s) {
    auto pos = std::find(s.begin(), s.end(), ',');
    std::string tmp;
    if(pos == s.end()) { 
        tmp = s;
        s.clear();
        return tmp;
    }
    std::string res;
    std::copy(pos+1, s.end(), tmp.begin());
    std::copy(s.begin(), pos, res.begin());
    s = tmp;
    return res; 
}
}

namespace boltzi{
SampleGeometry load_geometry(std::string fname, double scale) {
	objl::Loader loader;

	bool loadout = loader.LoadFile(fname.c_str());
    if (!loadout) throw BoltziException("SampleGeometry file not found.");

    std::vector<std::shared_ptr<Boundary>> bounds;

    for (int i = 0; i < loader.LoadedMeshes.size(); i++) {
        //each mesh becomes one Boundary
        objl::Mesh cur_mesh = loader.LoadedMeshes[i];

        std::vector<Boundary::Primitive> primitives;
        std::string bound_name{cur_mesh.MeshName};

        for (int j = 0; j < cur_mesh.Indices.size(); j += 3) {
            auto v0 = cur_mesh.Vertices[cur_mesh.Indices[j]];
            auto v1 = cur_mesh.Vertices[cur_mesh.Indices[j+1]];
            auto v2 = cur_mesh.Vertices[cur_mesh.Indices[j+2]];
            auto normal = cur_mesh.Vertices[cur_mesh.Indices[j]].Normal;

            vec3 p_normal = arma::normalise(vec3{normal.X, normal.Y, normal.Z});  // Assume all have the same normal.
            auto e_1 = v1.Position - v0.Position;
            vec3 p_e_1 = arma::normalise(vec3{e_1.X,e_1.Y,e_1.Z});
            vec3 p_e_2 = arma::cross(p_normal, p_e_1);

            util::Triangle p_tri{
                vec3{v0.Position.X, v0.Position.Y, v0.Position.Z}*scale,
                vec3{v1.Position.X, v1.Position.Y, v1.Position.Z}*scale,
                vec3{v2.Position.X, v2.Position.Y, v2.Position.Z}*scale
            };

            primitives.emplace_back(
                p_tri,
                p_normal,
                p_e_1,
                p_e_2
            );
        }
        
        // at least 5 matches is necessary since regex doesn't ignore non capturing groups for some reason

        // is SpecularVaccumBoundary 
        if(std::smatch m;std::regex_search(bound_name, m, vacuum_boundary)) {
            if(m.size()!=2 || m[1].str().length()!=0) throw BoltziException("Boundary of type vacuum_specular does not take any parameters.");
            bounds.push_back(
                std::shared_ptr<Boundary>(new boundaries::SpecularVacuumBoundary{primitives})
            );
        }
        // is SimpleIsothermalBoundary 
        else if(std::smatch m;std::regex_search(bound_name, m, isothermal_boundary)) {
            if(m.size()!=2) throw BoltziException("Boundary of type isothermal_simple does take exactly one parameter (temperature).");
            std::string parameters = m[1].str();
            double param = std::stod(next_param(parameters));
            auto next = next_param(parameters);
            if(next.length() != 0) 
                throw BoltziException("Boundary of type isothermal_simple does take exactly one parameter (temperature).");
            bounds.push_back(
                std::shared_ptr<Boundary>(new boundaries::SimpleIsothermalBoundary{primitives, param})
            );
        }
        else if(std::smatch m;std::regex_search(bound_name, m, rough_bk)) {
            if(m.size()!=2) throw BoltziException("Boundary of type rough_bk does take exactly one parameter (standard deviation).");
            std::string parameters = m[1].str();
            double param = std::stod(next_param(parameters));
            auto next = next_param(parameters);
            if(next.length() != 0) 
                throw BoltziException("Boundary of type rough_bk does take exactly one parameter (standard deviation).");
            bounds.push_back(
                std::shared_ptr<Boundary>(new boundaries::BKVacuumBoundary{primitives, param})
            );
        }
        else if(std::smatch m;std::regex_search(bound_name, m, rough_lambert)) {
            if(m.size()!=2 || m[1].str().length()!=0) throw BoltziException("Boundary of type rough_lambert does not take any parameters.");
            bounds.push_back(
                std::shared_ptr<Boundary>(new boundaries::LambertVacuumBoundary{primitives})
            );
        }
        else if(std::smatch m;std::regex_search(bound_name, m, rough_even)) {
            if(m.size()!=2 || m[1].str().length()!=0) throw BoltziException("Boundary of type rough_even does not take any parameters.");
            bounds.push_back(
                std::shared_ptr<Boundary>(new boundaries::EvenVacuumBoundary{primitives})
            );
        }
        else if(std::smatch m;std::regex_search(bound_name, m, measure_boundary)) {
            if(m.size()!=2 || m[1].str().length()!=0) throw BoltziException("Boundary of type measure does not take any parameters.");
            bounds.push_back(
                std::shared_ptr<Boundary>(new boundaries::MeasureBoundary{primitives})
            );
        }
        else {
            throw BoltziException("Unknown or unsupported boundary type found.");
        }
    }
	return SampleGeometry{bounds, scale};
}
}
