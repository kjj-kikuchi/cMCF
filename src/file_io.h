//
//  file_io.h
//  cMCF
//
//  Created by 菊池祐作 on 2025/05/16.
//

#ifndef file_io_h
#define file_io_h

#include <fstream>
#include "mesh.h"

inline void read_obj(std::string const& filename, Mesh &mesh)
{
    std::ifstream ifs(filename);
    if (ifs.fail()) {
        std::cerr << "Failed to open file." << "\n";
        std::exit(1);
    }
    std::string line;
    std::vector<Eigen::Vector3d> vertices;
    while (std::getline(ifs, line)){
        if (line.empty() || line[0] == '#') {
            // Do nothing
        }
        else if (line[0] == 'v' && line[1] != 'n') {
            Eigen::Vector3d v;
            std::istringstream string_in{line.substr(1)};
            string_in >> v(0) >> v(1) >> v(2);
            vertices.push_back(v);
        }
        else if (line[0] == 'v' && line[1] == 'n') {
            Eigen::Vector3d vn;
            std::istringstream string_in{line.substr(2)};
            string_in >> vn(0) >> vn(1) >> vn(2);
            mesh.vertex_normal.push_back(vn);
        }
        else if (line[0] == 'f') {
            Eigen::Vector3i f;
            std::replace(line.begin(), line.end(), '/', ' ');
            std::istringstream string_in{line.substr(1)};
            string_in >> f(0) >> f(0) >> f(1) >> f(1) >> f(2) >> f(2);
            f -= Eigen::Vector3i{1, 1, 1};
            mesh.F.push_back(f);
        }
    }

    mesh.V.resize(vertices.size(), 3);
    for (int i = 0; i < vertices.size(); ++i) {
        mesh.V.row(i) = vertices[i];
    }
}

std::string remove_extension(std::string const& filename) {
    std::filesystem::path p{filename};
    p.replace_extension("");
    return p.string();
}

inline void write_obj(std::string filename, Mesh const& mesh)
{
    std::ofstream of;
    filename = remove_extension(filename) + "_cmfc.obj";

    of.open(filename, std::ios::out);
    for (int i = 0; i < mesh.V.rows(); ++i) {
        of << "v " << mesh.V.row(i) << std::endl;
    }

    for (const auto &face : mesh.F) {
        Eigen::Vector3i f = face + Eigen::Vector3i::Ones();
        of << "f " << f.transpose() << std::endl;
    }
    of.close();
}

inline void write_vtk(std::string filename, Mesh const& mesh, Eigen::MatrixXd const& vertices, int idx)
{
    std::ofstream of;
    filename = remove_extension(filename) + "_cmfc" + std::to_string(idx) +  ".vtk";

    of.open(filename, std::ios::out);
    of << "# vtk DataFile Version 3.0\n" << "vtk output\n" << "ASCII\n" << "DATASET POLYDATA\n";

    of << "POINTS " << vertices.rows() << " float\n";
    of << vertices << std::endl;

    of << "POLYGONS " << mesh.F.size() << " " << 4*mesh.F.size() << std::endl;
    for (const auto &f : mesh.F) {
        of << "3 " << f.transpose() << std::endl;
    }

    of << "POINT_DATA " << vertices.rows() << std::endl;
    of << "NORMALS normal float\n";
//    of << "COLOR_SCALARS normal 3\n";
    for(const auto &n : mesh.vertex_normal) {
        of << n.transpose() << "\n";
    }
//    of << "SCALARS NormalMagnitude float 1\n";
//    of << "LOOKUP_TABLE default\n";
//    for(const auto &n : mesh.vertex_normal) {
//        of << n.norm() << "\n";
//    }
//    of.close();
}


#endif /* file_io_h */
