//
//  main.cpp
//  cMCF
//
//  Created by 菊池祐作 on 2025/05/16.
//

#include <iostream>
#include <chrono>
#include <thread>
#include "mesh.h"
#include "cmcf.h"
#include "file_io.h"


int main(int argc, const char * argv[])
{
    if (argc != 2) {
        std::cout << "wrong command line argument" << std::endl;
        std::exit(1);
    }
    std::string filename;
    filename = std::string(argv[1]);


    Mesh mesh;
    read_obj(filename, mesh);
    mesh.make_halfedge_list();
    mesh.normalize();

    std::vector<Eigen::MatrixXd> mesh_list;

    auto start = std::chrono::system_clock::now();

    CMCF::compute(mesh, mesh_list);

    auto end = std::chrono::system_clock::now();

    for (int i = 0; i < mesh_list.size(); ++i) {
        write_vtk(filename, mesh, mesh_list[i], i);
    }

    using namespace std::chrono_literals;
    std::cout << "time : " << (end - start) / 1.0s << " s\n" << std::endl;

    return 0;
}
