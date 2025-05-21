//
//  cmcf.h
//  cMCF
//
//  Created by 菊池祐作 on 2025/05/16.
//

#ifndef cmcf_h
#define cmcf_h

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/SparseLU>
#include "mesh.h"

namespace CMCF
{
using SparseMatrix = Eigen::SparseMatrix<double>;
using SparseVector = Eigen::SparseVector<double>;
using SparseSolver = Eigen::SparseLU<SparseMatrix>;
using Triplet = Eigen::Triplet<double>;

SparseMatrix cotangent_laplacian(Mesh const& mesh)
{
    SparseMatrix result(mesh.V.rows(), mesh.V.rows());
    std::vector<Triplet> triplets;

    for (int i = 0; i < mesh.V.rows(); i++) {
        int h_temp = mesh.h_out[i];
        int h_end = h_temp;
        double sum_weight = 0;

        do {
            int v_i = i;
            int v_j = mesh.halfedges[h_temp].v_tgt(mesh.F);
            int v_alpha = mesh.halfedges[mesh.h_ccw(h_temp)].v_tgt(mesh.F);
            int v_beta  = mesh.halfedges[mesh.h_cw(h_temp)].v_tgt(mesh.F);

            Eigen::Vector3d e_alphai = mesh.V.row(v_i) - mesh.V.row(v_alpha);
            Eigen::Vector3d e_alphaj = mesh.V.row(v_j) - mesh.V.row(v_alpha);
            Eigen::Vector3d e_betai  = mesh.V.row(v_i) - mesh.V.row(v_beta);
            Eigen::Vector3d e_betaj  = mesh.V.row(v_j) - mesh.V.row(v_beta);

            double alpha = acos(e_alphai.dot(e_alphaj) / (e_alphai.norm()*e_alphaj.norm()));
            double beta = acos(e_betaj.dot(e_betai) / (e_betaj.norm()*e_betai.norm()));

            double cot_weight = (1/tan(alpha) + 1/tan(beta)) / 2.0;
            if(isnan(cot_weight)) std::cout << "error : cotan weight is not correct.\n";
            triplets.push_back(Triplet(v_i, v_j, cot_weight));
            sum_weight += cot_weight;

            h_temp = mesh.h_ccw(h_temp);
        } while (h_temp != h_end && h_temp >= 0);

        triplets.push_back(Triplet(i, i, -sum_weight));
    }
    result.setFromTriplets(triplets.begin(), triplets.end());

    return result;
}

SparseMatrix mass_matrix(Mesh const& mesh)
{
    SparseMatrix result(mesh.V.rows(), mesh.V.rows());
    std::vector<Triplet> triplets;
    std::vector<double> vertex_area(mesh.V.rows(), 0.0);

    for (int i = 0; i < mesh.halfedges.size(); ++i) {
        int v_i = mesh.halfedges[i].v_src(mesh.F);
        int v_j = mesh.halfedges[i].v_tgt(mesh.F);
        double face_area1 = mesh.triangle_area(mesh.h_vec(i), mesh.h_vec(mesh.h_ccw(i)));
        double face_area2 = mesh.triangle_area(mesh.h_vec(i), mesh.h_vec(mesh.h_cw(i)));

        // set the area to the corresponding edge
        triplets.push_back(Triplet(v_i, v_j, (face_area1 + face_area2) / 12.0));

        vertex_area[v_i] += face_area1;
    }

    for (int i = 0; i < mesh.V.rows(); ++i) {
        // set the area to the corresponding vertex
        triplets.push_back(Triplet(i, i, vertex_area[i] / 6.0));
    }
    result.setFromTriplets(triplets.begin(), triplets.end());

    return result;
}

double mean_vertex_difference(Eigen::MatrixXd const& vertices,
                               Eigen::MatrixXd const& new_vertices)
{
    Eigen::VectorXd diff_norm = (new_vertices - vertices).rowwise().norm();
    return diff_norm.mean();
}

void compute(Mesh& mesh, std::vector<Eigen::MatrixXd> &mesh_list)
{
    int iter_num = std::pow(2, 10) + 1;
    int snapshot_index = 1;

    mesh_list.push_back(mesh.V);
    Eigen::MatrixXd new_vertices;

    SparseMatrix cot_laplacian = cotangent_laplacian(mesh);
    double const dt = 1e-3;

    for (int i = 0; i < iter_num; ++i) {
        SparseMatrix mass = mass_matrix(mesh);
        SparseMatrix Lhs = mass - dt * cot_laplacian;
        Eigen::MatrixXd Rhs = mass * mesh.V;

        SparseSolver solver;
        solver.compute(Lhs);
        if (solver.info() != Eigen::Success) {
            std::cout << "decomposition failed\n";
        }

        new_vertices = solver.solve(Rhs);

        mesh.update(new_vertices);
        mesh.normalize();

        if (i == snapshot_index) {
            mesh_list.push_back(new_vertices);
            snapshot_index *= 2;
        }

        std::cout << "iter : " << i << "\n";
    }
}

}

#endif /* cmcf_h */
