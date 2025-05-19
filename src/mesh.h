//
//  mesh.h
//  cMCF
//
//  Created by 菊池祐作 on 2025/05/16.
//

#ifndef mesh_h
#define mesh_h

#include <map>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "halfedge.h"

struct Mesh
{
    Eigen::MatrixXd V;
    std::vector<Eigen::Vector3i> F;
    std::vector<Eigen::Vector3d> vertex_normal;
    std::vector<Halfedge> halfedges;
    std::vector<int> h_out;

    void update(Eigen::MatrixXd const& new_vertices)
    {
        V = new_vertices;
    }

    void make_halfedge_list()
    {
        std::map<std::pair<int, int>, int> h_map;
        for (int i = 0; i < F.size(); i++) {
            for (int j = 0; j < 3; j++) {
                Halfedge h;
                h.idx = 3*i+j;

                auto key = std::make_pair(h.v_src(F), h.v_tgt(F));
                auto keyswap = std::make_pair(h.v_tgt(F), h.v_src(F));
                if (h_map.contains(keyswap)) {
                    h.h_opp = h_map.at(keyswap);
                    halfedges[h_map.at(keyswap)].h_opp = 3*i+j;
                }
                else {
                    h.h_opp = -1;
                    h_map.emplace(key, 3*i+j);
                }
                halfedges.push_back(h);
            }
        }
        // h_out を計算
        h_out.resize(V.rows());
        for (int i = 0; i < halfedges.size(); i++) {
            //  h_out に境界半辺が保存されていない場合のみ更新
            if (halfedges[ h_out[halfedges[i].v_src(F)] ].h_opp != -1) {
                h_out[halfedges[i].v_src(F)] = i;
            }
        }
    }

    int h_ccw(int i) const
    {
        if(i < 0) return i;
        return halfedges[ halfedges[i].h_prev() ].h_opp;
    }

    int h_cw(int i) const
    {
        if (i < 0) return i;
        else if (halfedges[i].h_opp < 0) return halfedges[i].h_opp;
        return halfedges[ halfedges[i].h_opp ].h_next();
    }

    Eigen::Vector3d h_vec(int i) const
    {
        return V.row(halfedges[i].v_tgt(F)) - V.row(halfedges[i].v_src(F));
    }

    double triangle_area(Eigen::Vector3d v0, Eigen::Vector3d v1, Eigen::Vector3d v2) const
    {
        return ((v1-v0).cross(v2-v0)).norm() / 2.0;
    }

    double triangle_area(Eigen::Vector3d e0, Eigen::Vector3d e1) const
    {
        return (e0.cross(e1)).norm() / 2.0;
    }

    double triangle_area(int f_idx)
    {
        Eigen::Vector3d e0 = V.row(F[f_idx](1)) - V.row(F[f_idx](0));
        Eigen::Vector3d e1 = V.row(F[f_idx](2)) - V.row(F[f_idx](0));
        return e0.cross(e1).norm() / 2.0;
    }

    void normalize()
    {
        // 表面積を1に正規化
        double mesh_surface_area = 0.0;
        for (int i = 0; i < F.size(); ++i) {
            mesh_surface_area += triangle_area(i);
        }
        for (int i = 0; i < V.rows(); ++i) {
            V.row(i) /= sqrt(mesh_surface_area);
        }
        // 重心を原点に移動
        Eigen::RowVector3d centroid = V.colwise().mean();
        V.rowwise() -= centroid;
    }
};

#endif /* mesh_h */
