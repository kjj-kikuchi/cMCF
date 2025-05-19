//
//  halfedge.h
//  cMCF
//
//  Created by 菊池祐作 on 2025/05/16.
//

#ifndef halfedge_h
#define halfedge_h

struct Halfedge
{
    int idx;
    int h_opp;

    int face() const
    {
        return int(idx / 3);
    }

    int h_next() const
    {
        if(idx % 3 == 2) return idx - 2;
        else return idx + 1;
    }

    int h_prev() const
    {
        if(idx % 3 == 0) return idx + 2;
        else return idx - 1;
    }

    int v_src(std::vector<Eigen::Vector3i> const& f) const
    {
        return f[face()][(idx+1)%3];
    }

    int v_tgt(std::vector<Eigen::Vector3i> const& f) const
    {
        return f[face()][(idx+2)%3];
    }
};

#endif /* halfedge_h */
