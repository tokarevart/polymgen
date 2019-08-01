// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <list>
#include <vector>
#include "../polyhedron.h"
#include "face.h"
#include "../vert.h"
#include "vert.h"
#include "../relations.h"
#include "../../real-type.h"

#include "../../definitions.h"

// TODO: maybe rename front-edge to volume-front-edge and pmg::front to pmg::volume::front
namespace pmg::front {

class Edge {
    // TODO: pair -> array
    using pair_vv = std::pair<pmg::Vert*, pmg::Vert*>;
    using pair_ff = std::pair<front::Face*, front::Face*>;

public:
    pmg::Edge* x;
    std::array<front::Vert*, 2> front_verts = { nullptr, nullptr };

    void add_front_vert(const front::Vert* fvert);
    void remove_front_vert(const front::Vert* fvert);
    bool front_verts_full() const;

    void refresh_angle_data();

    real_t angle();
    real_t complexity();
    real_t compute_angle();
    real_t compute_complexity();

    // TODO: rewrite these methods based on relations::opposite function
    pair_vv opp_verts();
    // NOTE: maybe store oppEdge and information of its existance
    front::Edge* find_opp_edge();

    // Adj means adjacent.
    pair_ff adj_ffaces();

    bool add_adj_fface(const front::Face* fface);
    bool remove_adj_fface(const front::Face* fface);
    void clear_adj_ffaces();
    bool adj_ffaces_contains(const front::Face* fface) const;
    void fill_adj_ffaces(const front::Face* fFace0, const front::Face* fFace1);
    // TODO: add method std::size_t nAdjFFaces() const;

    bool contains(const front::Vert* fvert) const;

    Edge(const Polyhedron* related_polyhedron, const front::Vert* fvert0, const front::Vert* fvert1);
    Edge(const Polyhedron* related_polyhedron, const pmg::Edge* edge);


private:
    Polyhedron* m_related_polyhedron; // TODO: it's not good to store it
    pair_ff m_adj_ffaces = { nullptr, nullptr };

    real_t m_angle = static_cast<real_t>(0);
    real_t m_complexity = static_cast<real_t>(0);

    bool m_need_angle_processing = true;
    bool m_need_complexity_processing = true;

    bool adj_faces_full();
    // TODO: remove this method and then before calling adj_ffaces() adjacent faces must be assigned manually
    pair_ff fill_adj_ffaces();
};

} // namespace pmg::front
