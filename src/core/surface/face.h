// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <array>
#include <list>
#include <vector>
#include "../edge.h"
#include "../vert.h"
#include "front/edge.h"
#include "front/vert.h"
#include "../../real-type.h"
#include "../genparams.h"

#include "../../definitions.h"


namespace pmg::surface {

// TODO: separate Mesher<2>
class Face {
public:
    // TODO: maybe use std::reference_wrapper instead of pointer
    // TODO: represent surface::Face as any number of surface::Edges instead
    std::array<surface::Edge*, 3> edges;

    real_t preferred_length() const;

    const std::list<pmg::Face*>& inner_faces() const;
    const std::list<pmg::Edge*>& inner_edges() const;
    const std::list<pmg::Vert*>& inner_verts() const;

    const std::list<front::Edge*>& front_edges() const;

    surface::Vert* find_vert_not(const surface::Edge* sEdge) const;
    surface::Edge* find_surface_edge_containing(const pmg::Edge* edge) const;

    // Needs surface edges to be already segmentized.
    void triangulate(real_t preferredLen, genparams::Surface gen_params = genparams::Surface());
    void smooth_mesh(std::size_t nIters);
    void delaunay_postp();
    void optimize_mesh(std::size_t nSmoothIters = 20, std::size_t nDelaunaySmoothIters = 3);

    bool contains(const surface::Edge* sEdge) const;
    bool contains(const surface::Vert* sVert) const;

    Face(const surface::Edge* sEdge0, const surface::Edge* sEdge1, const surface::Edge* sEdge2);


private:
    enum class ExhaustType {
        DontExhaust,
        WithoutNewVert,
        WithNewVert
    };

    using pair_rr = std::pair<real_t, real_t>;
    using pair_ff = std::pair<pmg::Face*, pmg::Face*>;
    using pair_ee = std::pair<pmg::Edge*, pmg::Edge*>;
    using vec3 = spt::vec<3, real_t>;

    real_t m_prefLen = static_cast<real_t>(0.0);

    std::list<pmg::Face*> m_inner_faces;
    std::list<pmg::Edge*> m_inner_edges;
    std::list<pmg::Vert*> m_inner_verts;

    std::list<front::Edge*> m_front_edges;
    std::list<front::Vert*> m_front_verts;

    front::Vert* find_front_vert(const pmg::Vert* vert) const;

    front::Edge* add_to_front(const pmg::Edge* edge);
    front::Vert* add_to_front(const pmg::Vert* vert);

    void remove_from_front(front::Edge* front_edge);
    void remove_from_front(front::Vert* fVert);

    bool any_vert_inside_potential_triangle_check(front::Vert* fVert) const;
    bool does_segment_intersects_with_front(const vec3& v0, const vec3& v1) const;
    bool does_segment_intersects_with_front(const pmg::Vert* v0, const vec3& v1) const;

    vec3 normal_in_triangle(front::Edge* front_edge, const vec3& opp_vert_pos);
    vec3 normal_in_triangle(front::Edge* front_edge, pmg::Edge* one_of_remaining_edges); // Do i need it?

    bool tryComputeNewVertPosType2(front::Edge* front_edge, vec3& out_pos);
    bool tryComputeNewVertPosType1(front::Edge* front_edge, vec3& out_pos, std::size_t smallAngleIdx);
    bool tryComputeNewVertPosType0(front::Edge* front_edge, vec3& out_pos);
    bool tryComputeNewVertPos(front::Edge* front_edge, vec3& out_pos);

    static pair_rr computeMinMaxEdgesLengths(const vec3& p0, const vec3& p1, const vec3& p2);
    static pair_rr computeMinMaxEdgesSqrLengths(const vec3& p0, const vec3& p1, const vec3& p2);
    static real_t  computeTriangSimpleQuality(const vec3& p0, const vec3& p1, const vec3& p2);
    static real_t  computeTriangSimpleSqrQuality(const vec3& p0, const vec3& p1, const vec3& p2);

    front::Edge* chooseEdgeForExhaustionWithNewVert(front::Vert* fVert);
    void exhaustWithNewVert(front::Edge* front_edge, const vec3& vertPos);
    void exhaustWithoutNewVert(front::Vert* fVert);

    bool tryExhaustWithoutNewVert(front::Vert* fVert);
    bool tryExhaustWithNewVert(front::Vert* fVert);

    bool globalIntersectionCheck() const;

    front::Vert* currentFrontVert(real_t maxCompl) const;
    bool exhaustWithoutNewVertPriorityPredicate(front::Vert* front_edge);
    bool exhaustWithNewVertPriorityPredicate(front::Vert* front_edge);
    ExhaustType computeExhaustionTypeQualityPriority(
        front::Vert* fVert,
        front::Edge*& out_withNWFrontEdge, vec3*& out_withNWNewVertPos);

    void processLastFace();
    void processAngles();

    pair_ff find2AdjFaces(pmg::Edge* edge) const;
    bool flipIfNeeded(pmg::Edge* edge);

    void computeFrontNormals() const;
    void initializeFront();
};

} // namespace pmg::surface
