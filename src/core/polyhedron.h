// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <array>
#include <list>
#include <vector>
#include "front/face.h"
#include "front/edge.h"
#include "shell/face.h"
#include "shell/edge.h"
#include "shell/vert.h"
#include "tetr.h"
#include "face.h"
#include "edge.h"
#include "vert.h"
#include "../helpers/spatial/vec.h"
#include "relations.h"
#include "../real-type.h"
#include "genparams.h"

#include "polyhedral-set.h"

#include "../definitions.h"


namespace pmg {
// refactor class
class Polyhedron { // TODO: make something like pmg::mesher<spt::polytope<3>>, spt::polytope<3> is pmg::polyhedron
    using pair_rr = std::pair<real_t, real_t>;

public:
    real_t preferred_length() const;

    // TODO: remove these methods
    void add_to_shell(const shell::Face* shellFace);
    void add_to_shell(const shell::Edge* shellEdge);
    void add_to_shell(const shell::Vert* shellVert);

    // TODO: remove these methods
    bool shell_contains(const shell::Face* shellFace) const;
    bool shell_contains(const shell::Edge* shellEdge) const;
    bool shell_contains(const shell::Vert* shellVert) const;

    // TODO: remove "inner" in names and return not only inner things but on shell too
    const std::vector<Tetr*>& inner_tetrs() const;
    const std::vector<Face*>& inner_faces() const;
    const std::vector<Vert*>& inner_verts() const;

    // NOTE: for debug
    const std::vector<front::Face*>& front_faces() const;
    const std::vector<front::Edge*>& front_edges() const;

    // Returns minimum and average tetrahedrons quality or absGrad.
    pair_rr analyze_mesh_quality(std::vector<Tetr*>::iterator* out_minQualityTetr = nullptr);
    pair_rr analyze_mesh_abs_grad();
    void tetrahedralize(real_t preferred_length, genparams::Volume gen_params = genparams::Volume());
    void smooth_mesh(std::size_t nIters);

    Polyhedron() {} // TODO: delete default constructor (Polyhedron() = delete;)
    // TODO: Polyhedron( Shell shell ) : Mesher(shell) {};
    Polyhedron(PolyhedralSet* polyhset); // NOTE: for debug
    ~Polyhedron();


private:
    enum class ExhaustType {
        DontExhaust,
        WithoutNewVert,
        WithNewVert
    };

    using pair_ff = std::pair<front::Face*, front::Face*>;
    using vec3 = spt::vec<3, real_t>;

    genparams::Volume m_genparams;

    pair_rr m_meshQuality;
    pair_rr m_meshAbsGrad;

    real_t m_prefLen = static_cast<real_t>(0.0);

    // TODO: create and use separate Mesher (maybe template <std::size_t Dim>)
    // TODO: replace raw pointers with smart pointers and then make benchmark
    PolyhedralSet* m_polyhset = nullptr;

    // TODO: use Shell class instead
    std::vector<shell::Face*> m_shell_faces;
    std::vector<shell::Edge*> m_shell_edges;
    std::vector<shell::Vert*> m_shell_verts;

    // TODO: create and use Volume class instead
    std::vector<pmg::Tetr*> m_inner_tetrs;
    std::vector<pmg::Face*> m_inner_faces;
    std::vector<pmg::Edge*> m_inner_edges;
    std::vector<pmg::Vert*> m_inner_verts;

    // TODO: create and use Front class instead
    // TODO: add front::Vert class
    std::vector<front::Face*> m_front_faces; // TODO: use std::unordered_set
    std::vector<front::Edge*> m_front_edges; // TODO: use std::set sorting by angle (first larger) instead

    bool shell_contains(const pmg::Vert* vert) const;

    shell::Edge* find_shell_edge(const shell::Vert* v0, const shell::Vert* v1) const;
    front::Face* find_front_face(const pmg::Face* face) const;
    std::vector<front::Edge*> find_front_edge(const pmg::Vert* v0, const pmg::Vert* v1) const;
    std::vector<front::Edge*> find_front_edge(const pmg::Edge* edge) const;

    // Adds new front Face and corresponding Face.
    front::Face* add_to_front(const pmg::Face* face, bool add_inner = true);
    front::Edge* add_to_front(const pmg::Edge* edge, bool add_inner = true);

    void remove_from_front(front::Face* fface);
    void remove_from_front(front::Edge* fedge);

    // This section is about various front intersection checks
    // TODO: move that and other methods to relations.h later
    // TODO: change methods names
    bool segmentIntersectMesh(const vec3& v0, const vec3& v1) const;
    bool segmentIntersectFront(const vec3& v0, const vec3& v1) const;
    bool edge_intersect_front(const pmg::Vert* v0, const vec3& v1) const;
    bool edge_intersect_front(const pmg::Vert* v0, const pmg::Vert* v1) const;
    bool edgeIntersectAnyFace(const pmg::Edge* edge) const;
    bool potentialEdgeIntersectFront(front::Edge* fedge) const;
    bool any_edge_intersect_face(const pmg::Vert* v0, const pmg::Vert* v1, const vec3& v2) const;
    bool any_edge_intersect_face(const pmg::Vert* v0, const pmg::Vert* v1, const pmg::Vert* v2) const;
    bool anyEdgeIntersectPotentialFaces(front::Edge* fedge) const;
    bool anyVertInsidePotentialTetrCheck(front::Edge* fedge) const;
    bool parallelFacesCheck(front::Edge* fedge) const;
    // TODO: add more proximity checks
    bool does_front_intersect_sphere(const vec3& center, real_t radius) const;
    bool frontSplitCheck(front::Edge* fedge, front::Edge* opp_fedge = nullptr) const;
    bool frontCollapseCheck(front::Edge* fedge, front::Edge* opp_fedge = nullptr) const;

    static pair_rr min_max_edges_lengths(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3);
    static pair_rr min_max_edges_sqr_lengths(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3);
    static real_t  tetr_simple_quality(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3);
    static real_t  tetr_simple_sqr_quality(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3);

    front::Edge* current_front_edge(real_t maxCompl) const; // TODO: use std::map::upper_bound(...) method instead
    bool exhaust_without_new_vert_priority_predicate(front::Edge* fedge);
    bool exhaust_with_new_vert_priority_predicate(front::Edge* fedge);
    ExhaustType exhaustion_type_quality_priority(
        front::Edge* fedge,
        front::Face*& out_withNWFFace, vec3*& out_with_nv_new_vert_pos);

    vec3 computeNormalInTetr(const front::Face* fface, const vec3& opp_vert_pos) const;
    vec3 computeNormalInTetr(const front::Face* fface, const pmg::Edge* one_of_remaining_edges) const;
    vec3 computeNormalInTetr(const vec3& fFacePos0, const vec3& fFacePos1, const vec3& fFacePos2, const vec3& opp_vert_pos) const;

    // TODO: add known neighbors initializations while creating new Tetr
    void setFEdgesInFrontSplit(const front::Edge* fedge, std::array<front::Edge*, 2> newOppFEdges, std::array<front::Face*, 2> newFFaces, pair_ff oppFFaces) const;
    void exhaustFrontCollapse(front::Edge* fedge, front::Edge* opp_fedge);
    void exhaustFrontSplit(front::Edge* fedge, front::Edge* opp_fedge);
    void exhaustWithoutNewVertOppEdgeExists(front::Edge* fedge, front::Edge* opp_fedge);
    void exhaustWithoutNewVertOppEdgeDontExists(front::Edge* fedge);
    void exhaust_without_new_vert(front::Edge* fedge, bool doesOppEdgeExists = true, front::Edge* opp_fedge = nullptr);


    bool tryComputeNewVertPosType3(
        front::Face* fface, vec3& out_pos);
    bool try_compute_new_vert_pos_type2(
        front::Face* frontFace, vec3& out_pos,
        std::size_t smallAngleIdx0, std::size_t smallAngleIdx1);
    bool try_compute_new_vert_pos_type1(
        front::Face* fface, vec3& out_pos,
        std::size_t smallAngleIdx);
    bool try_compute_new_vert_pos_type0(
        front::Face* fface, vec3& out_pos);
    bool try_compute_new_vert_pos(front::Face* fface, vec3& out_pos);

    real_t sqr4FaceArea(const front::Face* fface) const;
    front::Face* chooseFaceForExhaustionWithNewVert(front::Edge* fedge);
    void         exhaust_with_new_vert(front::Face* fface, const vec3& vertPos);

    bool try_exhaust_without_new_vert(front::Edge* fedge);
    bool try_exhaust_with_new_vert(front::Edge* fedge);

    bool global_intersection();

    bool isFrontExhausted();
    void process_angles();

    void debug();

    void compute_front_normals();
    void initializeFFaceFEdges(front::Face* fface) const;
    void initialize_front();
};

} // namespace pmg
