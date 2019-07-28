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
    enum class exhaust_type {
        dont_exhaust,
        without_new_vert,
        with_new_vert
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
    bool edge_intersect_any_face(const pmg::Edge* edge) const;
    bool will_edge_intersect_front(front::Edge* fedge) const;
    bool any_edge_intersect_face(const pmg::Vert* v0, const pmg::Vert* v1, const vec3& v2) const;
    bool any_edge_intersect_face(const pmg::Vert* v0, const pmg::Vert* v1, const pmg::Vert* v2) const;
    bool will_any_edge_intersect_faces(front::Edge* fedge) const;
    bool will_any_vert_inside_tetr(front::Edge* fedge) const;
    bool will_parallel_faces(front::Edge* fedge) const;
    // TODO: add more proximity checks
    bool does_front_intersect_sphere(const vec3& center, real_t radius) const;
    bool will_front_split(front::Edge* fedge, front::Edge* opp_fedge = nullptr) const;
    bool will_front_collapse(front::Edge* fedge, front::Edge* opp_fedge = nullptr) const;

    static pair_rr min_max_edges_lengths(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3);
    static pair_rr min_max_edges_sqr_lengths(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3);
    static real_t  tetr_simple_quality(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3);
    static real_t  tetr_simple_sqr_quality(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3);

    bool exhaust_without_new_vert_priority_predicate(front::Edge* fedge);
    bool exhaust_with_new_vert_priority_predicate(front::Edge* fedge);
    exhaust_type exhaustion_type_quality_priority(
        front::Edge* fedge,
        front::Face*& out_withNWFFace, vec3*& out_with_nv_new_vert_pos);

    vec3 compute_normal_in_tetr(const front::Face* fface, const vec3& opp_vert_pos) const;
    vec3 compute_normal_in_tetr(const front::Face* fface, const pmg::Edge* one_of_remaining_edges) const;
    vec3 compute_normal_in_tetr(const vec3& fface_pos0, const vec3& fface_pos1, const vec3& fface_pos2, const vec3& opp_vert_pos) const;

    // TODO: add known neighbors initializations while creating new Tetr
    void set_front_edges_in_front_split(
        const front::Edge* fedge, std::array<front::Edge*, 2> new_opp_fedges, 
        std::array<front::Face*, 2> new_ffaces, pair_ff opp_ffaces) const;
    void exhaust_front_collapse(front::Edge* fedge, front::Edge* opp_fedge);
    void exhaust_front_split(front::Edge* fedge, front::Edge* opp_fedge);
    void exhaust_without_new_vert_opp_edge_exists(front::Edge* fedge, front::Edge* opp_fedge);
    void exhaust_without_new_vert_opp_edge_doesnt_exist(front::Edge* fedge);
    void exhaust_without_new_vert(front::Edge* fedge, bool does_opp_edge_exists = true, front::Edge* opp_fedge = nullptr);

    bool try_compute_new_vert_pos_type3(
        front::Face* fface, vec3& out_pos);
    bool try_compute_new_vert_pos_type2(
        front::Face* fface, vec3& out_pos,
        std::size_t small_angle_idx0, std::size_t small_angle_idx1);
    bool try_compute_new_vert_pos_type1(
        front::Face* fface, vec3& out_pos,
        std::size_t small_angle_idx);
    bool try_compute_new_vert_pos_type0(
        front::Face* fface, vec3& out_pos);
    bool try_compute_new_vert_pos(front::Face* fface, vec3& out_pos);

    real_t sqr_face_4areas(const front::Face* fface) const;
    front::Face* choose_face_for_exhaustion_with_new_vert(front::Edge* fedge);
    void         exhaust_with_new_vert(front::Face* fface, const vec3& vert_pos);

    bool try_exhaust_without_new_vert(front::Edge* fedge);
    bool try_exhaust_with_new_vert(front::Edge* fedge);

    bool global_intersection();

    front::Edge* current_front_edge(real_t max_compl) const; // TODO: use std::map::upper_bound(...) method instead
    bool front_exhausted();
    void process_angles();

    void debug();
    real_t exhausted_volume() const;

    void compute_front_normals();
    void initialize_fface_fedges(front::Face* fface) const;
    void initialize_front();
};

} // namespace pmg
