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
    std::vector<front::Face*> m_frontFaces; // TODO: use std::unordered_set
    std::vector<front::Edge*> m_front_edges; // TODO: use std::set sorting by angle (first larger) instead

    bool shell_contains(const pmg::Vert* vert) const;

    shell::Edge* find_shell_edge(const shell::Vert* v0, const shell::Vert* v1) const;
    front::Face* findFrontFace(const pmg::Face* face) const;
    std::vector<front::Edge*> findFEdge(const pmg::Vert* v0, const pmg::Vert* v1) const;
    std::vector<front::Edge*> findFEdge(const pmg::Edge* edge) const;

    // Adds new front Face and corresponding Face.
    front::Face* add_to_front(const pmg::Face* face, bool addInner = true);
    front::Edge* add_to_front(const pmg::Edge* edge, bool addInner = true);

    void remove_from_front(front::Face* fFace);
    void remove_from_front(front::Edge* front_edge);

    // This section is about various front intersection checks
    // TODO: move that and other methods to relations.h later
    // TODO: change methods names
    bool segmentIntersectMesh(const vec3& v0, const vec3& v1) const;
    bool segmentIntersectFront(const vec3& v0, const vec3& v1) const;
    bool edgeIntersectFront(const pmg::Vert* v0, const vec3& v1) const;
    bool edgeIntersectFront(const pmg::Vert* v0, const pmg::Vert* v1) const;
    bool edgeIntersectAnyFace(const pmg::Edge* edge) const;
    bool potentialEdgeIntersectFront(front::Edge* front_edge) const;
    bool anyEdgeIntersectFace(const pmg::Vert* v0, const pmg::Vert* v1, const vec3& v2) const;
    bool anyEdgeIntersectFace(const pmg::Vert* v0, const pmg::Vert* v1, const pmg::Vert* v2) const;
    bool anyEdgeIntersectPotentialFaces(front::Edge* front_edge) const;
    bool anyVertInsidePotentialTetrCheck(front::Edge* front_edge) const;
    bool parallelFacesCheck(front::Edge* front_edge) const;
    // TODO: add more proximity checks
    bool doesFrontIntersectSphere(const vec3& center, real_t radius) const;
    bool frontSplitCheck(front::Edge* front_edge, front::Edge* oppFEdge = nullptr) const;
    bool frontCollapseCheck(front::Edge* front_edge, front::Edge* oppFEdge = nullptr) const;

    static pair_rr computeMinMaxEdgesLengths(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3);
    static pair_rr computeMinMaxEdgesSqrLengths(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3);
    static real_t  computeTetrSimpleQuality(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3);
    static real_t  computeTetrSimpleSqrQuality(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3);

    front::Edge* currentFrontEdge(real_t maxCompl) const; // TODO: use std::map::upper_bound(...) method instead
    bool exhaustWithoutNewVertPriorityPredicate(front::Edge* front_edge);
    bool exhaustWithNewVertPriorityPredicate(front::Edge* front_edge);
    ExhaustType computeExhaustionTypeQualityPriority(
        front::Edge* front_edge,
        front::Face*& out_withNWFFace, vec3*& out_withNWNewVertPos);

    vec3 computeNormalInTetr(const front::Face* fFace, const vec3& opp_vert_pos) const;
    vec3 computeNormalInTetr(const front::Face* fFace, const pmg::Edge* one_of_remaining_edges) const;
    vec3 computeNormalInTetr(const vec3& fFacePos0, const vec3& fFacePos1, const vec3& fFacePos2, const vec3& opp_vert_pos) const;

    // TODO: add known neighbors initializations while creating new Tetr
    void setFEdgesInFrontSplit(const front::Edge* front_edge, std::array<front::Edge*, 2> newOppFEdges, std::array<front::Face*, 2> newFFaces, pair_ff oppFFaces) const;
    void exhaustFrontCollapse(front::Edge* front_edge, front::Edge* oppFEdge);
    void exhaustFrontSplit(front::Edge* front_edge, front::Edge* oppFEdge);
    void exhaustWithoutNewVertOppEdgeExists(front::Edge* front_edge, front::Edge* oppFEdge);
    void exhaustWithoutNewVertOppEdgeDontExists(front::Edge* front_edge);
    void exhaustWithoutNewVert(front::Edge* front_edge, bool doesOppEdgeExists = true, front::Edge* oppFEdge = nullptr);


    bool tryComputeNewVertPosType3(
        front::Face* fFace, vec3& out_pos);
    bool tryComputeNewVertPosType2(
        front::Face* frontFace, vec3& out_pos,
        std::size_t smallAngleIdx0, std::size_t smallAngleIdx1);
    bool tryComputeNewVertPosType1(
        front::Face* fFace, vec3& out_pos,
        std::size_t smallAngleIdx);
    bool tryComputeNewVertPosType0(
        front::Face* fFace, vec3& out_pos);
    bool tryComputeNewVertPos(front::Face* fFace, vec3& out_pos);

    real_t sqr4FaceArea(const front::Face* fFace) const;
    front::Face* chooseFaceForExhaustionWithNewVert(front::Edge* front_edge);
    void         exhaustWithNewVert(front::Face* fFace, const vec3& vertPos);

    bool tryExhaustWithoutNewVert(front::Edge* front_edge);
    bool tryExhaustWithNewVert(front::Edge* front_edge);

    bool globalIntersectionCheck();

    bool isFrontExhausted();
    void processAngles();

    void debug();

    void computeFrontNormals();
    void initializeFFaceFEdges(front::Face* fFace) const;
    void initializeFront();
};

} // namespace pmg
