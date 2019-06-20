// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cmath>
#include <algorithm>
#include "vec.h"

// TODO: experiment with pass by value instead of reference and make benchmark
namespace spt {
       
vec<3> project(
    const vec<3>& point,
    const vec<3>& line_p0, const vec<3>& line_p1);

// TODO: return std::optional<vec<3>>
bool project(
    vec<3>& out,
    const vec<3>& point,
    const vec<3>& segm_p0, const vec<3>& segm_p1);

vec<3> project(
    const vec<3>& point,
    const vec<3>& plane_p0, const vec<3>& plane_p1, const vec<3>& plane_p2);

// TODO: return std::optional<vec<3>>
bool project(
    vec<3>& out,
    const vec<3>& point,
    const vec<3>& trngl_p0, const vec<3>& trngl_p1, const vec<3>& trngl_p2);

bool does_ray_intersect_plane(
    const vec<3>& dir,
    const vec<3>& plane_p0, const vec<3>& plane_p1, const vec<3>& plane_p2);

// TODO: return std::optional<vec<3>>
bool ray_intersect_plane(
    vec<3>& out_intersectPoint,
    const vec<3>& origin, const vec<3>& dir,
    const vec<3>& plane_p0, const vec<3>& plane_p1, const vec<3>& plane_p2);

bool does_ray_intersect_triangle(
    const vec<3>& origin, const vec<3>& dir,
    const vec<3>& trngl_p0, const vec<3>& trngl_p1, const vec<3>& trngl_p2);

// TODO: return std::optional<vec<3>>
bool ray_intersect_triangle(
    vec<3>& out_intersectPoint,
    const vec<3>& origin, const vec<3>& dir,
    const vec<3>& trngl_p0, const vec<3>& trngl_p1, const vec<3>& trngl_p2);

bool does_line_intersect_plane(
    const vec<3>& line_point, const vec<3>& line_dir,
    const vec<3>& plane_p0, const vec<3>& plane_p1, const vec<3>& plane_p2);

// TODO: return std::optional<vec<3>>
bool line_intersect_plane(
    vec<3>& out_intersectPoint,
    const vec<3>& line_point, const vec<3>& line_dir,
    const vec<3>& plane_p0, const vec<3>& plane_p1, const vec<3>& plane_p2);

vec<3> line_intersect_plane(
    const vec<3>& line_point, const vec<3>& line_dir,
    const vec<3>& plane_p0, const vec<3>& plane_p1, const vec<3>& plane_p2);

bool does_segment_intersect_triangle(
    const vec<3>& segm_p0, const vec<3>& segm_p1,
    const vec<3>& trngl_p0, const vec<3>& trngl_p1, const vec<3>& trngl_p2);

// TODO: return std::optional<vec<3>>
bool segment_intersect_triangle(
    vec<3>& out_intersectPoint,
    const vec<3>& segm_p0, const vec<3>& segm_p1,
    const vec<3>& trngl_p0, const vec<3>& trngl_p1, const vec<3>& trngl_p2);

bool does_segment_intersect_plane(
    const vec<3>& segm_p0, const vec<3>& segm_p1,
    const vec<3>& plane_p0, const vec<3>& plane_p1, const vec<3>& plane_p2);

// TODO: return std::optional<vec<3>>
bool segment_intersect_plane(
    vec<3>& out_intersectPoint,
    const vec<3>& segm_p0, const vec<3>& segm_p1,
    const vec<3>& plane_p0, const vec<3>& plane_p1, const vec<3>& plane_p2);

vec<3> segment_intersect_plane(
    const vec<3>& segm_p0, const vec<3>& segm_p1,
    const vec<3>& plane_p0, const vec<3>& plane_p1, const vec<3>& plane_p2);

bool does_triangle_intersect_sphere(
    const vec<3>& trngl_p0, const vec<3>& trngl_p1, const vec<3>& trngl_p2,
    const vec<3>& center, vec<3>::real_type radius);

vec<3>::real_type sqrs_sum(
    const vec<3>& point,
    const vec<3>& trngl_p0, const vec<3>& trngl_p1, const vec<3>& trngl_p2);

vec<3>::real_type max_sqrs_sum(
    const vec<3>& trngl_p0, const vec<3>& trngl_p1, const vec<3>& trngl_p2);

bool is_point_on_plane(
    const vec<3>& point,
    const vec<3>& plane_p0, const vec<3>& plane_p1, const vec<3>& plane_p2);

bool is_point_on_triangle(
    const vec<3>& point,
    const vec<3>& trngl_p0, const vec<3>& trngl_p1, const vec<3>& trngl_p2);

bool is_point_on_triangle(
    const vec<3>& point,
    const vec<3>& trngl_p0, const vec<3>& trngl_p1, const vec<3>& trngl_p2,
    vec<3>::real_type max_sqrs_sum);

bool is_point_in_tetrahedron(
    const vec<3>& point,
    const vec<3>& tetr_p0, const vec<3>& tetr_p1, const vec<3>& tetr_p2, const vec<3>& tetr_p3);

vec<3> closest_segment_point_to_point(
    const vec<3>& point,
    const vec<3>& segm_p0, const vec<3>& segm_p1);

vec<3> closest_triangle_point_to_point_on_plane(
    const vec<3>& point,
    const vec<3>& trngl_p0, const vec<3>& trngl_p1, const vec<3>& trngl_p2);

vec<3>::real_type distance_point_to_line(
    const vec<3>& point,
    const vec<3>& line_p0, const vec<3>& line_p1);

vec<3>::real_type distance_point_to_segment(
    const vec<3>& point,
    const vec<3>& segm_p0, const vec<3>& segm_p1);

vec<3>::real_type distance_point_to_plane(
    const vec<3>& point,
    const vec<3>& plane_p0, const vec<3>& plane_p1, const vec<3>& plane_p2);

vec<3>::real_type distance_point_to_triangle(
    const vec<3>& point,
    const vec<3>& trngl_p0, const vec<3>& trngl_p1, const vec<3>& trngl_p2);

vec<3>::real_type distance_point_to_triangle_on_plane(
    const vec<3>& point,
    const vec<3>& trngl_p0, const vec<3>& trngl_p1, const vec<3>& trngl_p2);

vec<3> lines_closest_point(
    const vec<3>& line0_p0, const vec<3>& line0_p1,
    const vec<3>& line1_p0, const vec<3>& line1_p1);

vec<3>::real_type lines_distance(
    const vec<3>& line0_p0, const vec<3>& line0_p1,
    const vec<3>& line1_p0, const vec<3>& line1_p1);

vec<3>::real_type segments_distance(
    const vec<3>& segm0_p0, const vec<3>& segm0_p1,
    const vec<3>& segm1_p0, const vec<3>& segm1_p1);

// CPA - Closest Point of Approach.
vec<3>::real_type cpa_time(
    const vec<3>& start0, const vec<3>& vel0,
    const vec<3>& start1, const vec<3>& vel1);

vec<3>::real_type cpa_distance(
    const vec<3>& start0, const vec<3>& vel0,
    const vec<3>& start1, const vec<3>& vel1);

} // namespace spt::algs
