#include <Eigen/Geometry>
#include <Eigen/LU>
#include <iostream>
#include <numeric>
#include <sps/plane.hpp>
#include <vector>

// This implementation is general but slow because it needs to solve linear systems. It is possible to implement this
// function much faster for specific plane types.
Eigen::VectorXd sps::AbstractPlane::CalcParameters(const Eigen::Vector2d& coords) const
{
    const Eigen::Vector2d              g_c(+0.0, +0.0);
    const std::vector<Eigen::Vector2d> g = {
        {-1.0, -1.0},
        {+1.0, -1.0},
        {+1.0, +1.0},
        {-1.0, +1.0},
    };

    const auto vertex_indices = [&]() -> std::pair<unsigned, unsigned> {
        if (coords(0) >= coords(1) && coords(0) <= -coords(1))
        {
            return {0, 1};
        }
        if (coords(0) >= coords(1) && coords(0) >= -coords(1))
        {
            return {1, 2};
        }
        if (coords(0) <= coords(1) && coords(0) >= -coords(1))
        {
            return {2, 3};
        }
        if (coords(0) <= coords(1) && coords(0) <= -coords(1))
        {
            return {3, 0};
        }

        assert(false);
        return {0, 0};
    }();

    // Prepare to calculate the barycentric coordinates
    const Eigen::Matrix3d R = [&]() {
        Eigen::Matrix3d R = Eigen::Matrix3d::Ones();

        R.block<2, 1>(0, 0) = g_c;
        R.block<2, 1>(0, 1) = g[std::get<0>(vertex_indices)];
        R.block<2, 1>(0, 2) = g[std::get<1>(vertex_indices)];

        return R;
    }();
    const Eigen::Vector3d r(coords(0), coords(1), 1.0);

    // Find the barycentric coordinates
    const Eigen::Vector3d barycentric_coords = Eigen::FullPivLU<Eigen::Matrix3d>(R).solve(r);

    return barycentric_coords(0) * m_center + barycentric_coords(1) * m_vertices[std::get<0>(vertex_indices)] +
           barycentric_coords(2) * m_vertices[std::get<1>(vertex_indices)];
}

Eigen::VectorXd sps::AbstractPlane::CalcGridParameters(const GridCellIndex&              grid_cell,
                                                       const unsigned int                num_candidates,
                                                       const double                      inter_level_scale,
                                                       const std::vector<GridCellIndex>& prev_grid_cells) const
{
    assert(num_candidates % 2 == 1);

    const int    grid_radius         = (num_candidates - 1) / 2;
    const double grid_to_coord_scale = (1.0 / static_cast<double>(grid_radius));

    // Find the coordinates in the top-level grid coordinate system
    auto zoom_transform = Eigen::Affine2d::Identity();
    for (auto& prev_grid_cell : prev_grid_cells)
    {
        const auto prev_unscaled_grid_coords = Eigen::Vector2d(static_cast<double>(std::get<0>(prev_grid_cell)),
                                                               static_cast<double>(std::get<1>(prev_grid_cell)));
        const auto center                    = grid_to_coord_scale * prev_unscaled_grid_coords;
        zoom_transform =
            zoom_transform * Eigen::Translation2d(center) * Eigen::UniformScaling<double>(inter_level_scale);
    }
    const auto unscaled_grid_coords =
        Eigen::Vector2d(static_cast<double>(std::get<0>(grid_cell)), static_cast<double>(std::get<1>(grid_cell)));
    const Eigen::Vector2d finest_level_coords   = grid_to_coord_scale * unscaled_grid_coords;
    const Eigen::Vector2d coarsest_level_coords = zoom_transform * finest_level_coords;

    return CalcParameters(coarsest_level_coords);
}

double sps::AbstractPlane::CalcArea() const
{
    const auto triangle_area_evaluator =
        [](const Eigen::VectorXd& v_0, const Eigen::VectorXd& v_1, const Eigen::VectorXd& v_2) {
            const Eigen::VectorXd r_0      = v_1 - v_0;
            const Eigen::VectorXd r_1      = v_2 - v_0;
            const double          dot_prod = r_0.dot(r_1);

            return 0.5 * std::sqrt(r_0.squaredNorm() * r_1.squaredNorm() - dot_prod * dot_prod);
        };

    std::array<double, 4> areas;

    areas[0] = triangle_area_evaluator(m_center, m_vertices[0], m_vertices[1]);
    areas[1] = triangle_area_evaluator(m_center, m_vertices[1], m_vertices[2]);
    areas[2] = triangle_area_evaluator(m_center, m_vertices[2], m_vertices[3]);
    areas[3] = triangle_area_evaluator(m_center, m_vertices[3], m_vertices[0]);

    return std::accumulate(areas.begin(), areas.end(), 0.0);
}
