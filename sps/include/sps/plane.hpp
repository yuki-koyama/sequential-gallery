#ifndef SPS_PLANE_HPP
#define SPS_PLANE_HPP

#include <Eigen/Core>
#include <array>
#include <utility>
#include <vector>

namespace sps
{
    class AbstractPlane
    {
    public:
        /// \brief Index for grid cells on a plane.
        ///
        /// \details (0, 0) is the center.
        using GridCellIndex = std::pair<int, int>;

        const Eigen::VectorXd&                GetCenter() const { return m_center; }
        const std::array<Eigen::VectorXd, 4>& GetVertices() const { return m_vertices; }

        /// \brief Calculate N-d parameters corresponding to the specified grid coordinates (i.e., [0, 1]^2).
        Eigen::VectorXd CalcParameters(const Eigen::Vector2d& coords) const;

        /// \brief Calculate N-d parameters corresponding to the specified settings and user interactions.
        Eigen::VectorXd CalcGridParameters(const GridCellIndex&              grid_cell,
                                           const unsigned                    num_candidates,
                                           const double                      inter_level_scale,
                                           const std::vector<GridCellIndex>& prev_grid_cells) const;

        double CalcArea() const;

    protected:
        Eigen::VectorXd                m_center;
        std::array<Eigen::VectorXd, 4> m_vertices;
    };

    class SimplePlane : public AbstractPlane
    {
    public:
        SimplePlane(const Eigen::VectorXd& x_base, const Eigen::VectorXd& x_opposite, const Eigen::VectorXd& cross_vec)
        {
            m_center      = 0.5 * (x_base + x_opposite);
            m_vertices[0] = x_base;
            m_vertices[1] = m_center + cross_vec;
            m_vertices[2] = x_opposite;
            m_vertices[3] = m_center - cross_vec;
        }

        const Eigen::VectorXd& GetBaseVertex() const { return m_vertices[0]; }
        Eigen::VectorXd        GetCrossVec() const { return m_vertices[1] - m_center; }
    };

    class PiecewisePlanarPlane : public AbstractPlane
    {
    public:
        PiecewisePlanarPlane(const Eigen::VectorXd& center,
                             const Eigen::VectorXd& vertex_0,
                             const Eigen::VectorXd& vertex_1,
                             const Eigen::VectorXd& vertex_2,
                             const Eigen::VectorXd& vertex_3)
        {
            m_center      = center;
            m_vertices[0] = vertex_0;
            m_vertices[1] = vertex_1;
            m_vertices[2] = vertex_2;
            m_vertices[3] = vertex_3;
        }
    };
} // namespace sps

#endif // SPS_PLANE_HPP
