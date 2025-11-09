#include <gtest/gtest.h>
#include <sps/plane.hpp>
#include <sps/sps.hpp>

TEST(SpsTest, RandomSubmitTest)
{
    constexpr int num_dims = 6;

    auto optimizer = sps::Optimizer(num_dims);

    optimizer.SubmitData(Eigen::VectorXd::Random(num_dims),
                         {Eigen::VectorXd::Random(num_dims), Eigen::VectorXd::Random(num_dims)});

    const auto x = Eigen::VectorXd::Random(num_dims);

    const auto mean = optimizer.GetCurrentPredictiveMean(x);
    const auto var  = optimizer.GetCurrentPredictiveVar(x);

    EXPECT_FALSE(std::isnan(mean));
    EXPECT_FALSE(std::isnan(var));
    EXPECT_TRUE(var >= 0.0);

    const auto acquisition = optimizer.GetCurrentExpectedImprovement(x);

    EXPECT_TRUE(acquisition >= 0.0);

    const auto plane = optimizer.RetrieveSearchPlane();

    EXPECT_TRUE(plane != nullptr);
    EXPECT_TRUE(plane->CalcArea() >= 0.0);

    const auto x_rand = plane->CalcParameters(Eigen::Vector2d::Random());

    EXPECT_TRUE(x_rand.size() == num_dims);
}

TEST(SpsPlaneTest, GridParametersStayWithinUnitCube)
{
    const Eigen::Vector2d x_base(0.25, 0.25);
    const Eigen::Vector2d x_opposite(0.75, 0.75);
    const Eigen::Vector2d cross_vec(0.0, 0.25);

    const auto plane = sps::SimplePlane(x_base, x_opposite, cross_vec);

    constexpr unsigned num_candidates      = 5;
    constexpr double   inter_level_scale   = 0.5;
    constexpr int      grid_radius         = (num_candidates - 1) / 2;
    const auto         selected_grid_cell  = sps::AbstractPlane::GridCellIndex(grid_radius, 0);

    std::vector<sps::AbstractPlane::GridCellIndex> prev_grid_cells;

    for (int depth = 0; depth < 5; ++depth)
    {
        const auto parameters =
            plane.CalcGridParameters(selected_grid_cell, num_candidates, inter_level_scale, prev_grid_cells);

        EXPECT_GE(parameters.minCoeff(), 0.0);
        EXPECT_LE(parameters.maxCoeff(), 1.0);

        prev_grid_cells.push_back(selected_grid_cell);
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
