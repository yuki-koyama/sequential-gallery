#include <gtest/gtest.h>
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

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
