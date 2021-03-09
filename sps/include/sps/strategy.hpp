#ifndef SPS_STRATEGY_HPP
#define SPS_STRATEGY_HPP

#include <Eigen/Core>
#include <memory>
#include <sps/plane.hpp>

namespace sequential_line_search
{
    class PreferenceRegressor;
}

namespace sps
{
    enum class PlaneCenterSelectionStrategy
    {
        LargestExpectValue,
        LastSelection
    };

    // -----------------------------------------------------------------------------------------------------------------

    class PlaneInitStrategy
    {
    public:
        virtual std::shared_ptr<AbstractPlane> Init(const unsigned num_dims) const = 0;
    };

    class FixedCenterFreeCrossVecsInitStrategy final : public PlaneInitStrategy
    {
    public:
        std::shared_ptr<AbstractPlane> Init(const unsigned num_dims) const override final;
    };

    // -----------------------------------------------------------------------------------------------------------------

    class PlaneConstructionStrategy
    {
    public:
        /// \brief Construct a plane based on a strategy implemented in the specific class.
        ///
        /// \param x_current_best The point that is expected to be the best one among the set of so-far observed points.
        /// This will be the center of the constructed plane. In the original SLS paper [SIGGRAPH 2017], this is set as
        /// the one that has the best latent value; however, it may be more reasonable to use the one the user selected
        /// in the last subtask.
        virtual std::shared_ptr<AbstractPlane>
        Construct(const std::shared_ptr<sequential_line_search::PreferenceRegressor> regressor,
                  const Eigen::VectorXd&                                             x_current_best) const = 0;
    };

    class JointCrossVecEstimationStrategy final : public PlaneConstructionStrategy
    {
    public:
        std::shared_ptr<AbstractPlane>
        Construct(const std::shared_ptr<sequential_line_search::PreferenceRegressor> regressor,
                  const Eigen::VectorXd& x_current_best) const override final;
    };

    class FixedCenterFreeCrossVecsStrategy final : public PlaneConstructionStrategy
    {
    public:
        std::shared_ptr<AbstractPlane>
        Construct(const std::shared_ptr<sequential_line_search::PreferenceRegressor> regressor,
                  const Eigen::VectorXd& x_current_best) const override final;
    };

    class BatchEiPiecewisePlanarStrategy final : public PlaneConstructionStrategy
    {
    public:
        std::shared_ptr<AbstractPlane>
        Construct(const std::shared_ptr<sequential_line_search::PreferenceRegressor> regressor,
                  const Eigen::VectorXd& x_current_best) const override final;
    };

    class ApproxIntegralJointCrossVecEstimationStrategy final : public PlaneConstructionStrategy
    {
    public:
        std::shared_ptr<AbstractPlane>
        Construct(const std::shared_ptr<sequential_line_search::PreferenceRegressor> regressor,
                  const Eigen::VectorXd& x_current_best) const override final;
    };

    class FirstEiThenPlaneIntegralStrategy final : public PlaneConstructionStrategy
    {
    public:
        std::shared_ptr<AbstractPlane>
        Construct(const std::shared_ptr<sequential_line_search::PreferenceRegressor> regressor,
                  const Eigen::VectorXd& x_current_best) const override final;
    };
} // namespace sps

#endif // SPS_STRATEGY_HPP
