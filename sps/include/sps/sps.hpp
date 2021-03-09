#ifndef SPS_SPS_HPP
#define SPS_SPS_HPP

#include <Eigen/Core>
#include <memory>
#include <sps/plane.hpp>
#include <sps/strategy.hpp>
#include <vector>

namespace sequential_line_search
{
    class PreferenceRegressor;
    class PreferenceDataManager;
} // namespace sequential_line_search

namespace sps
{
    class Optimizer
    {
    public:
        Optimizer(
            const unsigned                     num_dims,
            const bool                         use_map_hyperparams = true,
            std::shared_ptr<PlaneInitStrategy> init_strategy = std::make_shared<FixedCenterFreeCrossVecsInitStrategy>(),
            std::shared_ptr<PlaneConstructionStrategy> construction_strategy =
                std::make_shared<FirstEiThenPlaneIntegralStrategy>(),
            const PlaneCenterSelectionStrategy center_strategy = PlaneCenterSelectionStrategy::LastSelection);

        void SubmitData(const Eigen::VectorXd& x_preferred, const std::vector<Eigen::VectorXd>& x_others);
        std::shared_ptr<const AbstractPlane> RetrieveSearchPlane() const;

        double GetCurrentPredictiveMean(const Eigen::VectorXd& params) const;
        double GetCurrentPredictiveVar(const Eigen::VectorXd& params) const;
        double GetCurrentExpectedImprovement(const Eigen::VectorXd& params) const;

    private:
        const bool                                       m_use_map_hyperparams;
        const std::shared_ptr<PlaneInitStrategy>         m_init_strategy;
        const std::shared_ptr<PlaneConstructionStrategy> m_construction_strategy;
        const PlaneCenterSelectionStrategy               m_center_strategy;

        const unsigned m_num_dims;

        std::shared_ptr<AbstractPlane> m_plane;

        std::shared_ptr<sequential_line_search::PreferenceRegressor>   m_regressor;
        std::shared_ptr<sequential_line_search::PreferenceDataManager> m_data;
    };
} // namespace sps

#endif // SPS_SPS_HPP
