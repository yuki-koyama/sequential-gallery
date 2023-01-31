#include <iostream>
#include <mathtoolbox/acquisition-functions.hpp>
#include <sequential-line-search/preference-data-manager.hpp>
#include <sequential-line-search/preference-regressor.hpp>
#include <sps/sps.hpp>

// ---------------------------------------------------------------------------------------------------------------------

constexpr bool   k_is_verbose           = false;
constexpr double k_default_signal_level = 0.200;
constexpr double k_default_length_scale = 0.500;
constexpr double k_prior_variance       = 0.010;
constexpr double k_btl_scale            = 0.010;
constexpr double k_merge_threshold      = 1e-03;

// ---------------------------------------------------------------------------------------------------------------------

sps::Optimizer::Optimizer(const unsigned                             num_dims,
                          const bool                                 use_map_hyperparams,
                          std::shared_ptr<PlaneInitStrategy>         init_strategy,
                          std::shared_ptr<PlaneConstructionStrategy> construction_strategy,
                          const PlaneCenterSelectionStrategy         center_strategy)
    : m_use_map_hyperparams(use_map_hyperparams), m_init_strategy(init_strategy),
      m_construction_strategy(construction_strategy), m_center_strategy(center_strategy), m_num_dims(num_dims)
{
    m_plane     = m_init_strategy->Init(m_num_dims);
    m_data      = std::make_shared<sequential_line_search::PreferenceDataManager>();
    m_regressor = nullptr;
}

void sps::Optimizer::SubmitData(const Eigen::VectorXd& x_preferred, const std::vector<Eigen::VectorXd>& x_others)
{
    m_data->AddNewPoints(x_preferred, x_others, true, k_merge_threshold);

    m_regressor = std::make_shared<sequential_line_search::PreferenceRegressor>(m_data->GetX(),
                                                                                m_data->GetD(),
                                                                                m_use_map_hyperparams,
                                                                                k_default_signal_level,
                                                                                k_default_length_scale,
                                                                                0.0,
                                                                                k_prior_variance,
                                                                                k_btl_scale);

    if constexpr (k_is_verbose)
    {
        std::cout << "Estimated point values: " << std::endl;
        if (m_regressor->GetSmallY().size() > 8)
        {
            std::cout << m_regressor->GetSmallY().segment<8>(0).transpose().format(Eigen::IOFormat(3)) << " ..."
                      << std::endl;
        }
        else
        {
            std::cout << m_regressor->GetSmallY().transpose().format(Eigen::IOFormat(2)) << std::endl;
        }

        if (m_use_map_hyperparams)
        {
            std::cout << "Learned hyperparameters: " << std::endl;
            std::cout << m_regressor->GetKernelHyperparams().transpose().format(Eigen::IOFormat(3)) << std::endl;
        }
    }

    const Eigen::VectorXd x_current_best = [&]() {
        switch (m_center_strategy)
        {
        case PlaneCenterSelectionStrategy::LargestExpectValue:
            return m_regressor->FindArgMax();
        case PlaneCenterSelectionStrategy::LastSelection:
            return x_preferred;
        }
    }();

    m_plane = m_construction_strategy->Construct(m_regressor, x_current_best);
}

std::shared_ptr<const sps::AbstractPlane> sps::Optimizer::RetrieveSearchPlane() const { return m_plane; }

double sps::Optimizer::GetCurrentPredictiveMean(const Eigen::VectorXd& params) const
{
    return (m_regressor.get() == nullptr) ? 0.0 : m_regressor->PredictMu(params);
}

double sps::Optimizer::GetCurrentPredictiveVar(const Eigen::VectorXd& params) const
{
    return (m_regressor.get() == nullptr) ? 0.0 : m_regressor->PredictSigma(params);
}

double sps::Optimizer::GetCurrentExpectedImprovement(const Eigen::VectorXd& params) const
{
    if (m_regressor.get() == nullptr)
    {
        return 0.0;
    }

    const auto ac_func = [&](const Eigen::VectorXd& x) {
        const auto mu    = [&](const Eigen::VectorXd& x) { return m_regressor->PredictMu(x); };
        const auto sigma = [&](const Eigen::VectorXd& x) { return m_regressor->PredictSigma(x); };

        return mathtoolbox::GetExpectedImprovement(x, mu, sigma, m_regressor->FindArgMax());
    };

    return ac_func(params);
}
