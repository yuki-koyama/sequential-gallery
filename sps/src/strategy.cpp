#include <numeric>
#include <cassert>
#include <mathtoolbox/acquisition-functions.hpp>
#include <nlopt-util.hpp>
#include <parallel-util.hpp>
#include <sequential-line-search/acquisition-function.hpp>
#include <sequential-line-search/preference-regressor.hpp>
#include <sps/strategy.hpp>

using Eigen::VectorXd;

// ---------------------------------------------------------------------------------------------------------------------

constexpr bool k_verbose_nlopt                = false;
constexpr bool k_ensure_opposite_vertex_bound = true;

// ---------------------------------------------------------------------------------------------------------------------

namespace sps::internal
{
    template <class T, class BinaryOp>
    VectorXd TransformAndReduceVecs(const std::vector<VectorXd>& vec_inputs,
                                    const std::vector<T>&        other_inputs,
                                    BinaryOp                     binary_op)
    {
#if false
        // Sadly, std::transform_reduce is a C++17 feature that has not been shipped with GCC yet.
        return std::transform_reduce(
            vec_inputs.begin(),
            vec_inputs.end(),
            other_inputs.begin(),
            VectorXd::Zero(vec_inputs.front().size()).eval(),
            [](const VectorXd& acc, const VectorXd& vec) { return acc + vec; },
            binary_op);
#else
        // This is based on C++03 features, which is probably slower than the above.
        std::vector<VectorXd> transformed_values(vec_inputs.size());
        std::transform(
            vec_inputs.begin(), vec_inputs.end(), other_inputs.begin(), transformed_values.begin(), binary_op);
        return std::accumulate(
            transformed_values.begin(), transformed_values.end(), VectorXd::Zero(vec_inputs.front().size()).eval());
#endif
    }

    template <class UnaryOp> double TransformAndReduceVecs(const std::vector<VectorXd>& vec_inputs, UnaryOp unary_op)
    {
#if false
        // Sadly, std::transform_reduce is a C++17 feature that has not been shipped with GCC yet.
        return std::transform_reduce(
            vec_inputs.begin(),
            vec_inputs.end(),
            0.0,
            [](const double acc, const double scalar) { return acc + scalar; },
            unary_op);
#else
        // This is based on C++03 features, which is probably slower than the above.
        std::vector<double> transformed_values(vec_inputs.size());
        std::transform(vec_inputs.begin(), vec_inputs.end(), transformed_values.begin(), unary_op);
        return std::accumulate(transformed_values.begin(), transformed_values.end(), 0.0);
#endif
    }

    VectorXd GetRandomSample(const VectorXd& lower_bound, const VectorXd& upper_bound);

    std::pair<VectorXd, VectorXd> CreateCrossVecUpperAndLowerBoundsForSymmetricPlane(const VectorXd& center);

    double   CalcOrthogonalityCost(const VectorXd& cross_vec_u, const VectorXd& cross_vec_v);
    VectorXd CalcOrthogonalityCostVecUDerivative(const VectorXd& cross_vec_u, const VectorXd& cross_vec_v);
    VectorXd CalcOrthogonalityCostVecVDerivative(const VectorXd& cross_vec_u, const VectorXd& cross_vec_v);

    double   CalcLengthCost(const VectorXd& cross_vec, const double target_length);
    VectorXd CalcLengthCostDerivative(const VectorXd& cross_vec, const double target_length);

    inline bool CheckIfParamIsInsideBound(const VectorXd& x) { return x.minCoeff() >= 0.0 && x.maxCoeff() <= 1.0; }
} // namespace sps::internal

// ---------------------------------------------------------------------------------------------------------------------

VectorXd sps::internal::GetRandomSample(const VectorXd& lower_bound, const VectorXd& upper_bound)
{
    const auto num_dims      = lower_bound.size();
    const auto random_sample = 0.5 * (Eigen::ArrayXd::Random(num_dims) + Eigen::ArrayXd::Ones(num_dims));

    return lower_bound + ((upper_bound - lower_bound).array() * random_sample).matrix();
}

std::pair<VectorXd, VectorXd> sps::internal::CreateCrossVecUpperAndLowerBoundsForSymmetricPlane(const VectorXd& center)
{
    const unsigned num_dims = center.size();

    VectorXd upper_bound(num_dims * 2);
    VectorXd lower_bound(num_dims * 2);

    for (unsigned i = 0; i < num_dims; ++i)
    {
        // The following lines assume x in [0, 1]^n
        const double dist_to_closest_bound = std::min(std::abs(center(i)), std::abs(1.0 - center(i)));

        upper_bound(i + num_dims * 0) = dist_to_closest_bound;
        upper_bound(i + num_dims * 1) = dist_to_closest_bound;
        lower_bound(i + num_dims * 0) = -dist_to_closest_bound;
        lower_bound(i + num_dims * 1) = -dist_to_closest_bound;
    }

    return {upper_bound, lower_bound};
}

double sps::internal::CalcOrthogonalityCost(const VectorXd& cross_vec_u, const VectorXd& cross_vec_v)
{
    const double dot_prod = cross_vec_u.dot(cross_vec_v);

    return dot_prod * dot_prod;
}

VectorXd sps::internal::CalcOrthogonalityCostVecUDerivative(const VectorXd& cross_vec_u, const VectorXd& cross_vec_v)
{
    return 2.0 * cross_vec_v * cross_vec_v.transpose() * cross_vec_u;
}

VectorXd sps::internal::CalcOrthogonalityCostVecVDerivative(const VectorXd& cross_vec_u, const VectorXd& cross_vec_v)
{
    return 2.0 * cross_vec_u * cross_vec_u.transpose() * cross_vec_v;
}

double sps::internal::CalcLengthCost(const VectorXd& cross_vec, const double target_length)
{
    return (cross_vec.squaredNorm() - target_length * target_length) *
           (cross_vec.squaredNorm() - target_length * target_length);
}

VectorXd sps::internal::CalcLengthCostDerivative(const VectorXd& cross_vec, const double target_length)
{
    return 4.0 * (cross_vec.squaredNorm() - target_length * target_length) * cross_vec;
}

// ---------------------------------------------------------------------------------------------------------------------

std::shared_ptr<sps::AbstractPlane> sps::FixedCenterFreeCrossVecsInitStrategy::Init(const unsigned num_dims) const
{
    const VectorXd upper_bound = VectorXd::Constant(2 * num_dims, +0.5);
    const VectorXd lower_bound = VectorXd::Constant(2 * num_dims, -0.5);

    constexpr double target_length = 1.0;

    const auto objective_func = [&](const VectorXd& x) {
        const VectorXd cross_vec_u = x.segment(0 * num_dims, num_dims);
        const VectorXd cross_vec_v = x.segment(1 * num_dims, num_dims);

        return internal::CalcLengthCost(cross_vec_u, target_length) +
               internal::CalcLengthCost(cross_vec_v, target_length) +
               internal::CalcOrthogonalityCost(cross_vec_u, cross_vec_v);
    };

    const auto grad_func = [&](const VectorXd& x) {
        const VectorXd cross_vec_u = x.segment(0 * num_dims, num_dims);
        const VectorXd cross_vec_v = x.segment(1 * num_dims, num_dims);

        VectorXd grad(num_dims * 2);

        grad.segment(0 * num_dims, num_dims) = internal::CalcLengthCostDerivative(cross_vec_u, target_length) +
                                               internal::CalcOrthogonalityCostVecUDerivative(cross_vec_u, cross_vec_v);
        grad.segment(1 * num_dims, num_dims) = internal::CalcLengthCostDerivative(cross_vec_v, target_length) +
                                               internal::CalcOrthogonalityCostVecVDerivative(cross_vec_u, cross_vec_v);

        return grad;
    };

    const VectorXd initial_cross_vecs = internal::GetRandomSample(lower_bound, upper_bound);

    const VectorXd cross_vecs = nloptutil::unconstrained::grad_based::bounded::solve(initial_cross_vecs,
                                                                                     upper_bound,
                                                                                     lower_bound,
                                                                                     objective_func,
                                                                                     grad_func,
                                                                                     nlopt::LD_LBFGS,
                                                                                     false,
                                                                                     1000,
                                                                                     1e-06,
                                                                                     1e-06,
                                                                                     k_verbose_nlopt);

    const auto cross_vec_u = cross_vecs.segment(0 * num_dims, num_dims);
    const auto cross_vec_v = cross_vecs.segment(1 * num_dims, num_dims);

    const auto center = VectorXd::Constant(num_dims, 0.50);

    return std::make_shared<SimplePlane>(center - cross_vec_u, center + cross_vec_u, cross_vec_v);
}

// ---------------------------------------------------------------------------------------------------------------------

std::shared_ptr<sps::AbstractPlane> sps::JointCrossVecEstimationStrategy::Construct(
    const std::shared_ptr<sequential_line_search::PreferenceRegressor> regressor,
    const Eigen::VectorXd&                                             x_current_best) const
{
    const unsigned num_dims = regressor->GetNumDims();

    // This is different from x_current_best; x_plus is used for EI calculation only.
    const VectorXd x_plus = regressor->FindArgMax();

    const auto  bounds      = internal::CreateCrossVecUpperAndLowerBoundsForSymmetricPlane(x_current_best);
    const auto& upper_bound = std::get<0>(bounds);
    const auto& lower_bound = std::get<1>(bounds);

    const auto ac_func = [&](const VectorXd& x) {
        const auto mu    = [&](const VectorXd& x) { return regressor->PredictMu(x); };
        const auto sigma = [&](const VectorXd& x) { return regressor->PredictSigma(x); };

        return mathtoolbox::GetExpectedImprovement(x, mu, sigma, x_plus);
    };

    const auto ac_grad_func = [&](const VectorXd& x) {
        const auto mu               = [&](const VectorXd& x) { return regressor->PredictMu(x); };
        const auto sigma            = [&](const VectorXd& x) { return regressor->PredictSigma(x); };
        const auto mu_derivative    = [&](const VectorXd& x) { return regressor->PredictMuDerivative(x); };
        const auto sigma_derivative = [&](const VectorXd& x) { return regressor->PredictSigmaDerivative(x); };

        return mathtoolbox::GetExpectedImprovementDerivative(x, mu, sigma, x_plus, mu_derivative, sigma_derivative);
    };

    constexpr double alpha = 1.0;

    const auto objective_func = [&](const VectorXd& x) {
        const VectorXd cross_vec_u = x.segment(0 * num_dims, num_dims);
        const VectorXd cross_vec_v = x.segment(1 * num_dims, num_dims);

        const VectorXd y_p_u = x_current_best + cross_vec_u;
        const VectorXd y_m_u = x_current_best - cross_vec_u;
        const VectorXd y_p_v = x_current_best + cross_vec_v;
        const VectorXd y_m_v = x_current_best - cross_vec_v;

        const double orthogonality_cost = internal::CalcOrthogonalityCost(cross_vec_u, cross_vec_v);

        return ac_func(y_p_u) + ac_func(y_m_u) + ac_func(y_p_v) + ac_func(y_m_v) - alpha * orthogonality_cost;
    };

    const auto grad_func = [&](const VectorXd& x) {
        const VectorXd cross_vec_u = x.segment(0 * num_dims, num_dims);
        const VectorXd cross_vec_v = x.segment(1 * num_dims, num_dims);

        const VectorXd y_p_u = x_plus + cross_vec_u;
        const VectorXd y_m_u = x_plus - cross_vec_u;
        const VectorXd y_p_v = x_plus + cross_vec_v;
        const VectorXd y_m_v = x_plus - cross_vec_v;

        VectorXd grad(num_dims * 2);

        grad.segment(0 * num_dims, num_dims) =
            ac_grad_func(y_p_u) - ac_grad_func(y_m_u) -
            alpha * internal::CalcOrthogonalityCostVecUDerivative(cross_vec_u, cross_vec_v);
        grad.segment(1 * num_dims, num_dims) =
            ac_grad_func(y_p_v) - ac_grad_func(y_m_v) -
            alpha * internal::CalcOrthogonalityCostVecVDerivative(cross_vec_u, cross_vec_v);

        return grad;
    };

    constexpr unsigned num_trials = 10;

    std::vector<VectorXd> solutions(num_trials);
    std::vector<double>   objective_values(num_trials);

    const auto process_performer = [&](const unsigned trial) {
        const VectorXd initial_joint_cross_vecs = internal::GetRandomSample(lower_bound, upper_bound);

        assert(initial_joint_cross_vecs.array().min(lower_bound.array()).isApprox(lower_bound.array()));
        assert(initial_joint_cross_vecs.array().max(upper_bound.array()).isApprox(upper_bound.array()));

        const VectorXd joint_cross_vecs = nloptutil::unconstrained::grad_based::bounded::solve(initial_joint_cross_vecs,
                                                                                               upper_bound,
                                                                                               lower_bound,
                                                                                               objective_func,
                                                                                               grad_func,
                                                                                               nlopt::LD_LBFGS,
                                                                                               true,
                                                                                               1000,
                                                                                               1e-06,
                                                                                               1e-06,
                                                                                               k_verbose_nlopt);

        solutions[trial]        = joint_cross_vecs;
        objective_values[trial] = objective_func(joint_cross_vecs);
    };

    parallelutil::parallel_for(num_trials, process_performer);

    const auto max_index =
        std::distance(objective_values.begin(), std::max_element(objective_values.begin(), objective_values.end()));
    const auto& joint_cross_vecs = solutions[max_index];

    const auto cross_vec_u = joint_cross_vecs.segment(0 * num_dims, num_dims);
    const auto cross_vec_v = joint_cross_vecs.segment(1 * num_dims, num_dims);

    return std::make_shared<SimplePlane>(x_plus - cross_vec_u, x_plus + cross_vec_u, cross_vec_v);
}

std::shared_ptr<sps::AbstractPlane> sps::FixedCenterFreeCrossVecsStrategy::Construct(
    const std::shared_ptr<sequential_line_search::PreferenceRegressor> regressor,
    const Eigen::VectorXd&                                             x_current_best) const
{
    const unsigned num_dims = regressor->GetNumDims();
    const auto [upper_bound, lower_bound] =
        internal::CreateCrossVecUpperAndLowerBoundsForSymmetricPlane(x_current_best);

    const VectorXd initial_joint_cross_vecs = internal::GetRandomSample(lower_bound, upper_bound);

    assert(initial_joint_cross_vecs.array().min(lower_bound.array()).isApprox(lower_bound.array()));
    assert(initial_joint_cross_vecs.array().max(upper_bound.array()).isApprox(upper_bound.array()));

    constexpr double target_length = 1.0;

    const auto objective_func = [&](const VectorXd& x) {
        const VectorXd cross_vec_u = x.segment(0 * num_dims, num_dims);
        const VectorXd cross_vec_v = x.segment(1 * num_dims, num_dims);

        return internal::CalcLengthCost(cross_vec_u, target_length) +
               internal::CalcLengthCost(cross_vec_v, target_length) +
               internal::CalcOrthogonalityCost(cross_vec_u, cross_vec_v);
    };

    const auto grad_func = [&](const VectorXd& x) {
        const VectorXd cross_vec_u = x.segment(0 * num_dims, num_dims);
        const VectorXd cross_vec_v = x.segment(1 * num_dims, num_dims);

        VectorXd grad(num_dims * 2);

        grad.segment(0 * num_dims, num_dims) = internal::CalcLengthCostDerivative(cross_vec_u, target_length) +
                                               internal::CalcOrthogonalityCostVecUDerivative(cross_vec_u, cross_vec_v);
        grad.segment(1 * num_dims, num_dims) = internal::CalcLengthCostDerivative(cross_vec_v, target_length) +
                                               internal::CalcOrthogonalityCostVecVDerivative(cross_vec_u, cross_vec_v);

        return grad;
    };

    const VectorXd initial_cross_vecs = internal::GetRandomSample(lower_bound, upper_bound);

    const VectorXd cross_vecs = nloptutil::unconstrained::grad_based::bounded::solve(initial_cross_vecs,
                                                                                     upper_bound,
                                                                                     lower_bound,
                                                                                     objective_func,
                                                                                     grad_func,
                                                                                     nlopt::LD_LBFGS,
                                                                                     false,
                                                                                     1000,
                                                                                     1e-06,
                                                                                     1e-06,
                                                                                     k_verbose_nlopt);

    const auto cross_vec_u = cross_vecs.segment(0 * num_dims, num_dims);
    const auto cross_vec_v = cross_vecs.segment(1 * num_dims, num_dims);

    return std::make_shared<SimplePlane>(x_current_best - cross_vec_u, x_current_best + cross_vec_u, cross_vec_v);
}

std::shared_ptr<sps::AbstractPlane> sps::BatchEiPiecewisePlanarStrategy::Construct(
    const std::shared_ptr<sequential_line_search::PreferenceRegressor> regressor,
    const Eigen::VectorXd&                                             x_current_best) const
{
    // TODO: Check whether the number of iterations is appropriate, or not
    const auto xs_ei = sequential_line_search::acquisition_func::FindNextPoints(*regressor, 4);

    const auto planarity_cost_evaluator =
        [](const VectorXd& x_0, const VectorXd& x_1, const VectorXd& x_2, const VectorXd& x_3) {
            const double cos_sim_02 = x_0.normalized().dot(x_2.normalized());
            const double cos_sim_13 = x_1.normalized().dot(x_3.normalized());

            return cos_sim_02 + cos_sim_13;
        };

    const double order_0_cost = planarity_cost_evaluator(xs_ei[0], xs_ei[1], xs_ei[2], xs_ei[3]);
    const double order_1_cost = planarity_cost_evaluator(xs_ei[1], xs_ei[0], xs_ei[2], xs_ei[3]);

    if (order_0_cost < order_1_cost)
    {
        return std::make_shared<PiecewisePlanarPlane>(x_current_best, xs_ei[0], xs_ei[1], xs_ei[2], xs_ei[3]);
    }
    else
    {
        return std::make_shared<PiecewisePlanarPlane>(x_current_best, xs_ei[1], xs_ei[0], xs_ei[2], xs_ei[3]);
    }
}

std::shared_ptr<sps::AbstractPlane> sps::ApproxIntegralJointCrossVecEstimationStrategy::Construct(
    const std::shared_ptr<sequential_line_search::PreferenceRegressor> regressor,
    const Eigen::VectorXd&                                             x_current_best) const
{
    const unsigned num_dims = regressor->GetNumDims();

    // This is different from x_current_best; x_plus is used for EI calculation only.
    const VectorXd x_plus = regressor->FindArgMax();

    const auto  bounds      = internal::CreateCrossVecUpperAndLowerBoundsForSymmetricPlane(x_current_best);
    const auto& upper_bound = std::get<0>(bounds);
    const auto& lower_bound = std::get<1>(bounds);

    const auto ac_func = [&](const VectorXd& x) {
        const auto mu    = [&](const VectorXd& x) { return regressor->PredictMu(x); };
        const auto sigma = [&](const VectorXd& x) { return regressor->PredictSigma(x); };

        return mathtoolbox::GetExpectedImprovement(x, mu, sigma, x_plus);
    };

    const auto ac_grad_func = [&](const VectorXd& x) {
        const auto mu               = [&](const VectorXd& x) { return regressor->PredictMu(x); };
        const auto sigma            = [&](const VectorXd& x) { return regressor->PredictSigma(x); };
        const auto mu_derivative    = [&](const VectorXd& x) { return regressor->PredictMuDerivative(x); };
        const auto sigma_derivative = [&](const VectorXd& x) { return regressor->PredictSigmaDerivative(x); };

        return mathtoolbox::GetExpectedImprovementDerivative(x, mu, sigma, x_plus, mu_derivative, sigma_derivative);
    };

    constexpr double alpha = 1.0;

    constexpr unsigned sample_resolution = 5;
    constexpr double   weight_per_sample = 1.0 / static_cast<double>(sample_resolution * sample_resolution);

    std::vector<double> ws_u(sample_resolution * sample_resolution);
    std::vector<double> ws_v(sample_resolution * sample_resolution);
    for (unsigned i = 0; i < sample_resolution; ++i)
    {
        const double w_i = 2.0 * static_cast<double>(i) / static_cast<double>(sample_resolution - 1) - 1.0;

        for (unsigned j = 0; j < sample_resolution; ++j)
        {
            const double w_j = 2.0 * static_cast<double>(j) / static_cast<double>(sample_resolution - 1) - 1.0;

            ws_u[i * sample_resolution + j] = w_i;
            ws_v[i * sample_resolution + j] = w_j;
        }
    }

    const auto objective_func = [&](const VectorXd& x) {
        const VectorXd cross_vec_u = x.segment(0 * num_dims, num_dims);
        const VectorXd cross_vec_v = x.segment(1 * num_dims, num_dims);

        std::vector<VectorXd> sample_points(sample_resolution * sample_resolution);
        for (unsigned i = 0; i < sample_resolution; ++i)
        {
            for (unsigned j = 0; j < sample_resolution; ++j)
            {
                sample_points[i * sample_resolution + j] = x_current_best +
                                                           ws_u[i * sample_resolution + j] * cross_vec_u +
                                                           ws_v[i * sample_resolution + j] * cross_vec_v;
            }
        }

        const double ac_sum =
            internal::TransformAndReduceVecs(sample_points, [&](const VectorXd& x) { return ac_func(x); });

        const double orthogonality_cost = internal::CalcOrthogonalityCost(cross_vec_u, cross_vec_v);

        return weight_per_sample * ac_sum - alpha * orthogonality_cost;
    };

    const auto grad_func = [&](const VectorXd& x) {
        const VectorXd cross_vec_u = x.segment(0 * num_dims, num_dims);
        const VectorXd cross_vec_v = x.segment(1 * num_dims, num_dims);

        std::vector<VectorXd> sample_points(sample_resolution * sample_resolution);
        for (unsigned i = 0; i < sample_resolution; ++i)
        {
            for (unsigned j = 0; j < sample_resolution; ++j)
            {
                sample_points[i * sample_resolution + j] = x_current_best +
                                                           ws_u[i * sample_resolution + j] * cross_vec_u +
                                                           ws_v[i * sample_resolution + j] * cross_vec_v;
            }
        }

        const VectorXd ac_sum_vec_u_derivative = sps::internal::TransformAndReduceVecs(
            sample_points, ws_u, [&](const VectorXd& sample_point, const double w_u) -> VectorXd {
                return w_u * ac_grad_func(sample_point);
            });
        const VectorXd ac_sum_vec_v_derivative = sps::internal::TransformAndReduceVecs(
            sample_points, ws_v, [&](const VectorXd& sample_point, const double w_v) -> VectorXd {
                return w_v * ac_grad_func(sample_point);
            });

        VectorXd grad(num_dims * 2);

        grad.segment(0 * num_dims, num_dims) =
            weight_per_sample * ac_sum_vec_u_derivative -
            alpha * internal::CalcOrthogonalityCostVecUDerivative(cross_vec_u, cross_vec_v);
        grad.segment(1 * num_dims, num_dims) =
            weight_per_sample * ac_sum_vec_v_derivative -
            alpha * internal::CalcOrthogonalityCostVecVDerivative(cross_vec_u, cross_vec_v);

        return grad;
    };

    constexpr unsigned num_trials = 10;

    std::vector<VectorXd> solutions(num_trials);
    std::vector<double>   objective_values(num_trials);

    auto process_performer = [&](const unsigned trial) {
        const VectorXd initial_joint_cross_vecs = internal::GetRandomSample(lower_bound, upper_bound);

        assert(initial_joint_cross_vecs.array().min(lower_bound.array()).isApprox(lower_bound.array()));
        assert(initial_joint_cross_vecs.array().max(upper_bound.array()).isApprox(upper_bound.array()));

        const VectorXd joint_cross_vecs = nloptutil::unconstrained::grad_based::bounded::solve(initial_joint_cross_vecs,
                                                                                               upper_bound,
                                                                                               lower_bound,
                                                                                               objective_func,
                                                                                               grad_func,
                                                                                               nlopt::LD_LBFGS,
                                                                                               true,
                                                                                               1000,
                                                                                               1e-06,
                                                                                               1e-06,
                                                                                               k_verbose_nlopt);

        solutions[trial]        = joint_cross_vecs;
        objective_values[trial] = objective_func(joint_cross_vecs);
    };

    parallelutil::parallel_for(num_trials, process_performer);

    const auto max_index =
        std::distance(objective_values.begin(), std::max_element(objective_values.begin(), objective_values.end()));
    const auto& joint_cross_vecs = solutions[max_index];

    const auto cross_vec_u = joint_cross_vecs.segment(0 * num_dims, num_dims);
    const auto cross_vec_v = joint_cross_vecs.segment(1 * num_dims, num_dims);

    return std::make_shared<SimplePlane>(x_current_best - cross_vec_u, x_current_best + cross_vec_u, cross_vec_v);
}

std::shared_ptr<sps::AbstractPlane> sps::FirstEiThenPlaneIntegralStrategy::Construct(
    const std::shared_ptr<sequential_line_search::PreferenceRegressor> regressor,
    const Eigen::VectorXd&                                             x_current_best) const
{
    constexpr unsigned num_trials = 10;

    const unsigned num_dims = regressor->GetNumDims();

    // This is different from x_current_best; x_plus is used for EI calculation only.
    const VectorXd x_plus = regressor->FindArgMax();

    // TODO: Check whether the number of iterations is appropriate, or not
    const VectorXd x_ei = sequential_line_search::acquisition_func::FindNextPoint(*regressor, num_trials);

    const VectorXd cross_vec_u = x_ei - x_current_best;
    const auto     bounds      = internal::CreateCrossVecUpperAndLowerBoundsForSymmetricPlane(x_current_best);
    const VectorXd upper_bound = std::get<0>(bounds).segment(0, num_dims);
    const VectorXd lower_bound = std::get<1>(bounds).segment(0, num_dims);

    const auto ac_func = [&](const VectorXd& x) {
        const auto mu    = [&](const VectorXd& x) { return regressor->PredictMu(x); };
        const auto sigma = [&](const VectorXd& x) { return regressor->PredictSigma(x); };

        return mathtoolbox::GetExpectedImprovement(x, mu, sigma, x_plus);
    };

    const auto ac_grad_func = [&](const VectorXd& x) {
        const auto mu               = [&](const VectorXd& x) { return regressor->PredictMu(x); };
        const auto sigma            = [&](const VectorXd& x) { return regressor->PredictSigma(x); };
        const auto mu_derivative    = [&](const VectorXd& x) { return regressor->PredictMuDerivative(x); };
        const auto sigma_derivative = [&](const VectorXd& x) { return regressor->PredictSigmaDerivative(x); };

        return mathtoolbox::GetExpectedImprovementDerivative(x, mu, sigma, x_plus, mu_derivative, sigma_derivative);
    };

    constexpr double alpha = 1.0;

    constexpr unsigned sample_resolution = 5;
    constexpr double   weight_per_sample = 1.0 / static_cast<double>(sample_resolution * sample_resolution);

    std::vector<double> ws_u(sample_resolution * sample_resolution);
    std::vector<double> ws_v(sample_resolution * sample_resolution);
    for (unsigned i = 0; i < sample_resolution; ++i)
    {
        const double w_i = 2.0 * static_cast<double>(i) / static_cast<double>(sample_resolution - 1) - 1.0;

        for (unsigned j = 0; j < sample_resolution; ++j)
        {
            const double w_j = 2.0 * static_cast<double>(j) / static_cast<double>(sample_resolution - 1) - 1.0;

            ws_u[i * sample_resolution + j] = w_i;
            ws_v[i * sample_resolution + j] = w_j;
        }
    }

    const auto objective_func = [&](const VectorXd& x) {
        const VectorXd& cross_vec_v = x;

        std::vector<VectorXd> sample_points(sample_resolution * sample_resolution);
        for (unsigned i = 0; i < sample_resolution; ++i)
        {
            for (unsigned j = 0; j < sample_resolution; ++j)
            {
                sample_points[i * sample_resolution + j] = x_current_best +
                                                           ws_u[i * sample_resolution + j] * cross_vec_u +
                                                           ws_v[i * sample_resolution + j] * cross_vec_v;
            }
        }

        const double ac_sum =
            internal::TransformAndReduceVecs(sample_points, [&](const VectorXd& x) { return ac_func(x); });

        const double orthogonality_cost = internal::CalcOrthogonalityCost(cross_vec_u, cross_vec_v);

        return weight_per_sample * ac_sum - alpha * orthogonality_cost;
    };

    const auto grad_func = [&](const VectorXd& x) {
        const VectorXd& cross_vec_v = x;

        std::vector<VectorXd> sample_points(sample_resolution * sample_resolution);
        for (unsigned i = 0; i < sample_resolution; ++i)
        {
            for (unsigned j = 0; j < sample_resolution; ++j)
            {
                sample_points[i * sample_resolution + j] = x_current_best +
                                                           ws_u[i * sample_resolution + j] * cross_vec_u +
                                                           ws_v[i * sample_resolution + j] * cross_vec_v;
            }
        }

        const VectorXd ac_sum_vec_v_derivative = sps::internal::TransformAndReduceVecs(
            sample_points, ws_v, [&](const VectorXd& sample_point, const double w_v) -> VectorXd {
                return w_v * ac_grad_func(sample_point);
            });

        const VectorXd grad = weight_per_sample * ac_sum_vec_v_derivative -
                              alpha * internal::CalcOrthogonalityCostVecVDerivative(cross_vec_u, cross_vec_v);

        return grad;
    };

    std::vector<VectorXd> solutions(num_trials);
    std::vector<double>   objective_values(num_trials);

    auto process_performer = [&](const unsigned trial) {
        const VectorXd initial_cross_vec_u = internal::GetRandomSample(lower_bound, upper_bound);

        const VectorXd cross_vec_u = nloptutil::unconstrained::grad_based::bounded::solve(initial_cross_vec_u,
                                                                                          upper_bound,
                                                                                          lower_bound,
                                                                                          objective_func,
                                                                                          grad_func,
                                                                                          nlopt::LD_LBFGS,
                                                                                          true,
                                                                                          1000,
                                                                                          1e-06,
                                                                                          1e-06,
                                                                                          k_verbose_nlopt);

        solutions[trial]        = cross_vec_u;
        objective_values[trial] = objective_func(cross_vec_u);
    };

    parallelutil::parallel_for(num_trials, process_performer);

    const auto max_index =
        std::distance(objective_values.begin(), std::max_element(objective_values.begin(), objective_values.end()));
    const auto& cross_vec_v = solutions[max_index];

    // Calculate the position of the vertex opposite to x_ei while ensuring it is inside the bound
    const VectorXd x_ei_opposite = [&]() {
        constexpr unsigned max_num_iters = 100;
        constexpr double   factor        = 0.90;

        // This may be outside the bound
        VectorXd x = x_current_best - cross_vec_u;

        // Perform a naive line search
        unsigned iter = 0;
        while (k_ensure_opposite_vertex_bound)
        {
            if (internal::CheckIfParamIsInsideBound(x) || iter++ > max_num_iters)
            {
                break;
            }

            // Move towards the "safe" direction slightly
            x = factor * x + (1.0 - factor) * x_current_best;
        }

        // Now this is ensured to be inside the bound
        return x;
    }();

    return std::make_shared<PiecewisePlanarPlane>(
        x_current_best, x_ei, x_current_best + cross_vec_v, x_ei_opposite, x_current_best - cross_vec_v);
}
