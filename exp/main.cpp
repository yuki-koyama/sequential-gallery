#include <Eigen/LU>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <mathtoolbox/probability-distributions.hpp>
#include <optimization-test-functions.hpp>
#include <parallel-util.hpp>
#include <sequential-line-search/sequential-line-search.hpp>
#include <sps/sps.hpp>
#include <timer.hpp>

Eigen::MatrixXd AppendStats(const Eigen::MatrixXd& data)
{
    Eigen::MatrixXd appended_data(data.rows(), data.cols() + 5);

    appended_data.block(0, 1, data.rows(), data.cols()) = data;

    const unsigned iters_col = 0;
    const unsigned mean_col  = data.cols() + 1;
    const unsigned sd_col    = data.cols() + 2;
    const unsigned upper_col = data.cols() + 3;
    const unsigned lower_col = data.cols() + 4;

    // Iterations
    for (unsigned row = 0; row < data.rows(); ++row)
    {
        appended_data(row, iters_col) = row + 1;
    }

    // Mean
    for (unsigned row = 0; row < data.rows(); ++row)
    {
        appended_data(row, mean_col) = data.row(row).mean();
    }

    // Standard deviation (biased)
    for (unsigned row = 0; row < data.rows(); ++row)
    {
        const double         mean     = appended_data(row, mean_col);
        const Eigen::ArrayXd data_row = data.row(row).array();
        const double         var =
            (data_row - Eigen::ArrayXd::Constant(data.cols(), mean)).square().sum() / static_cast<double>(data.cols());

        appended_data(row, sd_col) = std::sqrt(var);
    }

    // Mean plus/minus 1.0 * SD
    appended_data.col(upper_col) = appended_data.col(mean_col) + appended_data.col(sd_col);
    appended_data.col(lower_col) = appended_data.col(mean_col) - appended_data.col(sd_col);

    return appended_data;
}

void ExportCsv(const Eigen::MatrixXd& data, const std::string& path)
{
    const unsigned num_trials = data.cols();

    const Eigen::IOFormat csv_format(Eigen::StreamPrecision, Eigen::DontAlignCols, ",");
    std::ofstream         residuals_stream(path);

    residuals_stream << "#Iters,";
    for (unsigned trial = 0; trial < num_trials; ++trial)
    {
        residuals_stream << "Trial #" << std::to_string(trial + 1) << ",";
    }
    residuals_stream << "Mean,SD,Upper,Lower\n";
    residuals_stream << AppendStats(data).format(csv_format);
}

// ---------------------------------------------------------------------------------------------------------------------

std::tuple<std::string, std::function<double(const Eigen::VectorXd&)>, Eigen::VectorXd>
GetRosenbrockSetting(const unsigned num_dims)
{
    const Eigen::VectorXd x_optimal   = 0.25 * otf::GetSolution(num_dims, otf::FunctionType::Rosenbrock);
    const auto            latent_func = [&](const Eigen::VectorXd& x) {
        return -otf::GetValue(4.0 * x, otf::FunctionType::Rosenbrock);
    };
    return {"rosenbrock-" + std::to_string(num_dims) + "d", latent_func, x_optimal};
}

std::tuple<std::string, std::function<double(const Eigen::VectorXd&)>, Eigen::VectorXd>
GetIsometricSetting(const unsigned num_dims)
{
    const Eigen::VectorXd x_optimal = Eigen::VectorXd::Constant(num_dims, 0.30);
    const auto latent_func = [=](const Eigen::VectorXd& x) { return std::exp(-(x - x_optimal).squaredNorm()); };
    return {"isometric-" + std::to_string(num_dims) + "d", latent_func, x_optimal};
}

// ---------------------------------------------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    // Constants
    constexpr bool export_areas        = false;
    constexpr bool use_map_hyperparams = true;
#if defined(ZIN_EXPERIMENT_SCALE_MINIMUM)
    constexpr unsigned num_trials = 1;
    constexpr unsigned num_iters  = 5;
#elif defined(ZIN_EXPERIMENT_SCALE_MASSIVE)
    constexpr unsigned num_trials          = 50;
    constexpr unsigned num_iters           = 25;
#else
    constexpr unsigned num_trials = 5;
    constexpr unsigned num_iters  = 10;
#endif
    constexpr unsigned num_candidates        = 5;
    constexpr unsigned num_zoom_levels       = 8;
    constexpr double   inter_level_scale     = 0.50;
    const double       finest_level_scale    = std::pow(inter_level_scale, num_zoom_levels - 1);
    const double       finest_level_interval = finest_level_scale / static_cast<double>(num_candidates);

#if defined(ZIN_EXPERIMENT_USE_LAST_SELECTION_STRATEGY)
    constexpr auto sps_center_strategy = sps::PlaneCenterSelectionStrategy::LastSelection;
    constexpr auto sls_slider_strategy = sequential_line_search::SliderEndSelectionStrategy::LastSelection;
#else
    constexpr auto     sps_center_strategy = sps::PlaneCenterSelectionStrategy::LargestExpectValue;
    constexpr auto     sls_slider_strategy = sequential_line_search::SliderEndSelectionStrategy::LargestExpectValue;
#endif

    // Random seed
    std::srand(std::time(nullptr));

    // Synthetic functions
    const std::vector<std::tuple<std::string, std::function<double(const Eigen::VectorXd& x)>, Eigen::VectorXd>>
        settings = {
#if defined(ZIN_EXPERIMENT_SCALE_MINIMUM)
            GetRosenbrockSetting(5),
            GetIsometricSetting(5),
#elif defined(ZIN_EXPERIMENT_SCALE_MASSIVE)
            GetRosenbrockSetting(5),
            GetRosenbrockSetting(10),
            GetRosenbrockSetting(15),
            GetRosenbrockSetting(20),
            GetIsometricSetting(5),
            GetIsometricSetting(10),
            GetIsometricSetting(15),
            GetIsometricSetting(20),
#else
            GetRosenbrockSetting(10),
            GetIsometricSetting(10),
#endif
        };

    // -----------------------------------------------------------------------------------------------------------------

    for (const auto& setting : settings)
    {
        const std::string&     dir_name    = std::get<0>(setting);
        const auto&            latent_func = std::get<1>(setting);
        const Eigen::VectorXd& x_optimal   = std::get<2>(setting);
        const unsigned         num_dims    = x_optimal.size();

        std::cout << dir_name << std::endl;

        // Create a directory for exporting reports
        if (int result = std::system(std::string("mkdir -p ./" + dir_name).c_str()); result < 0)
        {
            throw std::runtime_error("Failed to create a directory.");
        }

        if constexpr (true)
        {
            // Strategies
            const auto init_strategy = std::make_shared<sps::FixedCenterFreeCrossVecsInitStrategy>();
            const std::vector<std::pair<std::string, std::shared_ptr<sps::PlaneConstructionStrategy>>>
                construction_strategies = {
                    {"first-ei-then-plane-integral", std::make_shared<sps::FirstEiThenPlaneIntegralStrategy>()},
                    {"fixed-center-free-cross-vecs", std::make_shared<sps::FixedCenterFreeCrossVecsStrategy>()},
#if false
                    {"joint-cross-vec-estimation", std::make_shared<sps::JointCrossVecEstimationStrategy>()},
                    {"batch-ie-piecewise-planar", std::make_shared<sps::BatchEiPiecewisePlanarStrategy>()},
                    {"approx-integral-joint-cross-vec-estimation",
                     std::make_shared<sps::ApproxIntegralJointCrossVecEstimationStrategy>()},
#endif
                };

            for (const auto& [strategy_name, construction_strategy] : construction_strategies)
            {
                // Results
                Eigen::MatrixXd residuals(num_iters, num_trials);
                Eigen::MatrixXd values(num_iters, num_trials);
                Eigen::MatrixXd areas(num_iters, num_trials);

                // Trials
                const auto trial_process = [&,
                                            strategy_name         = strategy_name,
                                            construction_strategy = construction_strategy](const unsigned trial) {
                    timer::Timer t("sps (" + strategy_name + ") - trial #" + std::to_string(trial + 1));

                    sps::Optimizer optimizer(
                        num_dims, use_map_hyperparams, init_strategy, construction_strategy, sps_center_strategy);

                    for (unsigned iter = 0; iter < num_iters; ++iter)
                    {
                        const auto plane =
                            std::static_pointer_cast<const sps::AbstractPlane>(optimizer.RetrieveSearchPlane());

                        Eigen::VectorXd x_best;
                        double          f_best = std::numeric_limits<double>::lowest();

                        for (double w_i = 0.0; w_i <= 1.0; w_i += finest_level_interval)
                        {
                            for (double w_j = 0.0; w_j <= 1.0; w_j += finest_level_interval)
                            {
                                const Eigen::Vector2d grid_coords{2.0 * w_i - 1.0, 2.0 * w_j - 1.0};

                                const Eigen::VectorXd x = plane->CalcParameters(grid_coords);
                                const double          f = latent_func(x);

                                if (f_best < f)
                                {
                                    x_best = x;
                                    f_best = f;
                                }
                            }
                        }

                        residuals(iter, trial) = (x_best - x_optimal).norm();
                        values(iter, trial)    = f_best;
                        areas(iter, trial)     = plane->CalcArea();

                        optimizer.SubmitData(x_best,
                                             {plane->GetCenter(),
                                              plane->GetVertices()[0],
                                              plane->GetVertices()[1],
                                              plane->GetVertices()[2],
                                              plane->GetVertices()[3]});
                    }
                };

                parallelutil::queue_based_parallel_for(num_trials, trial_process);

                // File IO
                ExportCsv(residuals, "./" + dir_name + "/sps-" + strategy_name + "-residual.csv");
                ExportCsv(values, "./" + dir_name + "/sps-" + strategy_name + "-value.csv");
                if constexpr (export_areas)
                {
                    ExportCsv(areas, "./" + dir_name + "/sps-" + strategy_name + "-area.csv");
                }
            }
        }

        // -------------------------------------------------------------------------------------------------------------

        if constexpr (true)
        {
            // Strategies
            constexpr auto init_strategy          = sequential_line_search::GenerateCenteredFixedLengthRandomSliderEnds;
            constexpr bool use_slider_enlargement = false;
            constexpr bool use_map_hyperparams    = false;

            // Results
            Eigen::MatrixXd residuals(num_iters, num_trials);
            Eigen::MatrixXd values(num_iters, num_trials);

            // Trials
            const auto trial_process = [&](const unsigned trial) {
                timer::Timer t("sls - trial #" + std::to_string(trial + 1));

                sequential_line_search::SequentialLineSearchOptimizer optimizer(
                    num_dims,
                    use_slider_enlargement,
                    use_map_hyperparams,
                    sequential_line_search::KernelType::ArdMatern52Kernel,
                    sequential_line_search::AcquisitionFuncType::ExpectedImprovement,
                    init_strategy,
                    sls_slider_strategy);

                for (unsigned iter = 0; iter < num_iters; ++iter)
                {
                    const auto [x_s, x_d] = optimizer.GetSliderEnds();

                    double          w_best;
                    Eigen::VectorXd x_best;
                    double          f_best = std::numeric_limits<double>::lowest();

                    for (double w = 0.0; w <= 1.0; w += finest_level_interval)
                    {
                        const Eigen::VectorXd x = (1.0 - w) * x_s + w * x_d;
                        const double          f = latent_func(x);

                        if (f_best < f)
                        {
                            w_best = w;
                            x_best = x;
                            f_best = f;
                        }
                    }

                    residuals(iter, trial) = (x_best - x_optimal).norm();
                    values(iter, trial)    = f_best;

                    optimizer.SubmitLineSearchResult(w_best);
                }
            };

            parallelutil::queue_based_parallel_for(num_trials, trial_process);

            // File IO
            ExportCsv(residuals, "./" + dir_name + "/sls-residual.csv");
            ExportCsv(values, "./" + dir_name + "/sls-value.csv");
        }
    }

    // -----------------------------------------------------------------------------------------------------------------

    return 0;
}
