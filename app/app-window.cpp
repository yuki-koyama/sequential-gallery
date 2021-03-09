#include "app-window.hpp"
#include "preview-widget.hpp"
#include "recorder.hpp"
#include "zoomable-grid-widget.hpp"
#include <QDir>
#include <QFutureWatcher>
#include <QHBoxLayout>
#include <QKeyEvent>
#include <QLabel>
#include <QProgressDialog>
#include <QPushButton>
#include <QStatusBar>
#include <QtConcurrent>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sequential-line-search/sequential-line-search.hpp>
#include <sps/sps.hpp>
#include <string-util.hpp>
#include <thread>
#include <visoptslider/visoptslider.hpp>

constexpr bool k_is_verbose              = true;
constexpr bool k_show_sliders            = false;
constexpr bool k_show_additional_sliders = false;
#if defined(ZIN_APP_FOR_SMPL)
constexpr int k_default_width  = 620;
constexpr int k_default_height = 780;
#else
constexpr int  k_default_width      = 1080;
constexpr int  k_default_height     = 780;
#endif
#if defined(ZIN_APP_FOR_STUDY)
constexpr int  k_max_num_iterations = 15;
constexpr bool k_show_study_button  = true;
#else
constexpr int  k_max_num_iterations = 30;
constexpr bool k_show_study_button  = false;
#endif

QString GenerateStatusLabel(const unsigned zoom_level, const unsigned num_zoom_levels)
{
    QString message;

    message += R"(<div style="text-align: center;">)";

    for (unsigned i = 0; i < num_zoom_levels; ++i)
    {
        if (i == zoom_level)
        {
            message += R"(<span style="color: rgba(255, 255, 255, 0.8);">Zoom #)";
        }
        else
        {
            message += R"(<span style="color: rgba(255, 255, 255, 0.3);">Zoom #)";
        }

        message += QString::number(i + 1);
        message += R"(</span>)";

        if (i + 1 != num_zoom_levels)
        {
            message += R"(<span style="color: rgba(255, 255, 255, 0.3);"> &raquo; </span>)";
        }
    }

    return message;
}

// -----------------------------------------------------------------------------

AppWindow::AppWindow() : QMainWindow()
{
    // Set the central widget
    QWidget* central_widget = new QWidget();
    central_widget->setLayout(new QHBoxLayout());
    this->setCentralWidget(central_widget);

    // Set the grid widget
    m_zoomable_grid_widget = new ZoomableGridWidget(this);
    central_widget->layout()->addWidget(m_zoomable_grid_widget);

    // Set the visopt slider widget
    m_mu_sliders_widget = new visopt::SlidersWidget();
    m_mu_sliders_widget->initialize(
        PreviewWidget::GetNumDimensions(),
        [](const Eigen::VectorXd& x) { return 0.0; },
        PreviewWidget::GetUpperBound(),
        PreviewWidget::GetLowerBound(),
        1.0,
        0.0);
    m_mu_sliders_widget->setVisualizationMinimumSize(120, 12);
    m_mu_sliders_widget->setVisible(k_show_sliders);
    central_widget->layout()->addWidget(m_mu_sliders_widget);

    m_sigma_sliders_widget = new visopt::SlidersWidget();
    m_sigma_sliders_widget->initialize(
        PreviewWidget::GetNumDimensions(),
        [](const Eigen::VectorXd& x) { return 0.0; },
        PreviewWidget::GetUpperBound(),
        PreviewWidget::GetLowerBound(),
        1.0,
        0.0);
    m_sigma_sliders_widget->setVisualizationMinimumSize(120, 12);
    m_sigma_sliders_widget->setVisible(k_show_additional_sliders);
    central_widget->layout()->addWidget(m_sigma_sliders_widget);

    m_ei_sliders_widget = new visopt::SlidersWidget();
    m_ei_sliders_widget->initialize(
        PreviewWidget::GetNumDimensions(),
        [](const Eigen::VectorXd& x) { return 0.0; },
        PreviewWidget::GetUpperBound(),
        PreviewWidget::GetLowerBound(),
        1.0,
        0.0);
    m_ei_sliders_widget->setVisualizationMinimumSize(120, 12);
    m_ei_sliders_widget->setVisible(k_show_additional_sliders);
    central_widget->layout()->addWidget(m_ei_sliders_widget);

    if constexpr (k_show_study_button)
    {
        m_clock_reset_button = new QPushButton("Start the task!");
        m_study_button       = new QPushButton("I'm satisfied!");

        m_clock_reset_button->setFixedWidth(120);
        m_study_button->setFixedWidth(120);

        // Set a callback to the buttons
        QObject::connect(m_clock_reset_button, &QPushButton::clicked, [&]() {
            m_recorder->ResetClock();
            m_clock_reset_button->setDisabled(true);
        });
        QObject::connect(m_study_button, &QPushButton::clicked, [&]() {
            std::ofstream file_stream(m_working_dir_path + "/satisfied.txt");
            file_stream << std::to_string(m_recorder->GetNumRecords());
            m_study_button->setDisabled(true);
        });

        central_widget->layout()->addWidget(m_clock_reset_button);
        central_widget->layout()->addWidget(m_study_button);
    }

    // Set policies (c.f., https://doc.qt.io/qt-5/qsizepolicy.html)
    m_zoomable_grid_widget->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
    m_mu_sliders_widget->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    m_sigma_sliders_widget->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    m_ei_sliders_widget->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

    // Initialize an optimizer
    m_optimizer = std::make_shared<sps::Optimizer>(PreviewWidget::GetNumDimensions());

    // Set an initial search plane
    m_zoomable_grid_widget->SetPlane(m_optimizer->RetrieveSearchPlane());

    // Connet the optimizer and the visoptslider widget
    m_mu_sliders_widget->setTargetFunction(
        [this](const Eigen::VectorXd& x) { return m_optimizer->GetCurrentPredictiveMean(x); });
    m_sigma_sliders_widget->setTargetFunction(
        [this](const Eigen::VectorXd& x) { return m_optimizer->GetCurrentPredictiveVar(x); });
    m_ei_sliders_widget->setTargetFunction(
        [this](const Eigen::VectorXd& x) { return m_optimizer->GetCurrentExpectedImprovement(x); });

    // Set a status bar at the bottom of the window
    this->setStatusBar(new QStatusBar());
    m_status_label = new QLabel();
    this->statusBar()->addWidget(m_status_label, 1);
    UpdateStatusMessage(0, ZoomableGridWidget::GetNumZoomLevels());

    // Prepare a working directory
    m_working_dir_path = "./" + stringutil::GetCurrentTimeAsString();
    if (!QDir().mkdir(QString::fromStdString(m_working_dir_path)))
    {
        throw std::runtime_error("Could not create a working directory.");
    }

    // Prepare a recorder object
    m_recorder = std::make_shared<Recorder>("sequential-plane-search", m_working_dir_path);

    // Set the window size
    this->setGeometry(0, 0, k_default_width, k_default_height);
}

void AppWindow::UpdateStatusMessage(const unsigned zoom_level, const unsigned num_zoom_levels)
{
    m_status_label->setText(GenerateStatusLabel(zoom_level, num_zoom_levels));
    m_status_label->update();
}

void AppWindow::StepPlaneSearch(const Eigen::VectorXd& x_preferred, const std::vector<Eigen::VectorXd>& x_others)
{
    // Record the current plane
    m_recorder->AddPlaneSearchRecord(m_zoomable_grid_widget->GetPlane(), x_preferred);

    // Finish the app if the number of iterations reaches to the predefined number
    if (m_recorder->GetNumRecords() == k_max_num_iterations)
    {
        const auto export_process = [&]() { this->ExportRecords(); };

        QProgressDialog      dialog(QString("Exporting data..."), QString(), 0, 0, this);
        QFutureWatcher<void> watcher;
        QObject::connect(&watcher, SIGNAL(finished()), &dialog, SLOT(reset()));
        watcher.setFuture(QtConcurrent::run(export_process));
        dialog.exec();
        watcher.waitForFinished();

        exit(0);
    }

    // Note: This process will be run by another thread. If it takes only less than 500 ms, it will sleep for the
    // remaining time before it finishes.
    const auto search_process = [&]() {
        constexpr int minimum_duration = 500;

        const auto time_start = std::chrono::system_clock::now();

        m_optimizer->SubmitData(x_preferred, x_others);

        const auto time_end = std::chrono::system_clock::now();
        const auto elapsed  = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();

        std::cout << "AppWindow::StepPlaneSearch() ... data submission took " << elapsed << " ms" << std::endl;

        const int sleep_duration = std::max(minimum_duration - static_cast<int>(elapsed), 0);
        std::this_thread::sleep_for(std::chrono::milliseconds(sleep_duration));
    };

    // Detach the optimizer from the sliders widget to avoid conflict (TODO: Investigate better solutions)
    m_mu_sliders_widget->setTargetFunction([](const Eigen::VectorXd& x) { return 0.0; });
    m_sigma_sliders_widget->setTargetFunction([](const Eigen::VectorXd& x) { return 0.0; });
    m_ei_sliders_widget->setTargetFunction([](const Eigen::VectorXd& x) { return 0.0; });

    // Show a progress dialog during the data submission
    QProgressDialog      dialog(QString("Calculating..."), QString(), 0, 0, this);
    QFutureWatcher<void> watcher;
    QObject::connect(&watcher, SIGNAL(finished()), &dialog, SLOT(reset()));
    watcher.setFuture(QtConcurrent::run(search_process));
    dialog.exec();
    watcher.waitForFinished();

    const auto plane = m_optimizer->RetrieveSearchPlane();

    if constexpr (k_is_verbose)
    {
        const Eigen::IOFormat format(3);

        std::cout << "New plane: " << std::endl;
        std::cout << "- x0: " << plane->GetVertices()[0].transpose().format(format) << std::endl;
        std::cout << "- x1: " << plane->GetVertices()[1].transpose().format(format) << std::endl;
        std::cout << "- x2: " << plane->GetVertices()[2].transpose().format(format) << std::endl;
        std::cout << "- x3: " << plane->GetVertices()[3].transpose().format(format) << std::endl;
    }

    m_zoomable_grid_widget->SetPlane(plane);

    // Attach the optimizer again to the sliders widget
    m_mu_sliders_widget->setTargetFunction(
        [this](const Eigen::VectorXd& x) { return m_optimizer->GetCurrentPredictiveMean(x); });
    m_sigma_sliders_widget->setTargetFunction(
        [this](const Eigen::VectorXd& x) { return m_optimizer->GetCurrentPredictiveVar(x); });
    m_ei_sliders_widget->setTargetFunction(
        [this](const Eigen::VectorXd& x) { return m_optimizer->GetCurrentExpectedImprovement(x); });

    m_mu_sliders_widget->setArgumentAndUpdateSliders(plane->GetCenter());
    m_mu_sliders_widget->setMaximumValue(m_optimizer->GetCurrentPredictiveMean(plane->GetCenter()));
    m_mu_sliders_widget->setMinimumValue(0.0);

    m_sigma_sliders_widget->setArgumentAndUpdateSliders(plane->GetCenter());
    m_sigma_sliders_widget->setMaximumValue(m_optimizer->GetCurrentPredictiveVar(plane->GetVertices()[0]));
    m_sigma_sliders_widget->setMinimumValue(0.0);

    m_ei_sliders_widget->setArgumentAndUpdateSliders(plane->GetCenter());
    m_ei_sliders_widget->setMaximumValue(m_optimizer->GetCurrentExpectedImprovement(plane->GetVertices()[0]));
    m_ei_sliders_widget->setMinimumValue(0.0);
}

void AppWindow::keyPressEvent(QKeyEvent* event)
{
    switch (event->key())
    {
    case Qt::Key_Space:
    {
        const std::string output_path = m_working_dir_path + "/output.png";
        const unsigned    center      = (m_zoomable_grid_widget->GetGridResolution() - 1) / 2;
        const QImage      image       = m_zoomable_grid_widget->RenderImage(center, center);

        image.save(QString::fromStdString(output_path));

        break;
    }
    case Qt::Key_A:
    {
        for (unsigned row = 0; row < m_zoomable_grid_widget->GetGridResolution(); ++row)
        {
            for (unsigned col = 0; col < m_zoomable_grid_widget->GetGridResolution(); ++col)
            {
                const std::string output_path =
                    m_working_dir_path + "/output_" + std::to_string(row) + "_" + std::to_string(col) + ".png";
                const QImage image = m_zoomable_grid_widget->RenderImage(row, col);

                image.save(QString::fromStdString(output_path));
            }
        }

        break;
    }
    case Qt::Key_S:
    {
        if (event->modifiers() & Qt::ControlModifier)
        {
            this->ExportRecords();
        }

        break;
    }
    case Qt::Key_W:
    {
        this->setGeometry(0, 0, k_default_width, k_default_height);
    }
    default:
        break;
    }
}

void AppWindow::ExportRecords() const
{
    m_recorder->ExportCsv();
    m_recorder->ExportImages();
}
