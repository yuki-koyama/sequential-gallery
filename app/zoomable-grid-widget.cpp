#include "zoomable-grid-widget.hpp"
#include "animated-zoomable-grid-layout.hpp"
#include "app-window.hpp"
#include <QTimer>
#include <iostream>

#if defined(ZIN_APP_FOR_SMPL)
constexpr unsigned k_grid_resolution = 3;
#else
constexpr unsigned k_grid_resolution = 5;
#endif
constexpr unsigned k_num_zoom_levels    = 4;
constexpr double   k_inter_level_scale  = 0.50;
constexpr double   k_animation_duration = 0.75;

unsigned ZoomableGridWidget::GetGridResolution() { return k_grid_resolution; }
unsigned ZoomableGridWidget::GetNumZoomLevels() { return k_num_zoom_levels; }

// ---------------------------------------------------------------------------------------------------------------------

ZoomableGridWidget::ZoomableGridWidget(AppWindow* app_window, QWidget* parent)
    : QOpenGLWidget(parent), m_zoom_level(0), m_app_window(app_window)
{
    m_layout = new AnimatedZoomableGridLayout(k_grid_resolution, k_animation_duration);

    m_preview_widgets = std::vector<std::vector<PreviewWidget*>>(
        k_grid_resolution, std::vector<PreviewWidget*>(k_grid_resolution, nullptr));
    for (unsigned i = 0; i < k_grid_resolution; ++i)
    {
        for (unsigned j = 0; j < k_grid_resolution; ++j)
        {
            PreviewWidget* preview_widget = new PreviewWidget({i, j});

            preview_widget->SetMousePressCallback([&](const PreviewWidget::WidgetId& widget_id) {
                if (m_zoom_level + 1 < k_num_zoom_levels)
                {
                    // Zoom up the space and continue user interaction
                    this->SetNextZoomLevelParameters(widget_id);

                    // Update the state
                    m_zoom_level++;
                    assert(m_zoom_level == m_selected_grid_cell_indices.size());

                    // Update the app window
                    m_app_window->UpdateStatusMessage(m_zoom_level, k_num_zoom_levels);
                }
                else
                {
                    // Ask the app to proceed to the next step
                    const Eigen::VectorXd& x_preferred =
                        m_preview_widgets[std::get<0>(widget_id)][std::get<1>(widget_id)]->GetParameters();

                    const Eigen::VectorXd&                center   = m_plane->GetCenter();
                    const std::array<Eigen::VectorXd, 4>& vertices = m_plane->GetVertices();

                    m_app_window->StepPlaneSearch(x_preferred,
                                                  {center, vertices[0], vertices[1], vertices[2], vertices[3]});

                    // Update the state
                    m_zoom_level = 0;

                    // Remove the buffer for storing user selections
                    m_selected_grid_cell_indices.clear();

                    // Update the app window
                    m_app_window->UpdateStatusMessage(m_zoom_level, k_num_zoom_levels);
                }
            });

            m_layout->addWidget(preview_widget, i, j);

            m_preview_widgets[i][j] = preview_widget;
        }
    }

    // Note: setLayout() should be called after widgets are added (TODO: find the reference taking about this).
    this->setLayout(m_layout);
}

void ZoomableGridWidget::SetPlane(const std::shared_ptr<const sps::AbstractPlane> plane)
{
    const unsigned radius = (k_grid_resolution - 1) / 2;

    for (unsigned i = 0; i < k_grid_resolution; ++i)
    {
        for (unsigned j = 0; j < k_grid_resolution; ++j)
        {
            // Calculate the parameters at the point (i - radius, j - radius)
            const sps::AbstractPlane::GridCellIndex grid_cell = {
                static_cast<int>(i) - static_cast<int>(radius),
                static_cast<int>(j) - static_cast<int>(radius),
            };

            const Eigen::VectorXd x = plane->CalcGridParameters(grid_cell, k_grid_resolution, k_inter_level_scale, {});

            m_preview_widgets[i][j]->SetParameters(x);
        }
    }

    m_plane = plane;
}

QImage ZoomableGridWidget::RenderImage(const unsigned row, const unsigned col) const
{
    return PreviewWidget::RenderImage(m_preview_widgets[row][col]->GetParameters());
}

void ZoomableGridWidget::SetNextZoomLevelParameters(const PreviewWidget::WidgetId& center_widget_id)
{
    std::cout << "GridWidget::ZoomPlane(): Central widget = (";
    std::cout << std::get<0>(center_widget_id) << ", " << std::get<1>(center_widget_id);
    std::cout << ")" << std::endl;

    assert(k_grid_resolution % 2 == 1);

    const int radius = (k_grid_resolution - 1) / 2;

    const auto& [center_row, center_col] = center_widget_id;

    m_selected_grid_cell_indices.push_back({static_cast<int>(std::get<0>(center_widget_id)) - radius,
                                            static_cast<int>(std::get<1>(center_widget_id)) - radius});

    // Layout animation
    const Eigen::Vector2d center_position = m_layout->getItemWrapper(center_row, center_col).target_position;
    const Eigen::Vector2d global_offset   = m_layout->getItemWrapper(radius, radius).target_position - center_position;
    const Eigen::Vector2d unit_offset     = {
        m_layout->getItemWrapper(0, 1).target_position(0) - m_layout->getItemWrapper(0, 0).target_position(0),
        m_layout->getItemWrapper(1, 0).target_position(1) - m_layout->getItemWrapper(0, 0).target_position(1)};

    for (int row = 0; row < k_grid_resolution; ++row)
    {
        for (int col = 0; col < k_grid_resolution; ++col)
        {
            auto& item_wrapper = m_layout->getItemWrapper(row, col);

            const Eigen::Vector2d local_offset =
                (1.0 / k_inter_level_scale) *
                Eigen::Vector2d{(static_cast<double>(col) - static_cast<double>(center_col)) * unit_offset(0),
                                (static_cast<double>(row) - static_cast<double>(center_row)) * unit_offset(1)};

            item_wrapper.initial_position = item_wrapper.current_position;
            item_wrapper.target_position  = center_position + global_offset + local_offset;
            item_wrapper.initial_scale    = 1.0;
            item_wrapper.target_scale     = 1.0;
            item_wrapper.current_scale    = 1.0;
            item_wrapper.initial_time     = std::chrono::system_clock::now();
        }
    }

    const int duration_in_millisecond = static_cast<int>(m_layout->getDuration() * 1000.0);

    QTimer::singleShot(duration_in_millisecond, [&]() {
        const QRect rect = m_layout->geometry();

        std::vector<Eigen::Vector2d> inside_cell_positions;
        for (int row = 0; row < k_grid_resolution; ++row)
        {
            for (int col = 0; col < k_grid_resolution; ++col)
            {
                const Eigen::Vector2d& position = m_layout->getItemWrapper(row, col).target_position;

                // Determine if the grid cell is inside the viewing region (i.e., currently visible)
                const bool is_inside =
                    position.minCoeff() >= 0.0 && position(0) <= rect.width() && position(1) <= rect.height();

                // If it is inside, remember the position
                if (is_inside)
                {
                    inside_cell_positions.push_back(position);
                }
            }
        }

        // Reset all the widgets to the original grid positions
        m_layout->setWidgetTargetPositionsUsingGridIndex();
        m_layout->setWidgetCurrentPositionsToTargetPositions();

        // Set scaling animations
        for (int row = 0; row < k_grid_resolution; ++row)
        {
            for (int col = 0; col < k_grid_resolution; ++col)
            {
                auto& item_wrapper = m_layout->getItemWrapper(row, col);

                item_wrapper.initial_position = item_wrapper.current_position;
                item_wrapper.target_position  = item_wrapper.current_position;
                item_wrapper.target_scale     = 1.0;
                item_wrapper.initial_scale    = 0.0;
                item_wrapper.current_scale    = 0.0;
                item_wrapper.initial_time     = std::chrono::system_clock::now();

                // Unset animations if the corresponding visuals are originally visible
                for (const Eigen::Vector2d& position : inside_cell_positions)
                {
                    if (m_layout->getItemWrapper(row, col).target_position.isApprox(position))
                    {
                        m_layout->getItemWrapper(row, col).initial_scale = 1.0;
                        m_layout->getItemWrapper(row, col).current_scale = 1.0;
                    }
                }
            }
        }

        // Preview widgets
        for (int i = 0; i < k_grid_resolution; ++i)
        {
            for (int j = 0; j < k_grid_resolution; ++j)
            {
                const sps::AbstractPlane::GridCellIndex grid_cell = {
                    static_cast<int>(i) - static_cast<int>(radius),
                    static_cast<int>(j) - static_cast<int>(radius),
                };

                // Calculate the parameters at the point (i - radius, j - radius)
                const Eigen::VectorXd x = m_plane->CalcGridParameters(
                    grid_cell, k_grid_resolution, k_inter_level_scale, m_selected_grid_cell_indices);

                m_preview_widgets[i][j]->SetParameters(x);
            }
        }

        // Note: It is important to call this update() here to avoid visual clutter. When the preview widgets are based
        // on QOpenGLWidget, this way seems to work appropriately. When not, it is probably necessary to handle the
        // update order in a different way.
        m_layout->update();
    });
}

void ZoomableGridWidget::resizeEvent(QResizeEvent* event) { m_layout->setWidgetTargetPositionsUsingGridIndex(); }
