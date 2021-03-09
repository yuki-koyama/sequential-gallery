#ifndef ZOOMABLE_GRID_WIDGET_HPP
#define ZOOMABLE_GRID_WIDGET_HPP

#include "preview-widget.hpp"
#include <QOpenGLWidget>
#include <memory>
#include <sps/plane.hpp>
#include <vector>

class AppWindow;
class PreviewWidget;
class AnimatedZoomableGridLayout;

class ZoomableGridWidget : public QOpenGLWidget
{
public:
    ZoomableGridWidget(AppWindow* app_window, QWidget* parent = nullptr);

    void                                      SetPlane(const std::shared_ptr<const sps::AbstractPlane> plane);
    std::shared_ptr<const sps::AbstractPlane> GetPlane() const { return m_plane; }

    QImage RenderImage(const unsigned row, const unsigned col) const;

    static unsigned GetGridResolution();
    static unsigned GetNumZoomLevels();

private:
    unsigned                                       m_zoom_level;
    std::vector<sps::AbstractPlane::GridCellIndex> m_selected_grid_cell_indices;
    std::shared_ptr<const sps::AbstractPlane>      m_plane;

    AppWindow*                               m_app_window;
    std::vector<std::vector<PreviewWidget*>> m_preview_widgets;
    AnimatedZoomableGridLayout*              m_layout;

    void SetNextZoomLevelParameters(const PreviewWidget::WidgetId& center_widget_id);

    void resizeEvent(QResizeEvent* event) override;
};

#endif // ZOOMABLE_GRID_WIDGET_HPP
