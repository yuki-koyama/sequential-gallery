#ifndef PREVIEW_WIDGET_HPP
#define PREVIEW_WIDGET_HPP

#include <Eigen/Core>
#include <QWidget>
#include <functional>
#include <utility>

namespace enhancer
{
    class EnhancerWidget;
}
class OverlayWidget;

class PreviewWidget : public QWidget
{
public:
    using WidgetId = std::pair<unsigned, unsigned>;

    PreviewWidget(const WidgetId& widget_id, QWidget* parent = nullptr);

    const Eigen::VectorXd& GetParameters() const;
    void                   SetParameters(const Eigen::VectorXd& parameters, const bool update_widget = true);
    void                   SetMousePressCallback(const std::function<void(const WidgetId&)>& callback);

    static QImage RenderImage(const Eigen::VectorXd& parameters);

    static unsigned               GetNumDimensions();
    static const Eigen::VectorXd& GetLowerBound();
    static const Eigen::VectorXd& GetUpperBound();

protected:
    void resizeEvent(QResizeEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseReleaseEvent(QMouseEvent* event) override;

private:
    const WidgetId m_widget_id;

    Eigen::VectorXd m_parameters;

    std::function<void(const WidgetId&)> m_mouse_press_callback;

    enhancer::EnhancerWidget* m_widget;

    OverlayWidget* m_overlay;
};

#endif // PREVIEW_WIDGET_HPP
