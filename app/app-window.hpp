#ifndef APP_WINDOW_HPP
#define APP_WINDOW_HPP

#include <Eigen/Core>
#include <QMainWindow>
#include <memory>
#include <vector>

namespace sps
{
    class Optimizer;
}
namespace visopt
{
    class SlidersWidget;
}
class PreviewWidget;
class QLabel;
class QPushButton;
class ZoomableGridWidget;
class Recorder;

class AppWindow : public QMainWindow
{
public:
    AppWindow();

    void UpdateStatusMessage(const unsigned zoom_level, const unsigned num_zoom_levels);
    void StepPlaneSearch(const Eigen::VectorXd& x_preferred, const std::vector<Eigen::VectorXd>& x_others);

protected:
    void keyPressEvent(QKeyEvent* event) override;

private:
    void ExportRecords() const;

    std::shared_ptr<sps::Optimizer> m_optimizer;
    std::string                     m_working_dir_path;
    std::shared_ptr<Recorder>       m_recorder;

    ZoomableGridWidget*         m_zoomable_grid_widget;
    std::vector<PreviewWidget*> m_preview_widgets;
    QLabel*                     m_status_label;
    QPushButton*                m_study_button;
    QPushButton*                m_clock_reset_button;

    visopt::SlidersWidget* m_mu_sliders_widget;
    visopt::SlidersWidget* m_sigma_sliders_widget;
    visopt::SlidersWidget* m_ei_sliders_widget;
};

#endif // APP_WINDOW_HPP
