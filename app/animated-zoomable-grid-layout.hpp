#ifndef ANIMATED_ZOOMABLE_GRID_LAYOUT_HPP
#define ANIMATED_ZOOMABLE_GRID_LAYOUT_HPP

#include <Eigen/Core>
#include <QLayout>
#include <chrono>
#include <memory>
#include <vector>

class QTimer;

class AnimatedZoomableGridLayout : public QLayout
{
public:
    struct ItemWrapper
    {
        QLayoutItem*                          item = nullptr;
        unsigned                              row;
        unsigned                              col;
        Eigen::Vector2d                       target_position;
        Eigen::Vector2d                       current_position;
        Eigen::Vector2d                       initial_position;
        double                                target_scale;
        double                                current_scale;
        double                                initial_scale;
        std::chrono::system_clock::time_point initial_time;
    };

    AnimatedZoomableGridLayout(const unsigned grid_resolution, const double duration);
    ~AnimatedZoomableGridLayout();

    void addWidget(QWidget* widget, const unsigned row, const unsigned col);
    void setWidgetTargetPositionsUsingGridIndex();
    void setWidgetCurrentPositionsToTargetPositions();
    void setRandomTargetPositions();

    ItemWrapper& getItemWrapper(const unsigned row, const unsigned col);

    void         addItem(QLayoutItem* item) override;
    int          count() const override;
    QLayoutItem* itemAt(int index) const override;
    QLayoutItem* takeAt(int index) override;

    QSize sizeHint() const override;
    void  setGeometry(const QRect& rect) override;

    void play();
    void stop();

    double getDuration() const;

private:
    std::vector<std::vector<ItemWrapper>> m_item_wrappers;

    std::chrono::system_clock::time_point m_previous_time;
    std::chrono::system_clock::time_point m_origin_time;

    std::unique_ptr<QTimer> m_timer;

    const unsigned m_grid_resolution;
    const double   m_duration;
};

#endif // ANIMATED_ZOOMABLE_GRID_LAYOUT_HPP
