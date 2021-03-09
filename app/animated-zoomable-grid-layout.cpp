#include "animated-zoomable-grid-layout.hpp"
#include <QTimer>
#include <QWidget>
#include <cmath>

constexpr int    k_target_fps      = 60;
constexpr int    k_target_interval = 1000 / k_target_fps;
constexpr double k_padding         = 10.0;

// ---------------------------------------------------------------------------------------------------------------------

namespace internal
{
    inline std::pair<double, double>
    CalcSizePerWidget(const QRect& rect, const double padding, const unsigned grid_resolution)
    {
        const double padding_sum = padding * static_cast<double>(grid_resolution + 1);
        const double width_per_widget =
            (static_cast<double>(rect.width()) - padding_sum) / static_cast<double>(grid_resolution);
        const double height_per_widget =
            (static_cast<double>(rect.height()) - padding_sum) / static_cast<double>(grid_resolution);

        return {width_per_widget, height_per_widget};
    }

    inline Eigen::Vector2d CalcTargetCenterPosition(const double   padding,
                                                    const double   width_per_widget,
                                                    const double   height_per_widget,
                                                    const unsigned grid_resolution,
                                                    const int      row,
                                                    const int      col)
    {
        return Eigen::Vector2d{(col + 1) * padding + (col + 0.5) * width_per_widget,
                               (row + 1) * padding + (row + 0.5) * height_per_widget};
    }

    template <typename T> T ApplyEasing(const T& init, const T& target, const double t, const double duration)
    {
        const double weight        = std::min(t / duration, 1.0);
        const double warped_weight = 1.0 - std::pow(1.0 - weight, 5.0);

        return warped_weight * target + (1.0 - warped_weight) * init;
    }
} // namespace internal

// ---------------------------------------------------------------------------------------------------------------------

AnimatedZoomableGridLayout::AnimatedZoomableGridLayout(const unsigned grid_resolution, const double duration)
    : m_previous_time(std::chrono::system_clock::now()), m_origin_time(std::chrono::system_clock::now()),
      m_grid_resolution(grid_resolution), m_duration(duration)
{
    m_timer = std::make_unique<QTimer>();

    // Note: With Qt 5.12+, we can also implement the same effect by
    // m_timer->callOnTimeout([=]() { this->update(); });
    QObject::connect(m_timer.get(), &QTimer::timeout, [=]() { this->update(); });

    // Set the animation resolution
    m_timer->setInterval(k_target_interval);

    // Make the animation active by default
    this->play();

    // Allocate grid items
    m_item_wrappers =
        std::vector<std::vector<ItemWrapper>>(m_grid_resolution, std::vector<ItemWrapper>(m_grid_resolution));
}

AnimatedZoomableGridLayout::~AnimatedZoomableGridLayout()
{
    QLayoutItem* item = takeAt(0);
    while (item != nullptr)
    {
        delete item;

        item = takeAt(0);
    }
}

// ---------------------------------------------------------------------------------------------------------------------

void AnimatedZoomableGridLayout::addWidget(QWidget* widget, const unsigned row, const unsigned col)
{
    m_item_wrappers[row][col] = {new QWidgetItem(widget),
                                 row,
                                 col,
                                 Eigen::Vector2d(),
                                 Eigen::Vector2d(),
                                 Eigen::Vector2d(),
                                 1.0,
                                 1.0,
                                 1.0,
                                 std::chrono::system_clock::now()};
}

void AnimatedZoomableGridLayout::setWidgetTargetPositionsUsingGridIndex()
{
    const QRect rect = this->geometry();

    const double padding_sum = k_padding * static_cast<double>(m_grid_resolution + 1);
    const double width_per_widget =
        (static_cast<double>(rect.width()) - padding_sum) / static_cast<double>(m_grid_resolution);
    const double height_per_widget =
        (static_cast<double>(rect.height()) - padding_sum) / static_cast<double>(m_grid_resolution);

    for (auto& col : m_item_wrappers)
    {
        for (auto& item_wrapper : col)
        {
            item_wrapper.target_position = internal::CalcTargetCenterPosition(
                k_padding, width_per_widget, height_per_widget, m_grid_resolution, item_wrapper.row, item_wrapper.col);
        }
    }
}

void AnimatedZoomableGridLayout::setWidgetCurrentPositionsToTargetPositions()
{
    for (auto& col : m_item_wrappers)
    {
        for (auto& item_wrapper : col)
        {
            item_wrapper.current_position = item_wrapper.target_position;
        }
    }
}

void AnimatedZoomableGridLayout::setRandomTargetPositions()
{
    const auto random_position_generator = [](const int width, const int height) -> Eigen::Vector2d {
        return (0.5 * (Eigen::Array2d::Random() + Eigen::Array2d::Ones()) * Eigen::Array2d(width, height)).matrix();
    };

    for (auto& col : m_item_wrappers)
    {
        for (auto& item_wrapper : col)
        {
            item_wrapper.target_position =
                random_position_generator(this->geometry().width(), this->geometry().height());
        }
    }
}

AnimatedZoomableGridLayout::ItemWrapper& AnimatedZoomableGridLayout::getItemWrapper(const unsigned row,
                                                                                    const unsigned col)
{
    assert(m_item_wrappers[row][col].item != nullptr);
    assert(row < m_grid_resolution && col < m_grid_resolution);

    return m_item_wrappers[row][col];
}

// ---------------------------------------------------------------------------------------------------------------------

void AnimatedZoomableGridLayout::addItem(QLayoutItem* item)
{
    throw std::runtime_error("AnimatedZoomableGridLayout::addItem() should not be called.");
}

int AnimatedZoomableGridLayout::count() const
{
    unsigned count = 0;
    for (const auto& col : m_item_wrappers)
    {
        for (const auto& item_wrapper : col)
        {
            if (item_wrapper.item != nullptr)
            {
                ++count;
            }
        }
    }
    return count;
}

QLayoutItem* AnimatedZoomableGridLayout::itemAt(int index) const
{
    if (index < 0 || index >= this->count())
    {
        return nullptr;
    }

    unsigned count = 0;
    for (const auto& col : m_item_wrappers)
    {
        for (const auto& item_wrapper : col)
        {
            if (item_wrapper.item != nullptr)
            {
                if (count == index)
                {
                    return item_wrapper.item;
                }
                ++count;
            }
        }
    }

    assert(false);
    return nullptr;
}

QLayoutItem* AnimatedZoomableGridLayout::takeAt(int index)
{
    if (index < 0 || index >= this->count())
    {
        return nullptr;
    }

    unsigned count = 0;
    for (auto& col : m_item_wrappers)
    {
        for (auto& item_wrapper : col)
        {
            if (item_wrapper.item != nullptr)
            {
                if (count == index)
                {
                    QLayoutItem* item = item_wrapper.item;

                    item_wrapper.item = nullptr;

                    return item;
                }

                ++count;
            }
        }
    }

    assert(false);
    return nullptr;
}

// ---------------------------------------------------------------------------------------------------------------------

QSize AnimatedZoomableGridLayout::sizeHint() const { return QSize(600, 400); }

void AnimatedZoomableGridLayout::setGeometry(const QRect& rect)
{
    QLayout::setGeometry(rect);

    const std::chrono::system_clock::time_point current_time = std::chrono::system_clock::now();

    m_previous_time = current_time;

    const auto& [width_per_widget, height_per_widget] = internal::CalcSizePerWidget(rect, k_padding, m_grid_resolution);

    for (auto& col : m_item_wrappers)
    {
        for (auto& item_wrapper : col)
        {
            const double t =
                static_cast<std::chrono::duration<double>>(current_time - item_wrapper.initial_time).count();

            const Eigen::Vector2d new_position =
                internal::ApplyEasing(item_wrapper.initial_position, item_wrapper.target_position, t, m_duration);
            const double new_scale =
                internal::ApplyEasing(item_wrapper.initial_scale, item_wrapper.target_scale, t, m_duration);

            const double clamped_new_scale = std::max(new_scale, 0.0);

            const QRect geom(rect.x() - 0.5 * clamped_new_scale * width_per_widget + new_position(0),
                             rect.y() - 0.5 * clamped_new_scale * height_per_widget + new_position(1),
                             clamped_new_scale * width_per_widget,
                             clamped_new_scale * height_per_widget);

            item_wrapper.item->setGeometry(geom);
            item_wrapper.current_position = new_position;
            item_wrapper.current_scale    = new_scale;
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------------

void AnimatedZoomableGridLayout::play() { m_timer->start(); }

void AnimatedZoomableGridLayout::stop() { m_timer->stop(); }

// ---------------------------------------------------------------------------------------------------------------------

double AnimatedZoomableGridLayout::getDuration() const { return m_duration; }
