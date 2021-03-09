#include "preview-widget.hpp"
#include <QHBoxLayout>
#include <QPainter>
#include <QResizeEvent>
#include <QTimer>
#include <enhancer/enhancerwidget.hpp>

constexpr bool k_display_some_values = false;

namespace internal
{
    const QImage& GetImageSingleton()
    {
        constexpr const char* paths[] = {
            "./resources/scaled/20170923-DSC06498.jpg", // 0
            "./resources/scaled/20180721-DSC01398.jpg", // 1
            "./resources/scaled/DSC03039.JPG",          // 2
        };
        constexpr auto path = paths[2]; // Edit the index to change the target photograph

        static QImage image(path);

        return image;
    }

    Eigen::VectorXd TransformParams(const Eigen::VectorXd& params)
    {
        constexpr double scale = 1.0 / 3.0;
        return Eigen::VectorXd::Constant(params.size(), 0.5 * (1.0 - scale)) + scale * params;
    }
} // namespace internal

class OverlayWidget : public QWidget
{
public:
    OverlayWidget(PreviewWidget* preview_widget, QWidget* parent = nullptr)
        : QWidget(parent), m_preview_widget(preview_widget), m_hovered(false), m_pressed(false),
          m_disallow_animation(false), m_border_color(20, 20, 20, 160), m_border_width(5.0)
    {
        this->setAttribute(Qt::WA_Hover, true);
    }

protected:
    void paintEvent(QPaintEvent* event)
    {
        const QPen pen = QPen(QBrush(m_border_color), m_border_width);

        QPainter painter;

        painter.begin(this);

        painter.setPen(pen);
        painter.drawRect(this->rect());

        if constexpr (k_display_some_values)
        {
            // Display some parameter values for debug
            const QString message = QString("x[0] = ") + QString::number(m_preview_widget->GetParameters()(0), 'f', 3) +
                                    QString("\n") + QString("x[1] = ") +
                                    QString::number(m_preview_widget->GetParameters()(1), 'f', 3);

            painter.drawText(this->rect(), Qt::AlignCenter, message);
        }

        painter.end();
    }

    bool event(QEvent* e)
    {
        switch (e->type())
        {
        case QEvent::HoverEnter:
        {
            m_hovered            = true;
            m_disallow_animation = true;
            m_border_color       = QColor(20, 250, 100, 200);

            this->update();

            break;
        }
        case QEvent::HoverLeave:
        {
            m_hovered            = false;
            m_disallow_animation = false;

            const auto update_procedure = [&]() {
                if (m_disallow_animation)
                {
                    return;
                }

                constexpr double speed = 0.2;

                m_border_color = QColor((1.0 - speed) * m_border_color.red() + speed * 20,
                                        (1.0 - speed) * m_border_color.green() + speed * 20,
                                        (1.0 - speed) * m_border_color.blue() + speed * 20,
                                        (1.0 - speed) * m_border_color.alpha() + speed * 160);
                this->update();
            };

            for (unsigned i = 0; i < 20; ++i)
            {
                QTimer::singleShot(20 * i, update_procedure);
            }

            break;
        }
        case QEvent::MouseButtonPress:
        {
            m_pressed            = true;
            m_border_width       = 20.0;
            m_disallow_animation = false;

            this->update();

            const auto update_procedure = [&]() {
                if (m_disallow_animation)
                {
                    return;
                }

                constexpr double speed = 0.4;

                m_border_width = (1.0 - speed) * m_border_width + speed * 5.0;

                this->update();
            };

            for (unsigned i = 1; i < 10; ++i)
            {
                QTimer::singleShot(20 * i, update_procedure);
            }
            break;
        }
        case QEvent::MouseButtonRelease:
        {
            m_pressed = false;

            break;
        }
        default:
            break;
        }
        return QWidget::event(e);
    }

private:
    PreviewWidget* m_preview_widget;
    bool           m_hovered;
    bool           m_pressed;

    bool m_disallow_animation;

    QColor m_border_color;
    double m_border_width;
};

// -----------------------------------------------------------------------------

PreviewWidget::PreviewWidget(const std::pair<unsigned, unsigned>& widget_id, QWidget* parent)
    : QWidget(parent), m_widget_id(widget_id)
{
    m_widget = new enhancer::EnhancerWidget(enhancer::EnhancerWidget::Policy::AspectFill);
    m_widget->setImage(internal::GetImageSingleton());

    QHBoxLayout* layout = new QHBoxLayout();
    this->setLayout(layout);

    layout->setMargin(0);
    layout->addWidget(m_widget);

    m_overlay = new OverlayWidget(this, m_widget);
}

// -----------------------------------------------------------------------------

const Eigen::VectorXd& PreviewWidget::GetParameters() const { return m_parameters; }

void PreviewWidget::SetParameters(const Eigen::VectorXd& parameters, const bool update_widget)
{
    m_parameters = parameters;

    const Eigen::VectorXd transformed_parameters = internal::TransformParams(parameters);

    m_widget->setParameters(transformed_parameters.data());

    if (update_widget)
    {
        m_widget->update();
    }
}

void PreviewWidget::SetMousePressCallback(const std::function<void(const WidgetId&)>& callback)
{
    m_mouse_press_callback = callback;
}

// ---------------------------------------------------------------------------------------------------------------------

QImage PreviewWidget::RenderImage(const Eigen::VectorXd& parameters)
{
    const Eigen::VectorXd transformed_parameters = internal::TransformParams(parameters);

    auto widget = new enhancer::EnhancerWidget(enhancer::EnhancerWidget::Policy::AspectFill);
    widget->setImage(internal::GetImageSingleton());
    widget->setParameters(transformed_parameters.data());
    widget->setFixedSize(1800, 1200);

    return widget->grabFramebuffer();
}

// -----------------------------------------------------------------------------

void PreviewWidget::resizeEvent(QResizeEvent* event)
{
    m_overlay->resize(event->size());
    event->accept();
}

void PreviewWidget::mousePressEvent(QMouseEvent* event) {}

void PreviewWidget::mouseReleaseEvent(QMouseEvent* event)
{
    if (m_mouse_press_callback)
    {
        m_mouse_press_callback(m_widget_id);
    }
}

// -----------------------------------------------------------------------------

unsigned PreviewWidget::GetNumDimensions() { return enhancer::NUM_PARAMETERS; }

const Eigen::VectorXd& PreviewWidget::GetLowerBound()
{
    static Eigen::VectorXd bound = Eigen::VectorXd::Zero(GetNumDimensions());

    return bound;
}

const Eigen::VectorXd& PreviewWidget::GetUpperBound()
{
    static Eigen::VectorXd bound = Eigen::VectorXd::Ones(GetNumDimensions());

    return bound;
}
