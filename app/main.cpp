#include "app-window.hpp"
#include <QApplication>
#include <QOpenGLFunctions_3_2_Core>

int main(int argc, char* argv[])
{
    QApplication app(argc, argv);

#if !defined(ZIN_APP_FOR_SMPL)
    Q_INIT_RESOURCE(enhancer_resources);
#endif

    constexpr int num_samples = 4;

    QSurfaceFormat format;
    format.setVersion(3, 2);
    format.setProfile(QSurfaceFormat::CoreProfile);
    format.setSamples(num_samples);
    QSurfaceFormat::setDefaultFormat(format);

    AppWindow app_window;
    app_window.show();

    return app.exec();
}
