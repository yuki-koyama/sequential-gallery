#include "recorder.hpp"
#include "preview-widget.hpp"
#include <fstream>
#include <sps/plane.hpp>

Recorder::Recorder(const std::string& task_type, const std::string& dir_path)
    : m_task_type(task_type), m_dir_path(dir_path)
{
    this->ResetClock();
}

void Recorder::ResetClock() { m_base_time_point = std::chrono::system_clock::now(); }

void Recorder::AddPlaneSearchRecord(const std::shared_ptr<const sps::AbstractPlane> plane,
                                    const Eigen::VectorXd&                          x_chosen)
{
    auto record = PlaneSearchRecord{};

    record.plane    = plane;
    record.x_chosen = x_chosen;

    m_records.push_back(record);
}

void Recorder::ExportCsv() const
{
    std::ofstream file_stream(m_dir_path + "/record.csv");

    PrintHeaderRow(&file_stream);
    for (int row = 0; row < m_records.size(); ++row)
    {
        PrintRow(row, &file_stream);
    }
}

void Recorder::ExportImages() const
{
    for (int i = 0; i < m_records.size(); ++i)
    {
        const auto& record = m_records[i];

        const QImage center_image = PreviewWidget::RenderImage(record.plane->GetCenter());
        const QImage chosen_image = PreviewWidget::RenderImage(record.x_chosen);

        center_image.save(QString::fromStdString(m_dir_path + "/center_" + std::to_string(i) + ".png"));
        chosen_image.save(QString::fromStdString(m_dir_path + "/chosen_" + std::to_string(i) + ".png"));

#if defined(ZIN_APP_FOR_SMPL)
        PreviewWidget::ExportObj(record.plane->GetCenter(), m_dir_path + "/shape_center_" + std::to_string(i) + ".obj");
        PreviewWidget::ExportObj(record.x_chosen, m_dir_path + "/shape_chosen_" + std::to_string(i) + ".obj");
#endif
    }
}

void Recorder::PrintHeaderRow(std::ostream* stream) const
{
    (*stream) << "index";
    (*stream) << ",file name";
    (*stream) << ",elapsed duration [ms]";
    (*stream) << ",x_chosen";
    (*stream) << ",x_center";
    (*stream) << ",x_0";
    (*stream) << ",x_1";
    (*stream) << ",x_2";
    (*stream) << ",x_3";
    (*stream) << std::endl;
}

void Recorder::PrintRow(const int row, std::ostream* stream) const
{
    assert(row >= 0 && row < m_records.size());

    const auto& record = m_records[row];

    const auto elapsed_duration = record.t_created - m_base_time_point;
    const auto format           = Eigen::IOFormat{4, Eigen::DontAlignCols, ", ", ", ", "", "", "", ""};

    (*stream) << std::to_string(row);
    (*stream) << "," << m_task_type;
    (*stream) << "," << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_duration).count();
    (*stream) << ",\"" << record.x_chosen.format(format) << "\"";
    (*stream) << ",\"" << record.plane->GetCenter().format(format) << "\"";
    (*stream) << ",\"" << record.plane->GetVertices()[0].format(format) << "\"";
    (*stream) << ",\"" << record.plane->GetVertices()[1].format(format) << "\"";
    (*stream) << ",\"" << record.plane->GetVertices()[2].format(format) << "\"";
    (*stream) << ",\"" << record.plane->GetVertices()[3].format(format) << "\"";
    (*stream) << std::endl;
}
