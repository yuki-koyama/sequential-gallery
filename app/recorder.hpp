#ifndef RECORDER_HPP
#define RECORDER_HPP

#include <Eigen/Core>
#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace sps
{
    class AbstractPlane;
}

class Recorder
{
public:
    Recorder(const std::string& task_type, const std::string& dir_path);

    void ResetClock();
    void AddPlaneSearchRecord(const std::shared_ptr<const sps::AbstractPlane> plane, const Eigen::VectorXd& x_chosen);

    unsigned GetNumRecords() const { return m_records.size(); }

    void ExportCsv() const;
    void ExportImages() const;

private:
    struct PlaneSearchRecord
    {
        PlaneSearchRecord() : t_created(std::chrono::system_clock::now()) {}

        std::shared_ptr<const sps::AbstractPlane> plane;
        Eigen::VectorXd                           x_chosen;
        std::chrono::system_clock::time_point     t_created;
    };

    void PrintHeaderRow(std::ostream* stream) const;
    void PrintRow(const int row, std::ostream* stream) const;

    std::chrono::system_clock::time_point m_base_time_point;

    const std::string m_task_type;
    const std::string m_dir_path;

    std::vector<PlaneSearchRecord> m_records;
};

#endif // RECORDER_HPP
