// This library needs to be built with C++17. The following three lines are an ad-hoc solution to ignore register
// specifiers in the dependency headers, which are deprecated since C++17.
#if (__cplusplus - 0) >= 201703L
#define register
#endif

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sps/sps.hpp>

namespace py = pybind11;
using namespace py::literals;

void SetSeed(const unsigned seed) { srand(seed); }

PYBIND11_MODULE(pysps, m)
{
    m.def("set_seed", &SetSeed, py::arg("seed"));

    // sps::AbstractPlane

    py::class_<sps::AbstractPlane, std::shared_ptr<sps::AbstractPlane>> abstract_plane_class(m, "AbstractPlane");

    abstract_plane_class.def("get_center", &sps::AbstractPlane::GetCenter);

    abstract_plane_class.def("get_vertices", &sps::AbstractPlane::GetVertices);

    abstract_plane_class.def("calc_parameters", &sps::AbstractPlane::CalcParameters, "coords"_a);

    abstract_plane_class.def("calc_grid_parameters",
                             &sps::AbstractPlane::CalcGridParameters,
                             "grid_cell"_a,
                             "num_candidates"_a,
                             "inter_level_scale"_a,
                             "prev_grid_cells"_a);

    // sps::Optimizer

    py::class_<sps::Optimizer> optimizer_class(m, "Optimizer");

    optimizer_class.def(py::init<const unsigned, const bool>(), "num_dimensions"_a, "use_map_hyperparams"_a = true);

    optimizer_class.def("submit_data", &sps::Optimizer::SubmitData, "x_preferred"_a, "x_others"_a);

    optimizer_class.def("retrieve_search_plane", &sps::Optimizer::RetrieveSearchPlane);
}
