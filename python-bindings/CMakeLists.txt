pybind11_add_module(pysps pysps.cpp)
target_link_libraries(pysps PRIVATE sps pybind11::module)
