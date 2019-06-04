#include <pybind11/pybind11.h>
#include <Sequence/SimData.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/PolyTableSlice.hpp>

namespace py = pybind11;

void init_windows(py::module & m)
{
    using SimDataWindows = Sequence::PolyTableSlice<Sequence::SimData>;
    using PolySitesWindows = Sequence::PolyTableSlice<Sequence::PolySites>;

#define MAKE_WINDOWS_BACKEND(CPPNAME, PYNAME)                                 \
    py::class_<CPPNAME>(m, PYNAME)                                            \
        .def("__init__",                                                      \
             [](CPPNAME& sdw, const Sequence::PolyTable& p,                   \
                const double window_size, const double step_len,              \
                const double starting_pos, const double ending_pos) {         \
                 new (&sdw) CPPNAME(p.sbegin(), p.send(), window_size,        \
                                    step_len, starting_pos, ending_pos);      \
             })                                                               \
        .def("__getitem__",                                                   \
             [](const CPPNAME& sdw, const unsigned i) { return sdw[i]; })     \
        .def("__len__", [](const CPPNAME& swd) { return swd.size(); });

    MAKE_WINDOWS_BACKEND(SimDataWindows, "SimDataWindows")
    MAKE_WINDOWS_BACKEND(PolySitesWindows, "PolySitesWindows")
}
