#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_PolyTable(py::module & m);
void init_fst(py::module & m);
void init_VariantMatrix(py::module &);
void init_summstats(py::module & );
void init_windows(py::module & );

PYBIND11_MODULE(_libsequence, m)
{
    init_PolyTable(m);
    init_fst(m);
    init_VariantMatrix(m);
    init_summstats(m);
    init_windows(m);
}
