#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_VariantMatrix(py::module &);

PYBIND11_MODULE(_libsequence, m)
{
    init_VariantMatrix(m);
}
