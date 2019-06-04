#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <Sequence/FST.hpp>
#include <Sequence/PolyTable.hpp>
#include <stdexcept>

namespace py = pybind11;

void init_fst(py::module & m)
{
    //py::object polytable
    //    = (py::object)py::module::import("libsequence.polytable")
    //          .attr("PolyTable");

    py::class_<Sequence::FST>(m, "Fst", "Fst")
        .def("__init__",
             [](Sequence::FST& fst, const Sequence::PolyTable& pt,
                const std::vector<unsigned>& config, const bool haveOutgroup,
                const unsigned outgroup) {
                 try
                     {
                         new (&fst)
                             Sequence::FST(&pt, config.size(), &config[0],
                                           nullptr, haveOutgroup, outgroup);
                     }
                 catch (std::runtime_error&)
                     {
                         throw std::runtime_error(
                             "exception caught from Fst constructor");
                     }
             },
             py::arg("polytable"), py::arg("config"),
             py::arg("haveOutgroup") = false, py::arg("outgroup") = 0)
        .def("__init__",
             [](Sequence::FST& fst, const Sequence::PolyTable& pt,
                const std::vector<unsigned>& config,
                const std::vector<double>& weights, const bool haveOutgroup,
                const unsigned outgroup) {
                 try
                     {
                         new (&fst) Sequence::FST(&pt, config.size(),
                                                  &config[0], &weights[0],
                                                  haveOutgroup, outgroup);
                     }
                 catch (std::runtime_error&)
                     {
                         throw std::runtime_error(
                             "exception caught from Fst constructor");
                     }
             },
             py::arg("polytable"), py::arg("config"), py::arg("weights"),
             py::arg("haveOutgroup") = false, py::arg("outgroup") = 0)
        .def("hsm", &Sequence::FST::HSM,
             "Hudson, Slatkin, Maddison (1992) Genetics")
        .def("slatkin", &Sequence::FST::Slatkin, "Slatkin (1991) Genetics")
        .def("hbk", &Sequence::FST::HBK, "Hudson, Boos, Kaplan (1992) MBE")
        .def("piB", &Sequence::FST::piB,
             ":math:`\\pi_B= \\frac{\\sum_{i<j}w_i "
             "w_j \\pi_{ij}}{\\sum_{i<j}w_i w_j}`")
        .def("piS", &Sequence::FST::piS,
             ":math:`\\pi_S = \\frac{\\sum_i w_i^2 \\pi_{ii}}{\\sum_i w_i^2}`")
        .def("piT", &Sequence::FST::piT,
             ":math:`\\pi_T = \\sum_i w_i^2 "
             "\\pi_{ii} + 2\\sum_{i<j}w_i w_j "
             "\\pi_{ij}`")
        .def("piD", &Sequence::FST::piD,
             ":math:`\\pi_D = \\frac{\\pi_T - \\pi_S}{2 \\sum_{i<j}w_i w_j}`")
        .def("shared", &Sequence::FST::shared)
        .def("fixed", &Sequence::FST::fixed)
        .def("priv", &Sequence::FST::Private);
}
