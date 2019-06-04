//
// Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of pylibseq.
//
// pylibseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pylibseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pylibseq.  If not, see <http://www.gnu.org/licenses/>.
//

#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <Sequence/PolyTable.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/SimData.hpp>
#include <Sequence/polySiteVector.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <Sequence/stateCounter.hpp>

namespace py = pybind11;

void init_PolyTable(py::module & m)
{
    py::class_<Sequence::PolyTable>(m, "PolyTable",
                                    "Base class for polymorphism tables")
        .def("data", &Sequence::PolyTable::GetData)
        .def("pos", &Sequence::PolyTable::GetPositions)
        .def("assign",
             [](Sequence::PolyTable& p, const Sequence::polySiteVector& v) {
                 auto rv = p.assign(v.cbegin(), v.cend());
                 if (!rv)
                     throw std::runtime_error("assignment failure");
             },
             "Assign data from a tuple of (position,data), where data "
             "are "
             "encoded as an str.")
        .def("assign",
             [](Sequence::PolyTable& p, const std::vector<double>& pos,
                const std::vector<std::string>& d) {
                 auto rv = p.assign(pos, d);
                 if (!rv)
                     throw std::runtime_error("assignment failure");
             },
             "Assign data from a list of positions and a  list of "
             "strings "
             "representing variable positions.")
        .def("empty", &Sequence::PolyTable::empty,
             "Return whether or not container is empty.")
        .def("numsites", &Sequence::PolyTable::numsites,
             "Return number of variable sites (columns).")
        .def("size", &Sequence::PolyTable::size,
             "Return number of elements (rows).")
        .def("position", &Sequence::PolyTable::position,
             "Return value of the i-th position.", py::arg("i"))
        .def("swap", &Sequence::PolyTable::swap,
             "Swap data with another polymorphism table")
        .def("__str__",
             [](const Sequence::PolyTable& pt) {
                 std::ostringstream o;
                 o << pt;
                 return o.str();
             })
        .def("__getitem__",
             [](const Sequence::PolyTable& pt, const std::size_t i) {
                 if (i >= pt.size())
                     {
                         throw py::index_error("index out of range");
                     }
                 return pt[i];
             },
             "Return i-th row.")
        .def("__len__", [](const Sequence::PolyTable& p) { return p.size(); },
             "Return number of elements (rows).")
        .def("__getstate__", [](const Sequence::PolyTable& p) {
            return py::make_tuple(p.GetPositions(), p.GetData());
        });

    py::class_<Sequence::PolySites, Sequence::PolyTable>(
        m, "PolySites",
        "A polymorphism table for Sequence data.  "
        "A/G/C/T/N/0/1 is the alphabet for analysis.")
        .def(py::init<>())
        .def("__init__",
             [](Sequence::PolySites& p, const Sequence::polySiteVector& v) {
                 new (&p) Sequence::PolySites(v.cbegin(), v.cend());
             })
        .def("__setstate__", [](Sequence::PolySites& p, py::tuple t) {
            if (t.size() != 2)
                {
                    throw std::runtime_error("incorrect tuple length");
                }
            new (&p)
                Sequence::PolySites(t[0].cast<std::vector<double>>(),
                                    t[1].cast<std::vector<std::string>>());
        });

    py::class_<Sequence::SimData, Sequence::PolyTable>(
        m, "SimData",
        "A polymorphism table for "
        "binary data.  0/1 = "
        "ancestral/derived.")
        .def(py::init<std::vector<double>, std::vector<std::string>>())
        .def(py::init<>())
        .def("from_stdin",
             [](Sequence::SimData& d) -> bool {
                 if (std::cin.eof())
                     return false;
                 std::cin >> d >> std::ws;
                 return true;
             },
             "Create an instance of this type by reading data from stdin.  "
             "Useful for reading input from msprime or ms.")
        .def("__init__",
             [](Sequence::SimData& d, const Sequence::polySiteVector& p) {
                 new (&d) Sequence::SimData(p.cbegin(), p.cend());
             })
        .def("__setstate__", [](Sequence::SimData& p, py::tuple t) {
            if (t.size() != 2)
                {
                    throw std::runtime_error("incorrect tuple length");
                }
            new (&p) Sequence::SimData(t[0].cast<std::vector<double>>(),
                                       t[1].cast<std::vector<std::string>>());
        });

// Expose functions. We use macros to avoid tedium
#define MAKE_POLYTABLE_MANIP_FUNCTION(NAME, FXN, TYPE)                        \
    m.def(NAME, &Sequence::FXN<TYPE>, py::arg("polytable"),                   \
          py::arg("skip_ancestral") = false, py::arg("ancestral") = 0,        \
          py::arg("gapchar") = "-");

    MAKE_POLYTABLE_MANIP_FUNCTION("removeMono", removeInvariantPos,
                                  Sequence::SimData);
    MAKE_POLYTABLE_MANIP_FUNCTION("removeGaps", removeGaps, Sequence::SimData);
    MAKE_POLYTABLE_MANIP_FUNCTION("removeAmbiguous", removeAmbiguous,
                                  Sequence::SimData);
    MAKE_POLYTABLE_MANIP_FUNCTION("removeMultiHits", removeMultiHits,
                                  Sequence::SimData);
    MAKE_POLYTABLE_MANIP_FUNCTION("removeMissing", removeMissing,
                                  Sequence::SimData);

    MAKE_POLYTABLE_MANIP_FUNCTION("removeMono", removeInvariantPos,
                                  Sequence::PolySites);
    MAKE_POLYTABLE_MANIP_FUNCTION("removeGaps", removeGaps,
                                  Sequence::PolySites);
    MAKE_POLYTABLE_MANIP_FUNCTION("removeAmbiguous", removeAmbiguous,
                                  Sequence::PolySites);
    MAKE_POLYTABLE_MANIP_FUNCTION("removeMultiHits", removeMultiHits,
                                  Sequence::PolySites);
    MAKE_POLYTABLE_MANIP_FUNCTION("removeMissing", removeMissing,
                                  Sequence::PolySites);

    m.def("isValid",
          [](const Sequence::PolyTable& p) {
              return Sequence::polyTableValid(&p);
          },
          "Return True/False if a polymorphism table is/isn't valid for "
          "further analysis.");

    m.def("removeColumns",
          [](const Sequence::SimData& d,
             std::function<bool(const Sequence::stateCounter&)> f) {
              return Sequence::removeColumns(d, f);
          },
          R"delim(
		  Remove columns from a :class:`libsequence.polytable.SimData` 
		  based on a function taking a :class:`libsequence.summstats.StateCounter`
		  as an argument.

		  :param p: A :class:`libsequence.polytable.SimData`.
		  :param f: A callable taking a :class:`libsequence.summstats.StateCounter` and returning True/False.

		  :rtype: A :class:`libsequence.polytable.SimData`.

		  .. note:: Sites are kept for which the callable returns False
		  )delim");

    m.def(
        "removeColumns",
        [](const Sequence::PolySites& d,
           std::function<bool(const Sequence::stateCounter&)> f,
           const bool skip_ancestral, const unsigned anc, const char gapchar) {
            return Sequence::removeColumns(d, f, skip_ancestral, anc, gapchar);
        },
        R"delim(
		  Remove columns from a :class:`libsequence.polytable.PolySites` 
		  based on a function taking a :class:`libsequence.summstats.StateCounter`
		  as an argument.

		  :param p: A :class:`libsequence.polytable.PolySites`.
		  :param f: A callable taking a :class:`libsequence.summstats.StateCounter` and returning True/False.
		  :param skip_ancestral: (False) Whether to skip ancestral sequence.
		  :param anc: (0) Index of ancestral sequence.
		  :param gapchar: ('-') The gap character.

		  :rtype: A :class:`libsequence.polytable.PolySites`.

		  .. note:: Sites are kept for which the callable returns False
		  )delim",
        py::arg("d"), py::arg("fxn"), py::arg("skip_ancestral") = false,
        py::arg("anc") = 0, py::arg("gapchar") = '-');
}
