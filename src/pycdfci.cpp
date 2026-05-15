#include "lib/pybind11/include/pybind11/pybind11.h"
#include "lib/pybind11/include/pybind11/stl.h"
#include "python_api.h"

namespace py = pybind11;

PYBIND11_MODULE(_cdfci, m) {
    m.doc() = "Python bindings for CDFCI";

    py::class_<Result>(m, "CDFCIResult")
        .def(py::init<>())
        .def_readonly("energy", &Result::energy)
        .def_readonly("state_energies", &Result::state_energies)
        .def_readonly("iterations", &Result::iterations);

    py::class_<CDFCIDriverFacade>(m, "CDFCI")
        .def(py::init<const std::string &>(), py::arg("fcidump_path"))
        .def("set_num_iterations", &CDFCIDriverFacade::set_num_iterations, py::arg("n"))
        .def("set_report_interval", &CDFCIDriverFacade::set_report_interval, py::arg("n"))
        .def("set_max_memory", &CDFCIDriverFacade::set_max_memory, py::arg("gb"))
        .def("set_max_load_factor", &CDFCIDriverFacade::set_max_load_factor, py::arg("x"))
        .def("run", &CDFCIDriverFacade::run, py::call_guard<py::gil_scoped_release>()); // GIL release

        m.def("run_cdfci_json",
            &CDFCIRunnerFacade::run_cdfci_json,
            py::arg("option_json"),
            py::call_guard<py::gil_scoped_release>());

        m.def("run_xcdfci_json",
            &CDFCIRunnerFacade::run_xcdfci_json,
            py::arg("option_json"),
            py::call_guard<py::gil_scoped_release>());

        m.def("run_cdfci_rhf_integrals_json",
            &CDFCIRunnerFacade::run_cdfci_rhf_integrals_json,
            py::arg("option_json"),
            py::arg("h1e_flat"),
            py::arg("eri_flat"),
            py::arg("norb"),
            py::arg("nelec"),
            py::arg("core_energy"),
            py::arg("ms2") = 0,
            py::call_guard<py::gil_scoped_release>());

        m.def("run_optorbfci_json",
            &CDFCIRunnerFacade::run_optorbfci_json,
            py::arg("option_json"),
            py::call_guard<py::gil_scoped_release>());

        m.def("tools_frozen_orbital_json",
            &CDFCIRunnerFacade::tools_frozen_orbital_json,
            py::arg("option_json"),
            py::call_guard<py::gil_scoped_release>());

        m.def("tools_symmetry_connection_json",
            &CDFCIRunnerFacade::tools_symmetry_connection_json,
            py::arg("option_json"),
            py::call_guard<py::gil_scoped_release>());
}