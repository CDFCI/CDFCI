"""Quick smoke tests for the Python CDFCI API.

Run from repository root after building the extension:
    make python-module
    PYTHONPATH=python python -m cdfci.test
"""

import cdfci


def test_legacy_facade():
    drv = cdfci.CDFCI("data/h2o_sto3g_psi4.FCIDUMP")
    drv.set_num_iterations(150000)
    drv.set_report_interval(10000)

    res = drv.run()
    ref = -75.7161
    assert isinstance(res.energy, float)
    assert abs(res.energy - ref) < 1e-3


def test_json_option_api():
    opt = {
        "hamiltonian": {
            "type": "molecule",
            "molecule": {
                "fcidump_path": "data/h2o_sto3g_psi4.FCIDUMP",
                "threshold": 0.0,
            },
        },
        "solver": {
            "type": "cdfci",
            "cdfci": {
                "num_iterations": 150000,
                "report_interval": 10000,
                "z_threshold": 0.0,
                "z_threshold_search": False,
            },
        },
        "max_memory": 0.5,
    }

    res = cdfci.run_cdfci(opt)
    ref = -75.7161
    assert isinstance(res.energy, float)
    assert abs(res.energy - ref) < 1e-3


def test_xcdfci_api_returns_states():
    opt = {
        "hamiltonian": {
            "type": "molecule",
            "molecule": {
                "fcidump_path": "data/h2o_sto3g_psi4.FCIDUMP",
                "threshold": 0.0,
            },
        },
        "solver": {
            "type": "xcdfci",
            "xcdfci": {
                "num_iterations": 150000,
                "report_interval": 10000,
                "z_threshold": 0.0,
                "z_threshold_search": False,
            },
        },
        "num_states": 2,
        "max_memory": 0.5,
    }

    res = cdfci.run_xcdfci(opt)
    assert isinstance(res.energy, float)
    assert isinstance(res.state_energies, list)
    assert len(res.state_energies) >= 2


def main():
    test_legacy_facade()
    test_json_option_api()
    test_xcdfci_api_returns_states()
    print("All Python API smoke tests passed.")


if __name__ == "__main__":
    main()
