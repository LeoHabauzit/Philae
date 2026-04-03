import pylab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from simcoon import simmit as sim
import os
from IPython.display import HTML
from pathlib import Path


basedir = str(Path(__file__).parent)


def umat_smadi(props, typesim):
    dir = os.path.dirname(os.path.realpath("__file__"))
    umat_name = (
        "SMADI"  # This is the 5 character code for the elastic-plastic subroutine
    )
    nstatev = 50  # The number of scalar variables required, only the initial temperature is stored here

    ##local orientation
    psi_rve = 0.0
    theta_rve = 0.0
    phi_rve = 0.0
    solver_type = 0
    corate_type = 3

    path_data = basedir + "/data"
    path_results = basedir + "/results_smadi/"

    # Run the simulation
    pathfile = f"path_{typesim}.txt"
    outputfile = f"results_{typesim}.txt"
    sim.solver(
        umat_name,
        props,
        nstatev,
        psi_rve,
        theta_rve,
        phi_rve,
        solver_type,
        corate_type,
        path_data,
        path_results,
        pathfile,
        outputfile,
    )
