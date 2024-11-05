from numpy import zeros, eye
from full_time_evolution_operator import *
from parameter_generator import *
from error import *


def error_vs_number_of_timestep(n, Omega, J, f_function, g_function, starttime, vectorofnumberoftimestep, terminatetime, exactU, method=None):
    errorvector = []

    if method == "trotter1_startingpoint":
        order = "1st Order Trotter"
    elif method == "strang2_midpoint":
        order = "2nd Order Strang"
    elif method in ["magnus2_scipy", "magnus2_midpoint"]:
        order = "2nd Order Strang"
    elif method in ["magnus4_scipy", "magnus4_GL2", "magnus4_GL3"]:
        order = "4th Order Yoshida"
    elif method == "CF42":
        order = "CF42"
    else:
        raise ValueError("Invalid method specified")

    for num_steps in vectorofnumberoftimestep:
        # Set up parameters based on the chosen method
        parameterstorer = parametergenerator(n, starttime, terminatetime, num_steps, Omega, f_function, g_function, integration_method=method)

        # Compute the full time evolution operator
        Circuit = full_time_evolution_operator(n, order, starttime, terminatetime, num_steps, parameterstorer, J)
        
        # Calculate and store the error for this timestep count
        error = l2error(Circuit, exactU)
        errorvector.append(error)

    return np.array(errorvector)



