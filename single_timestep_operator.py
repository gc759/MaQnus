# from ClassicallyDefineGateLayers import *
# import numpy as np
from circuit_structures import *
from numpy import ones, eye


def trotter1_1step_circuit(n,h,paraX,paraY,paraZ,RXXutwall_value_trotter_c1, RZZutwall_value_trotter_c2, RYYutwall_value_trotter_c3):
    """
    Calculate the single-step time evolution operator using 1st order Trotter splitting.

    Returns:
    - ndarray: The time evolution operator for a single step.
    """

    circuit1 = Xblock(n,h,1,paraX,RXXutwall_value_trotter_c1)
    circuit2 = Zblock(n,h,1,paraZ,RZZutwall_value_trotter_c2)
    circuit3 = Yblock(n,h,1,paraY,RYYutwall_value_trotter_c3)

    circuit = circuit3@circuit2@circuit1
    return circuit

def strang2_1step_circuit(n,h,paraX,paraY,paraZ,RXXutwall_value_strang_c1, RZZutwall_value_strang_c2, RYYutwall_value_strang_c3):

    circuit1 = Xblock(n,h,1/2,paraX,RXXutwall_value_strang_c1)
    circuit2 = Zblock(n,h,1/2,paraZ,RZZutwall_value_strang_c2)
    circuit3 = Yblock(n,h,1,paraY,RYYutwall_value_strang_c3)
    circuit4 = circuit2
    circuit5 = circuit1

    circuit = circuit5@circuit4@circuit3@circuit2@circuit1
    return circuit

def yoshida4_1step_circuit(n,h,paraX,paraY,paraZ,RXXutwall_value_yoshida4_c1, RZZutwall_value_yoshida4_c2, RYYutwall_value_yoshida4_c3, RXXutwall_value_yoshida4_c5, RZZutwall_value_yoshida4_c6, RYYutwall_value_yoshida4_c7):
    
    x0 = -2**(1/3)/(2-2**(1/3))
    x1 = 1/(2-2**(1/3))

    # First part
    circuit1 = Xblock(n,h,x1/2,paraX,RXXutwall_value_yoshida4_c1)
    circuit2 = Zblock(n,h,x1/2,paraZ,RZZutwall_value_yoshida4_c2) #RYwall_Omega(n,paraY,x1*h) @ RYYutwall(n,J,x1*h)
    circuit3 = Yblock(n,h,x1,paraY,RYYutwall_value_yoshida4_c3) #RZwall_Omega(n,paraZ,2*x1*h) @ RZZutwall(n,J,2*x1*h)
    circuit4 = circuit2

    # Combined part
    circuit5 = Xblock(n,h,(x1+x0)/2,paraX,RXXutwall_value_yoshida4_c5) #RXwall_Omega(n,paraX,(x1+x0)*h) @ RXXutwall(n,J,(x1+x0)*h)

    # Second part
    circuit6 = Zblock(n,h,x0/2,paraZ,RZZutwall_value_yoshida4_c6) #RYwall_Omega(n,paraY,x0*h) @ RYYutwall(n,J,x0*h)
    circuit7 = Yblock(n,h,x0,paraY,RYYutwall_value_yoshida4_c7) #RZwall_Omega(n,paraZ,2*x0*h) @ RZZutwall(n,J,2*x0*h)
    circuit8 = circuit6
    
    # Combined part
    circuit9 = circuit5
    
    # Third part
    circuit10 = circuit4
    circuit11 = circuit3
    circuit12 = circuit2
    circuit13 = circuit1

    circuit = circuit13@circuit12@circuit11@circuit10@circuit9@circuit8@circuit7@circuit6@circuit5@circuit4@circuit3@circuit2@circuit1
    return circuit