from classical_gate_simulations import *
from numpy import eye,kron

def RXwall_Omega(n, Omega, theta):
    circuit = RX(theta*Omega[0])  # Starting with the first element Omega[0]
    for i in range(1, n):  # Start the loop from the second element Omega[1]
        circuit = kron(circuit, RX(theta*Omega[i]))  # Append new gate to the right
    return circuit

def RYwall_Omega(n, Omega, theta):
    circuit = RY(theta*Omega[0])  # Starting with the first element Omega[0]
    for i in range(1, n):  # Start the loop from the second element Omega[1]
        circuit = kron(circuit, RY(theta*Omega[i]))  # Append new gate to the right
    return circuit

def RZwall_Omega(n, Omega, theta):
    circuit = RZ(theta*Omega[0])  # Starting with the first element Omega[0]
    for i in range(1, n):  # Start the loop from the second element Omega[1]
        circuit = kron(circuit, RZ(theta*Omega[i]))  # Append new gate to the right
    return circuit

# “upper” triangular wall of RXX gates
def RXXutwall(n,J,theta):
    if n ==1:
        circuit = I
    else:
        circuit = eye(2**n)
        j = 1
        while j<n:
            k = j+1
            while k<=n:
                operator = distantcouplingXX(n,j,k,theta* J[j-1,k-1])
                circuit = circuit @ operator
                k = k+1
            j = j+1
    return circuit

# “upper” triangular wall of RYY gates
def RYYutwall(n,J,theta):
    if n ==1:
        circuit = I
    else:
        circuit = eye(2**n)
        j = 1
        while j<n:
            k = j+1
            while k<=n:
                operator = distantcouplingYY(n,j,k,theta* J[j-1,k-1])
                circuit = circuit @ operator
                k = k+1
            j = j+1
    return circuit

# “upper” triangular wall of RZZ gates
def RZZutwall(n,J,theta):
    if n ==1:
        circuit = I
    else:
        circuit = eye(2**n)
        j = 1
        while j<n:
            k = j+1
            while k<=n:
                operator = distantcouplingZZ(n,j,k,theta* J[j-1,k-1])
                circuit = circuit @ operator
                k = k+1
            j = j+1
    return circuit


def Xblock(n,h,step_coefficient_in_splitting,paraX,RXXutwall_value): 
    """
    Constructs a quantum circuit block with:
    - a time-dependent single-gate layer and 
    - a pre-calculated, time-independent coupling gate wall.
    Allows flexibility in the splitting weight (step coefficient) for the time-dependent single-gate layer.

    Parameters:
    - paraX: Parameter for the RX gate construction.
    - RXXutwall_value: Predefined time-independent coupling gate wall.

    Returns:
    - circuit: The constructed quantum circuit block.
    """

    # step coefficient 1 (e^Ah) equivalent to RX(2*h), step 1/2 (e^Ah/2) equivalent to RX(h)
    # since Rxx(theta) = exp(-i*theta*X/2)
    circuit = RXwall_Omega(n,paraX,step_coefficient_in_splitting*2*h) @ RXXutwall_value
    return circuit

def Yblock(n,h,step_coefficient_in_splitting,paraY,RYYutwall_value):
    circuit = RYwall_Omega(n,paraY,step_coefficient_in_splitting*2*h) @ RYYutwall_value
    return circuit

def Zblock(n,h,step_coefficient_in_splitting,paraZ,RZZutwall_value):
    circuit = RZwall_Omega(n,paraZ,step_coefficient_in_splitting*2*h) @ RZZutwall_value
    return circuit