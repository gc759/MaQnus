from numpy import array,cos,sin,exp,eye,kron,pi

# Identity gate
I = eye(2)

# Pauli Gates. They are single-qubit gates/operators, ie 2by2 matrices
X = array([[0 ,1],[1, 0]])
Y = array([[0 ,-1j],[1j, 0]])
Z = array([[1 ,0],[0, -1]])

# # Hadmard gate
# Hadmard = array([[1,1],[1,-1]])/sqrt(2)

# NOT gate. This is a single-qubit gate/operator, ie a 2by2 matrix
NOT = X

# Phase shift gate
def P(phi):
    Pphi = array([[1,0],[0,exp(1j*phi)]])
    return Pphi

# Rotation gates. They are single-qubit gates/operators, ie 2by2 matrices
def RX(theta): # exp(-iXθ/2)
    RXtheta = array([[cos(theta/2), -1j * sin(theta/2)], [-1j*sin(theta/2) ,cos(theta/2)]])
    return RXtheta

def RY(theta): # exp(-iYθ/2)
    RYtheta = array([[cos(theta/2), -sin(theta/2)], [sin(theta/2) ,cos(theta/2)]])
    return RYtheta

def RZ(theta): # exp(-iZθ/2)
    RZtheta = array([[exp(-1j*theta/2) ,0],[0 ,exp(1j*theta/2)]])
    return RZtheta

def Single_Gate_Among_nQubits(n,index,gate):
    """
    Applies a specified quantum gate to a single target qubit in a circuit with multiple qubits.

    Parameters:
    - num_qubits (int): Total number of qubits in the circuit.
    - target_qubit (int): Index of the qubit where the gate will be applied.
    - gate (Gate): The quantum gate to apply (e.g., X, H, RX, etc.).

    Returns:
    - circuit: A quantum circuit with the specified gate applied to the target qubit.
    """

    if index<=n:
        i = 1
        op = gate
        while i<index:
            op = kron(op,I)
            i = i+1

        while i<n:
            op = kron(I,op)
            i = i+1

        return op
    else:
        print("invalid input of qubit number or gate index")


def CNOTdistant(n,a,b):
    """
    Constructs a layer with a controlled-NOT (CNOT) gate applied between two distant qubits.

    Parameters:
    - num_qubits (int): Total number of qubits in the circuit.
    - control_qubit (int): Index of the control qubit for the CNOT gate.
    - target_qubit (int): Index of the target qubit where the NOT operation is applied.

    Returns:
    - circuit: The quantum circuit layer with the distant CNOT gate applied.
    """

    if n>=a and n>=b and a!=b:
        if a<b:
            counter1 = 1
            counter2 = 1
            ctrl1 = (I+Z)/2
            ctrl2 = (I-Z)/2

            op1 = ctrl1
            while counter1<a:
                    op1 = kron(op1,I)
                    counter1 = counter1+1

            while counter1<n:
                op1 = kron(I,op1)
                counter1 = counter1+1

            op2 = ctrl2
            while counter2<a:
                    op2 = kron(op2,I)
                    counter2 = counter2+1
    
            counter2 = counter2+1
            while counter2<b:
                op2 = kron(I,op2)
                counter2 = counter2+1

            op2 = kron(NOT,op2)
            counter2 = counter2+1

            while counter2<n+1:
                op2 = kron(I,op2)
                counter2 = counter2+1

            OP = op1+op2

        

        else:
            counter1 = 1
            counter2 = 1
            ctrl1 = (I+Z)/2
            ctrl2 = (I-Z)/2

            op1 = ctrl1
            while counter1<a:
                    op1 = kron(op1,I)
                    counter1 = counter1+1

            while counter1<n:
                op1 = kron(I,op1)
                counter1 = counter1+1


            op2 = NOT
            while counter2<b:
                    op2 = kron(op2,I)
                    counter2 = counter2+1
            counter2 = counter2+1
            while counter2<a:
                op2 = kron(I,op2)
                counter2 = counter2+1

            op2 = kron(ctrl2,op2)

            counter2 = counter2+1

            while counter2<n+1:
                op2 = kron(I,op2)
                counter2 = counter2+1

            OP = op1+op2            

    else:
        print("invalid input of qubit number or gate index")

    return OP

# distant XX coupling gate (XX coupling of site a and b upon n qubits)
def distantcouplingXX(n,a,b,theta):
    l1 = CNOTdistant(n,a,b)
    l2 = Single_Gate_Among_nQubits(n,a,RX(theta))
    l3 = CNOTdistant(n,a,b)
    op = l3 @ l2 @ l1
    return op

# distant YY coupling gate (YY coupling of site a and b upon n qubits)
def distantcouplingYY(n,a,b,theta):
    l1 = Single_Gate_Among_nQubits(n,b,P(pi/2))
    l2 = CNOTdistant(n,a,b)
    l3 = Single_Gate_Among_nQubits(n,a,RY(-theta))
    l4 = CNOTdistant(n,a,b)
    l5 = Single_Gate_Among_nQubits(n,b,P(-pi/2))
    op = l5 @ l4 @ l3 @ l2 @l1
    return op

# distant ZZ coupling gate (ZZ coupling of site a and b upon n qubits)
def distantcouplingZZ(n,a,b,theta):
    l1 = CNOTdistant(n,a,b)
    l2 = Single_Gate_Among_nQubits(n,b,RZ(theta))
    l3 = CNOTdistant(n,a,b)
    op = l3 @ l2 @ l1
    return op