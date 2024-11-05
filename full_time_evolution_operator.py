from single_timestep_operator import *

def full_time_evolution_operator(n,splitting_type,starttime,terminatetime,numberoftimestep,parameterstorer,J):
    """
    Calculate the time evolution operator for an entire simulation from start_time to end_time 
    based on a specified splitting method.

    Parameters:
    - splitting_type: Type of splitting.
    - parameter_storer (ndarray): Parameters for X, Y, and Z evolution for each timestep.

    Returns:
    - ndarray: The resulting time evolution operator.
    """

    h = (terminatetime-starttime)/numberoftimestep

    if splitting_type == '1st Order Trotter':

        # Define Trotterized values with fixed splitting weight of 1
        RXX_time_indep_wall1 = RXXutwall(n,J,1*2*h)
        RZZ_time_indep_wall2 = RZZutwall(n,J,1*2*h)
        RYY_time_indep_wall3 = RYYutwall(n,J,1*2*h)
        
        op = eye(2**n)
        for i in range (numberoftimestep):
            paraX  = parameterstorer[i,0,:]
            paraY  = parameterstorer[i,1,:]
            paraZ  = parameterstorer[i,2,:]

            op =  trotter1_1step_circuit(n,h,paraX,paraY,paraZ,RXX_time_indep_wall1, RZZ_time_indep_wall2, RYY_time_indep_wall3)@op
    
    elif splitting_type == '2nd Order Strang':
        RXX_time_indep_wall1 = RXXutwall(n,J,1/2*2*h)
        RZZ_time_indep_wall2 = RZZutwall(n,J,1/2*2*h)
        RYY_time_indep_wall3 = RYYutwall(n,J,1*2*h)

        op = eye(2**n)
        for i in range (numberoftimestep):
            paraX  = parameterstorer[i,0,:]
            paraY  = parameterstorer[i,1,:]
            paraZ  = parameterstorer[i,2,:]

            op =  strang2_1step_circuit(n,h,paraX,paraY,paraZ,RXX_time_indep_wall1, RZZ_time_indep_wall2, RYY_time_indep_wall3)@op

    elif splitting_type == '4th Order Yoshida':
        x0 = -2**(1/3)/(2-2**(1/3))
        x1 = 1/(2-2**(1/3))

        RXX_time_indep_wall1 = RXXutwall(n,J,x1/2*2*h)
        RZZ_time_indep_wall2 = RZZutwall(n,J,x1/2*2*h)
        RYY_time_indep_wall3 = RYYutwall(n,J,x1*2*h)
        RXX_time_indep_wall5 = RXXutwall(n,J,(x1+x0)/2*2*h)
        RZZ_time_indep_wall6 = RZZutwall(n,J,x0/2*2*h)
        RYY_time_indep_wall7 = RYYutwall(n,J,x0*2*h)
    
        op = eye(2**n)
        for i in range (numberoftimestep):
            paraX  = parameterstorer[i,0,:]
            paraY  = parameterstorer[i,1,:]
            paraZ  = parameterstorer[i,2,:]

            op =  yoshida4_1step_circuit(n,h,paraX,paraY,paraZ,RXX_time_indep_wall1, RZZ_time_indep_wall2, RYY_time_indep_wall3, RXX_time_indep_wall5, RZZ_time_indep_wall6, RYY_time_indep_wall7)@op

    elif splitting_type == "CF42":
        x0 = -2**(1/3)/(2-2**(1/3))
        x1 = 1/(2-2**(1/3))

        RXX_time_indep_wall1 = RXXutwall(n,1/2*J,x1/2*2*h)
        RZZ_time_indep_wall2 = RZZutwall(n,1/2*J,x1/2*2*h)
        RYY_time_indep_wall3 = RYYutwall(n,1/2*J,x1*2*h)
        RXX_time_indep_wall5 = RXXutwall(n,1/2*J,(x1+x0)/2*2*h)
        RZZ_time_indep_wall6 = RZZutwall(n,1/2*J,x0/2*2*h)
        RYY_time_indep_wall7 = RYYutwall(n,1/2*J,x0*2*h)
        
        op = eye(2**n)
        for i in range (numberoftimestep):
            paraX1  = parameterstorer[i,0,:]
            paraY1  = parameterstorer[i,1,:]
            paraZ1  = parameterstorer[i,2,:]
            paraX2  = parameterstorer[i,3,:]
            paraY2  = parameterstorer[i,4,:]
            paraZ2  = parameterstorer[i,5,:]

            op = yoshida4_1step_circuit(n,h,paraX1,paraY1,paraZ1,RXX_time_indep_wall1, RZZ_time_indep_wall2, RYY_time_indep_wall3, RXX_time_indep_wall5, RZZ_time_indep_wall6, RYY_time_indep_wall7)@ yoshida4_1step_circuit(n,h,paraX2,paraY2,paraZ2,RXX_time_indep_wall1, RZZ_time_indep_wall2, RYY_time_indep_wall3, RXX_time_indep_wall5, RZZ_time_indep_wall6, RYY_time_indep_wall7) @ op

    else:
        print("invalid splitting type")

    return op