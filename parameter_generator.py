import numpy as np
from scipy import integrate
from cmath import sqrt


def integral1_scipy(h,t,function):
    integrand = lambda zeta: function(t+zeta)
    integral,_ = integrate.quad(integrand, 0, h)
    return (integral)

def integral2_scipy(h,t,function):
    integrand = lambda zeta: (2*zeta-h)*function(t+zeta)
    integral,_ = integrate.quad(integrand, 0, h)
    return (integral)

def integral3_scipy(h,t,functionf,functiong):
    integrand = lambda xi,zeta: functionf(t+xi)*functiong(t+zeta) - functionf(t+zeta)*functiong(t+xi)
    integral,_ = integrate.dblquad(integrand, 0, h  , lambda variable: 0, lambda variable: variable)
    return (integral)

# midpoint
def integral1_midpoint(h,t,function):
    return h*function(t+h/2)

# two Gauss–Legendre knots

def integral1_GL2(h,t,function):

    t1 = h/2*(1+1/sqrt(3))
    t2 = h/2*(1-1/sqrt(3))

    result = h/2*function(t+t1) + h/2*function(t+t2)

    return result

def integral2_GL2(h,t,function):

    t1 = h/2*(1+1/sqrt(3))
    t2 = h/2*(1-1/sqrt(3))

    result = h/2*(2*t1-h)*function(t+t1) + h/2*(2*t2-h)*function(t+t2)

    return result

def integral3_GL2(h,t,functionf,functiong):
    
    knots = np.array([h/2*(1-1/sqrt(3)),h/2*(1+1/sqrt(3))])
    wmatrix = (h**2)*np.array([[1/8,1/8 + 1/(4*sqrt(3))],[1/8 - 1/(4*sqrt(3)),1/8]])
    def L(xi,zeta):
        return functionf(t+xi)*functiong(t+zeta)-functionf(t+zeta)*functiong(t+xi)
    
    result = 0

    for i in range(0,2):
        for j in range (0,2):
            result = result + wmatrix[j,i] * L(knots[j], knots[i])

    return result

def integral1_GL3(h,t,function):

    t1 = h/2*(1-sqrt(3/5))
    t2 = h/2*(1)
    t3 = h/2*(1+sqrt(3/5))

    result = h/2*5/9*function(t+t1) + h/2*8/9*function(t+t2) + h/2*5/9*function(t+t3)

    return result

def integral2_GL3(h,t,function):

    t1 = h/2*(1-sqrt(3/5))
    t2 = h/2*(1)
    t3 = h/2*(1+sqrt(3/5))

    result = h/2*5/9*(2*t1-h)*function(t+t1) + h/2*8/9*(2*t2-h)*function(t+t2) + h/2*5/9*(2*t3-h)*function(t+t3)

    return result

def integral3_GL3(h,t,functionf,functiong):

    knots = np.array([h/2*(1-sqrt(3/5)),h/2*(1),h/2*(1+sqrt(3/5))])
    wmatrix = h**2*np.array([[25/648,1/162*(10+3*sqrt(15)),1/648*(25+6*sqrt(15))],[-1/162*(-10+3*sqrt(15)),8/81,1/162*(10+3*sqrt(15))],[1/648*(25-6*sqrt(15)),-1/162*(-10+3*sqrt(15)),25/648]])

    def L(xi,zeta):
        return functionf(t+xi)*functiong(t+zeta)-functionf(t+zeta)*functiong(t+xi)
    
    result = 0

    for i in range(0,3):
        for j in range (0,3):
            result = result + wmatrix[j,i] * L(knots[j], knots[i])

    return result


def integralA1_scipy(h,t,function):
        integrand = lambda zeta: function(t + zeta)
        integral, _ = integrate.quad(integrand, 0, h)
        return integral
    

def integralA2_scipy(h,t,function):
        integrand = lambda zeta: function(t + zeta) * (2 * (zeta) / h - 1)
        integral, _ = integrate.quad(integrand, 0, h)
        return 3 * integral


def parametergenerator(n, starttime, terminatetime, numberoftimestep, Omega, f_function, g_function, integration_method=None):
    timestep = (terminatetime - starttime) / numberoftimestep
    timegrid = np.linspace(starttime, terminatetime, numberoftimestep + 1)


    if integration_method == "CF42":
        parameterstorer = np.ones((numberoftimestep, 6, n))
        time_offset = 0
        integral1_func = integralA1_scipy
        integral2_func = integralA2_scipy
    else:
        parameterstorer = np.ones((numberoftimestep, 3, n))
        if integration_method == "trotter1_startingpoint":
            time_offset = 0
            integral1_func = None
            integral2_func = None
            integral3_func = None
        elif integration_method == "strang2_midpoint":
            time_offset = timestep / 2
            integral1_func = None
            integral2_func = None
            integral3_func = None
        elif integration_method == "magnus2_scipy":
            time_offset = 0
            integral1_func = integral1_scipy
            integral2_func = None
            integral3_func = None
        elif integration_method == "magnus2_midpoint":
            time_offset = 0
            integral1_func = integral1_midpoint
            integral2_func = None
            integral3_func = None
        elif integration_method == "magnus4_scipy":
            time_offset = 0
            integral1_func = integral1_scipy
            integral2_func = integral2_scipy
            integral3_func = integral3_scipy
        elif integration_method == "magnus4_GL2":
            time_offset = 0
            integral1_func = integral1_GL2
            integral2_func = integral2_GL2
            integral3_func = integral3_GL2
        elif integration_method == "magnus4_GL3":
            time_offset = 0
            integral1_func = integral1_GL3
            integral2_func = integral2_GL3
            integral3_func = integral3_GL3
        else:
            raise ValueError("Invalid integration_method specified")

    # Generate parameters based on the integration method
    for i in range(numberoftimestep):
        t = timegrid[i] + time_offset

        if integration_method == "CF42":
            # CF42-specific parameter calculation
            A1f = 1/2*integralA1_scipy(timestep,t,f_function)/timestep
            A1g = 1/2*integralA1_scipy(timestep,t,g_function)/timestep
            A1O = 1/2*Omega
            A2f = 1/3*integralA2_scipy(timestep,t,f_function)/timestep
            A2g = 1/3*integralA2_scipy(timestep,t,g_function)/timestep
            A2O = 0

            paraX1 = A1f + A2f
            paraY1 = A1g + A2g
            paraZ1 = A1O + A2O

            paraX2 = A1f - A2f
            paraY2 = A1g - A2g
            paraZ2 = A1O - A2O

            parameterstorer[i, 0, :] = paraX1
            parameterstorer[i, 1, :] = paraY1
            parameterstorer[i, 2, :] = paraZ1
            parameterstorer[i, 3, :] = paraX2
            parameterstorer[i, 4, :] = paraY2
            parameterstorer[i, 5, :] = paraZ2

        elif integration_method in ["trotter1_startingpoint", "strang2_midpoint"]:
            paraX = f_function(t) * np.ones(n)
            paraY = g_function(t) * np.ones(n)
            paraZ = Omega

            parameterstorer[i, 0, :] = paraX
            parameterstorer[i, 1, :] = paraY
            parameterstorer[i, 2, :] = paraZ
        else:
            paraX = integral1_func(timestep, t, f_function) / timestep
            paraY = integral1_func(timestep, t, g_function) / timestep
            paraZ = Omega

            # Adjust for fourth-order Magnus methods
            if integration_method.startswith("magnus4"):
                paraX += Omega * (integral2_func(timestep, t, g_function) / timestep)
                paraY -= Omega * (integral2_func(timestep, t, f_function) / timestep)
                paraZ = Omega - (integral3_func(timestep, t, f_function, g_function) / timestep)

            parameterstorer[i, 0, :] = paraX
            parameterstorer[i, 1, :] = paraY
            parameterstorer[i, 2, :] = paraZ

    return parameterstorer