from numpy import real,trace,conjugate,transpose
from numpy.linalg import norm
from cmath import sqrt


# trace error (overlap)
# this error is only valid for unitary matrices
def traceerror(Uapprox,exactU):
    error = real(sqrt(1-real(trace(conjugate(transpose(Uapprox))@exactU))/Uapprox.shape[0]))
    return error


# l2 error
def l2error(Uapprox,exactU):
    error = norm(exactU-Uapprox)
    return error

