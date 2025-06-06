##
# Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
# File:      qcqp_sdo_relaxation.py
#
# Purpose:   Demonstrate how to use SDP to solve
#            convex relaxation of a mixed-integer QCQP 
##

import math
import sys
import numpy
import mosek
from mosek.fusion import *
import mosek.fusion.pythonic

# The relaxed SDP model
def miqcqp_sdo_relaxation(n,P,q):
    M = Model()

    M.setLogHandler(sys.stdout) 

    Z = M.variable("Z", Domain.inPSDCone(n+1))
    X = Z[:-1, :-1]
    x = Z[:-1, -1:]

    M.constraint( X.diag() >= x )
    M.constraint( Z[n, n] == 1 )

    M.objective( ObjectiveSense.Minimize, 
        Expr.dot( P, X ) + 2.0 * Expr.dot(x, q) )

    return M

# A direct integer model for minimizing |Ax-b|
def int_least_squares(n, A, b):
    M = Model()

    M.setLogHandler(sys.stdout) 

    x = M.variable("x", n, Domain.integral(Domain.unbounded()))
    t = M.variable("t", 1, Domain.unbounded())

    M.constraint( Expr.vstack(t, A @ x - b), Domain.inQCone() )
    M.objective( ObjectiveSense.Minimize, t )

    return M

# problem dimensions
n = 20
m = 2*n

# problem data
A = numpy.reshape(numpy.random.normal(0., 1.0, n*m), (m,n))
c = numpy.random.uniform(0., 1.0, n)
P = A.transpose().dot(A)
q = - P.dot(c)
b = A.dot(c)

# solve the problems
M = miqcqp_sdo_relaxation(n, P, q)
Mint = int_least_squares(n, A, b)
M.solve()
Mint.solve()

M.writeTask('M.ptf')
Mint.writeTask('Mint.ptf')

# rounded and optimal solution
xRound = numpy.rint(M.getVariable("Z").slice([0,n], [n,n+1]).level())
xOpt = numpy.rint(Mint.getVariable("x").level())

print(M.getSolverDoubleInfo("optimizerTime"), Mint.getSolverDoubleInfo("optimizerTime"))
print(numpy.linalg.norm(A.dot(xRound)-b), numpy.linalg.norm(A.dot(xOpt)-b)) 
