####
##  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
##
##  File:      lo1.py
##
##  Purpose: Demonstrates how to solve the problem
##
##  maximize 3*x0 + 1*x1 + 5*x2 + x3
##  such that
##           3*x0 + 1*x1 + 2*x2        = 30,
##           2*x0 + 1*x1 + 3*x2 + 1*x3 > 15,
##                  2*x1 +      + 3*x3 < 25
##  and
##           x0,x1,x2,x3 > 0,
##           0 < x1 < 10
####

import sys
from mosek.fusion import *
import mosek.fusion.pythonic

def main(args):
    A = [[3.0, 1.0, 2.0, 0.0],
         [2.0, 1.0, 3.0, 1.0],
         [0.0, 2.0, 0.0, 3.0]]
    c = [3.0, 1.0, 5.0, 1.0]

    # Create a model with the name 'lo1'
    with Model("lo1") as M:

        # Create variable 'x' of length 4
        x = M.variable("x", 4, Domain.greaterThan(0.0))

        # Create constraints
        M.constraint(x[1] <= 10)
        M.constraint("c1", x.T @ A[0] == 30.0)
        M.constraint("c2", x.T @ A[1] >= 15.0)
        M.constraint("c3", x.T @ A[2] <= 25.0)

        # Set the objective function to (c^t * x)
        M.objective("obj", ObjectiveSense.Maximize, x.T @ c)

        # Solve the problem
        M.solve()

        # Get the solution values
        sol = x.level()
        print('\n'.join(["x[%d] = %f" % (i, sol[i]) for i in range(4)]))

if __name__ == '__main__':
    main(sys.argv[1:])
