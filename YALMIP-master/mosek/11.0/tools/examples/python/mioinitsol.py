##
#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File :      mioinitsol.py
#
#  Purpose :   Demonstrates how to solve a small mixed
#              integer linear optimization problem using the MOSEK Python API.
##
import sys
import mosek

# Since the actual value of Infinity is ignores, we define it solely
# for symbolic purposes:
inf = 0.0

# Define a stream printer to grab output from MOSEK
def streamprinter(text):
    sys.stdout.write(text)
    sys.stdout.flush()


def main():
    # Create a task
    with mosek.Task() as task:
        # Attach a printer to the task
        task.set_Stream(mosek.streamtype.log, streamprinter)

        bkc = [mosek.boundkey.up]
        blc = [-inf, ]
        buc = [2.5]

        bkx = [mosek.boundkey.lo,
               mosek.boundkey.lo,
               mosek.boundkey.lo,
               mosek.boundkey.lo]

        blx = [0.0, 0.0, 0.0, 0.0]
        bux = [inf, inf, inf, inf]

        c = [7.0, 10.0, 1.0, 5.0]

        asub = [0, 0, 0, 0]
        acof = [1.0, 1.0, 1.0, 1.0]

        ptrb = [0, 1, 2, 3]
        ptre = [1, 2, 3, 4]

        numvar = len(bkx)
        numcon = len(bkc)

        # Input linear data
        task.inputdata(numcon, numvar,
                       c, 0.0,
                       ptrb, ptre, asub, acof,
                       bkc, blc, buc,
                       bkx, blx, bux)

        # Input objective sense
        task.putobjsense(mosek.objsense.maximize)

        # Define variables to be integers
        task.putvartypelist([0, 1, 2],
                            [mosek.variabletype.type_int]*3)

        # Assign values to integer variables. 
        # (We only set a slice of xx)
        task.putxxslice(mosek.soltype.itg, 0, 3, [1.0, 1.0, 0.0])

        # Request constructing the solution from integer variable values
        task.putintparam(mosek.iparam.mio_construct_sol, mosek.onoffkey.on)

        # Optimize
        task.optimize()
        task.writedata("mioinitsol.ptf")

        task.solutionsummary(mosek.streamtype.log)

        if task.solutiondef(mosek.soltype.itg):
            # Output a solution
            xx = task.getxx(mosek.soltype.itg)

            print("Integer optimal solution")
            for j in range(0, numvar):
                print("\tx[%d] = %e" % (j, xx[j]))

            # Was the initial guess used?
            constr = task.getintinf(mosek.iinfitem.mio_construct_solution)
            constrVal = task.getdouinf(mosek.dinfitem.mio_construct_solution_obj)
            print("Construct solution utilization: {0}\nConstruct solution objective: {1:.3f}\n".format(constr, constrVal))
        else:
            print("No integer solution is available.")


# call the main function
try:
    main()
except mosek.MosekException as e:
    print("ERROR: %s" % str(e.errno))
    if e.msg is not None:
        print("\t%s" % e.msg)
    sys.exit(1)
except:
    import traceback
    traceback.print_exc()
    sys.exit(1)
