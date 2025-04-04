/*
  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.

  File :      sdo2.cs

  Purpose :   Solves the semidefinite problem with two symmetric variables:

                 min   <C1,X1> + <C2,X2>
                 st.   <A1,X1> + <A2,X2> = b
                             (X2)_{1,2} <= k
                
                 where X1, X2 are symmetric positive semidefinite,

                 C1, C2, A1, A2 are assumed to be constant symmetric matrices,
                 and b, k are constants.
*/
using System;

namespace mosek.example
{
  public class sdo1
  {
    public static void Main(string[] args)
    {
      /* Input data */
      int    numcon      = 2;              /* Number of constraints. */
      int    numbarvar   = 2;
      int[]    dimbarvar = {3, 4};         /* Dimension of semidefinite variables */

      /* Objective coefficients concatenated */
      int[]    Cj = { 0, 0, 1, 1, 1, 1 };   /* Which symmetric variable (j) */
      int[]    Ck = { 0, 2, 0, 1, 1, 2 };   /* Which entry (k,l)->v */
      int[]    Cl = { 0, 2, 0, 0, 1, 2 };
      double[] Cv = { 1.0, 6.0, 1.0, -3.0, 2.0, 1.0 };

      /* Equality constraints coefficients concatenated */
      int[]    Ai = { 0, 0, 0, 0, 0, 0 };   /* Which constraint (i = 0) */
      int[]    Aj = { 0, 0, 0, 1, 1, 1 };   /* Which symmetric variable (j) */
      int[]    Ak = { 0, 2, 2, 1, 1, 3 };   /* Which entry (k,l)->v */
      int[]    Al = { 0, 0, 2, 0, 1, 3 };
      double[] Av = { 1.0, 1.0, 2.0, 1.0, -1.0, -3.0 };

      /* The second constraint - one-term inequality */
      int[]    A2i = { 1 };                        /* Which constraint (i = 1) */
      int[]    A2j = { 1 };                        /* Which symmetric variable (j = 1) */
      int[]    A2k = { 1 };                        /* Which entry A(1,0) = A(0,1) = 0.5 */
      int[]    A2l = { 0 };
      double[] A2v = { 0.5 };

      mosek.boundkey[] bkc = { mosek.boundkey.fx,
                               mosek.boundkey.up
                             };
      double[]     blc     = { 23.0, 0.0 };
      double[]     buc     = { 23.0, -3.0 };

      // Create a task object.
      using (mosek.Task task = new mosek.Task())
      {
        // Directs the log task stream to the user specified
        // method task_msg_obj.stream
        task.set_Stream (mosek.streamtype.log, new msgclass (""));

        /* Append numcon empty constraints.
           The constraints will initially have no bounds. */
        task.appendcons(numcon);

        /* Append numbarvar semidefinite variables. */
        task.appendbarvars(dimbarvar);

        /* Set objective (6 nonzeros).*/
        task.putbarcblocktriplet(Cj, Ck, Cl, Cv);

        /* Set the equality constraint (6 nonzeros).*/
        task.putbarablocktriplet(Ai, Aj, Ak, Al, Av);

        /* Set the inequality constraint (1 nonzero).*/
        task.putbarablocktriplet(A2i, A2j, A2k, A2l, A2v);

        /* Set constraint bounds */
        task.putconboundslice(0, 2, bkc, blc, buc);


        /* Run optimizer */
        task.optimize();
        task.solutionsummary(mosek.streamtype.msg);

        mosek.solsta solsta = task.getsolsta(mosek.soltype.itr);

        switch (solsta) {
          case mosek.solsta.optimal:

            /* Retrieve the soution for all symmetric variables */
            Console.WriteLine("Solution (lower triangular part vectorized):");
            for(int i = 0; i < numbarvar; i++) {
              int dim = dimbarvar[i] * (dimbarvar[i] + 1) / 2;
              double[] barx = new double[dim];

              task.getbarxj(mosek.soltype.itr, i, barx);

              Console.Write("X" + (i+1) + ": ");
              for (int j = 0; j < dim; ++j)
                Console.Write(barx[j] + " ");
              Console.WriteLine();
            }

            break;
          case mosek.solsta.dual_infeas_cer:
          case mosek.solsta.prim_infeas_cer:
            Console.WriteLine("Primal or dual infeasibility certificate found.");
            break;
          case mosek.solsta.unknown:
            Console.WriteLine("The status of the solution could not be determined.");
            break;
          default:
            Console.WriteLine("Other solution status.");
            break;
        }
      }
    }
  }

  class msgclass : mosek.Stream
  {
    string prefix;
    public msgclass (string prfx)
    {
      prefix = prfx;
    }

    public override void streamCB (string msg)
    {
      Console.Write ("{0}{1}", prefix, msg);
    }
  }
}
