/*
  File : portfolio_1_basic.java

  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.

  Purpose :   Implements a basic portfolio optimization model.
*/
package com.mosek.fusion.examples;

import mosek.fusion.*;
import java.io.FileReader;
import java.io.BufferedReader;
import java.util.ArrayList;

public class portfolio_1_basic {
  public static double sum(double[] x) {
    double r = 0.0;
    for (int i = 0; i < x.length; ++i) r += x[i];
    return r;
  }

  public static double dot(double[] x, double[] y) {
    double r = 0.0;
    for (int i = 0; i < x.length; ++i) r += x[i] * y[i];
    return r;
  }

  /*
  Purpose:
      Computes the optimal portfolio for a given risk

  Input:
      n: Number of assets
      mu: An n dimmensional vector of expected returns
      GT: A matrix with n columns so (GT')*GT  = covariance matrix
      x0: Initial holdings
      w: Initial cash holding
      gamma: Maximum risk (=std. dev) accepted

  Output:
      Optimal expected return and the optimal portfolio
  */
  public static double BasicMarkowitz
  ( int n,
    double[] mu,
    double[][] GT,
    double[] x0,
    double   w,
    double   gamma)
  throws mosek.fusion.SolutionError {

    Model M = new Model("Basic Markowitz");
    try {
      // Redirect log output from the solver to stdout for debugging.
      // if uncommented.
      //M.setLogHandler(new java.io.PrintWriter(System.out));

      // Defines the variables (holdings). Shortselling is not allowed.
      Variable x = M.variable("x", n, Domain.greaterThan(0.0));

      //  Maximize expected return
      M.objective("obj", ObjectiveSense.Maximize, Expr.dot(mu, x));

      // The amount invested  must be identical to intial wealth
      M.constraint("budget", Expr.sum(x), Domain.equalsTo(w + sum(x0)));

      // Imposes a bound on the risk
      M.constraint("risk", Expr.vstack(gamma, Expr.mul(GT, x)), Domain.inQCone());

      // Solves the model.
      M.solve();

      // Check if the solution is an optimal point
      SolutionStatus solsta = M.getPrimalSolutionStatus();
      if (solsta != SolutionStatus.Optimal)
      {
          // See https://docs.mosek.com/latest/javafusion/accessing-solution.html about handling solution statuses.
          throw new SolutionError(String.format("Unexpected solution status: %s", solsta.toString()));
      }

      return dot(mu, x.level());
    } finally {
      M.dispose();
    }
  }

  /*
    The example. Reads in data and solves the portfolio models.
   */
  public static void main(String[] argv)
  throws java.io.IOException,
         java.io.FileNotFoundException,
         mosek.fusion.SolutionError {

    int        n      = 8;
    double     w      = 59.0;
    double[]   mu     = {0.07197349, 0.15518171, 0.17535435, 0.0898094 , 0.42895777, 0.39291844, 0.32170722, 0.18378628};
    double[]   x0     = {8.0, 5.0, 3.0, 5.0, 2.0, 9.0, 3.0, 6.0};
    double[]   gammas = {36};
    double[][] GT     = {
        {0.30758, 0.12146, 0.11341, 0.11327, 0.17625, 0.11973, 0.10435, 0.10638},
        {0.     , 0.25042, 0.09946, 0.09164, 0.06692, 0.08706, 0.09173, 0.08506},
        {0.     , 0.     , 0.19914, 0.05867, 0.06453, 0.07367, 0.06468, 0.01914},
        {0.     , 0.     , 0.     , 0.20876, 0.04933, 0.03651, 0.09381, 0.07742},
        {0.     , 0.     , 0.     , 0.     , 0.36096, 0.12574, 0.10157, 0.0571 },
        {0.     , 0.     , 0.     , 0.     , 0.     , 0.21552, 0.05663, 0.06187},
        {0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.22514, 0.03327},
        {0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.2202 }
    };

    System.out.println("\n-----------------------------------------------------------------------------------");
    System.out.println("Basic Markowitz portfolio optimization");
    System.out.println("-----------------------------------------------------------------------------------\n");
    for ( int i = 0; i < gammas.length; ++i) {
      double expret = BasicMarkowitz( n, mu, GT, x0, w, gammas[i]);
      System.out.format("Expected return: %.4e Std. deviation: %.4e\n",
                        expret,
                        gammas[i]);
    }
  }
}

