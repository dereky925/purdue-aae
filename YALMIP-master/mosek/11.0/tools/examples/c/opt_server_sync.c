/*
  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.

  File:    opt_server_sync.c

  Purpose :   Demonstrates how to use MOSEK OptServer
              to solve optimization problem synchronously
*/
#include "mosek.h"

static void MSKAPI printstr(void *handle, const char str[])
{
  printf("%s", str);
}

int main(int argc, const char * argv[])
{
  MSKenv_t    env  = NULL;
  MSKtask_t   task = NULL;
  MSKrescodee res  = MSK_RES_OK;
  MSKrescodee trm  = MSK_RES_OK;

  if (argc <= 3)
  {
    fprintf(stderr, "Syntax: opt_server_sync inputfile addr [certpath]\n");    
    return 0;
  }
  else
  {
    const char * taskfile = argv[1];
    const char * address  = argv[2];
    const char * certfile = argc > 3 ? argv[3] : NULL;

    // Create the mosek environment.
    // The `NULL' arguments here, are used to specify customized
    // memory allocators and a memory debug file. These can
    // safely be ignored for now.
    res = MSK_makeenv(&env, NULL);

    // Create a task object linked with the environment env.
    // We create it with 0 variables and 0 constraints initially,
    // since we do not know the size of the problem.
    if (res == MSK_RES_OK)
      res = MSK_maketask(env, 0, 0, &task);

    // Direct the task log stream to a user specified function
    if (res == MSK_RES_OK)
      res = MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);

    // We assume that a problem file was given as the first command
    // line argument (received in `argv')
    if (res == MSK_RES_OK)
      res = MSK_readdata(task, taskfile);

    // Set OptServer URL
    if (res == MSK_RES_OK)
      res = MSK_putoptserverhost(task, address);
    
    // Path to certificate, if any
    if (MSK_RES_OK == res && certfile)
      res = MSK_putstrparam(task, MSK_SPAR_REMOTE_TLS_CERT_PATH,certfile);

    // Solve the problem remotely
    if (res == MSK_RES_OK)
      res = MSK_optimizetrm(task, &trm);
    printf("%s:%d: res = %d\n",__FILE__,__LINE__,res);

    // Print a summary of the solution.
    if (res == MSK_RES_OK)
      res = MSK_solutionsummary(task, MSK_STREAM_LOG);

    // Delete task and environment
    MSK_deletetask(&task);
    MSK_deleteenv(&env);
  }
  return res;
}
