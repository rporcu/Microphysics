#ifndef ARMONIC_OSCILLATOR_H
#define ARMONIC_OSCILLATOR_H

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <react_util.H>

using namespace unit_test_rp;

void armonic_oscillator()
{
  amrex::Real solution[INT_NEQS] = {1., 0., 0.};

  const amrex::Real tstart(0.);
  const amrex::Real tend(10.0);
  const amrex::Real tinterval(tend - tstart);
  const int nsteps = 10000;

  std::cout << "Maximum Time (s): " << tend << std::endl;
  for (int n = 0; n < INT_NEQS-1; ++n) {
    std::cout << "Elem (" << n << "): " << solution[n] << std::endl;
  }

  burn_t burn_state;

  for (int n = 0; n < INT_NEQS-1; ++n) {
    burn_state.xn[n] = solution[n];
  }
  burn_state.rho = 0;
  burn_state.T = 0;
  burn_state.e = 0;

  // output just the instantaneous RHS
  RArray1D ydot;
  actual_rhs(burn_state, ydot);

  std::cout << "RHS at t = 0" << std::endl;
  for(int n = 0; n < INT_NEQS-1; ++n){
    std::cout << std::setw(6) << "Elem (" << n << "): " << ydot(n+1) << std::endl;
  }

  // loop over steps, burn, and output the current state
  unsigned long nstep_int = 0;
  amrex::Real t(tstart);

  std::ofstream output;
  output.open("solution.csv");
  output << "time,position,velocity,exact_solution,abs_error" << std::endl;

  amrex::Real exact = std::cos(t);
  amrex::Real error = std::abs(exact - burn_state.xn[0]);

  output << t << "," << burn_state.xn[0] << "," << burn_state.xn[1] << "," <<
    exact << "," << error << std::endl;

  for (int n = 0; n < nsteps; n++) {

    // compute the time we wish to integrate to
    amrex::Real dt = tinterval / amrex::Real(nsteps);

    burner(burn_state, dt);

    if (! burn_state.success) {
      amrex::Error("integration failed");
      amrex::Abort("Error");
    }

    nstep_int += burn_state.n_step;

    t += dt;

    exact = std::cos(t);
    error = std::abs(exact - burn_state.xn[0]);

    output << t << "," << burn_state.xn[0] << "," << burn_state.xn[1] << "," <<
      exact << "," << error << std::endl;
  }

  output.close();

  // output diagnostics to the terminal
  std::cout << "------------------------------------" << std::endl;
  std::cout << "successful? " << burn_state.success << std::endl;

  std::cout << "------------------------------------" << std::endl;
  std::cout << "new solution: " << std::endl;
  for (int n = 0; n < INT_NEQS-1; ++n) {
    std::cout << "Elem (" << n << "): " << std::setw(18) << std::setprecision(16) << burn_state.xn[n] << std::endl;
  }

  std::cout << "number of steps taken: " << nstep_int << std::endl;

}
#endif
