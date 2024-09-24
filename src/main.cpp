#ifdef OPENMC_MPI
#include <mpi.h>
#endif
#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/message_passing.h"
#include "openmc/particle_restart.h"
#include "openmc/random_ray/random_ray_simulation.h"
#include "openmc/settings.h"


// Enzyme
#include "enzyme/fprt/fprt.h"
#define FROM 16
#define TO 8

template <typename fty> fty *__enzyme_truncate_mem_func(fty *, int, int);
template <typename fty> fty *__enzyme_truncate_op_func(fty *, int, int, int);


void main_enzyme(int argc, char* argv[], int err) {
  using namespace openmc;

  // start problem based on mode
  switch (settings::run_mode) {
  case RunMode::FIXED_SOURCE:
  case RunMode::EIGENVALUE:
    switch (settings::solver_type) {
    case SolverType::MONTE_CARLO:
      err = openmc_run();
      break;
    case SolverType::RANDOM_RAY:
      openmc_run_random_ray();
      err = 0;
      break;
    }
    break;
  case RunMode::PLOTTING:
    err = openmc_plot_geometry();
    break;
  case RunMode::PARTICLE:
    if (mpi::master)
      run_particle_restart();
    err = 0;
    break;
  case RunMode::VOLUME:
    err = openmc_calculate_volumes();
    break;
  default:
    break;
  }
  if (err)
    fatal_error(openmc_err_msg);

  // Finalize and free up memory
  err = openmc_finalize();
  if (err)
    fatal_error(openmc_err_msg);
}

int main(int argc, char* argv[])
{
  using namespace openmc;
  int err;

  // Initialize run -- when run with MPI, pass communicator
#ifdef OPENMC_MPI
  MPI_Comm world {MPI_COMM_WORLD};
  err = openmc_init(argc, argv, &world);
#else
  err = openmc_init(argc, argv, nullptr);
#endif
  if (err == -1) {
    // This happens for the -h and -v flags
    return 0;
  } else if (err) {
    fatal_error(openmc_err_msg);
  }

  __enzyme_truncate_op_func(main_enzyme, FROM, 0, TO)(argc, argv, err);

    // If MPI is in use and enabled, terminate it
#ifdef OPENMC_MPI
  MPI_Finalize();
#endif
}
