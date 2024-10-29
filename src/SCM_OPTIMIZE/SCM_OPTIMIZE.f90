module SCM_OPTIMIZE
  !> Module containing minimization algorithms
  !> The interface require at each call an array of parameters "new_xval"
  !> The module will store the last array of parameters and will perform a linear mixing with the new one
  !> At convergence the logical variable "convg" will be set to .true.
  !> The mixing is set via the set_mixing subroutine
  !> The error tolerance is set via the set_error subroutine
  !> If old_val is present, the routine will return the old value inside the array
  !> the routine is available for rank 1 and rank 2 arrays in single and double precision

  use SCM_COMMON
  use SCM_CG
  use SCM_SD
  !use SCM_PRAXIS
  use iso_fortran_env, only: dp => real64, sp => real32
  implicit none
  
  interface conjugate_gradient
     module procedure fmin_cg_f, fmin_cg_df
  end interface conjugate_gradient

  interface steepest_descent
     module procedure fmin_sd_f
  end interface steepest_descent

  public :: conjugate_gradient, steepest_descent

  private
  

contains
end module SCM_OPTIMIZE
