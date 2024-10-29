module SCMethods

  USE SCM_COMMON, only: set_mixing, get_mixing, set_error, get_error
  USE SCM_LMIXING, only: linear_mixing
  USE SCM_DIIS, only: diis
  USE SCM_OPTIMIZE, only: conjugate_gradient, steepest_descent
  implicit none

  public :: set_mixing, get_mixing, set_error, get_error, linear_mixing, diis
  public :: conjugate_gradient, steepest_descent
  private
contains
  
end module SCMethods
