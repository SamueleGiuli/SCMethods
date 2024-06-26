module SCMethods

  USE SCM_COMMON, only: set_mixing, get_mixing, set_error, get_error
  USE SCM_LMIXING, only: linear_mixing
  USE SCM_DIIS, only: diis
  implicit none

  public :: set_mixing, get_mixing, set_error, get_error, linear_mixing, diis
  private
contains
  
end module SCMethods
