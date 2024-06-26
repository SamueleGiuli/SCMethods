module SCM_COMMON
  use iso_fortran_env, only: dp => real64, sp => real32
  implicit none
  
  public :: set_mixing, get_mixing, set_error, get_error

  real(dp) :: mixing=1.d0, error=1.d-4


  interface set_mixing
     module procedure set_mixing_dp, set_mixing_sp
  end interface set_mixing
  interface get_mixing
     module procedure get_mixing_dp, get_mixing_sp
  end interface get_mixing

  interface set_error
     module procedure set_error_dp, set_error_sp
  end interface set_error
  interface get_error
     module procedure get_error_dp, get_error_sp
  end interface get_error

contains
    
  subroutine set_mixing_dp(mix)
    real(dp),intent(in) :: mix
    mixing = mix
  end subroutine set_mixing_dp
  !
  subroutine set_mixing_sp(mix)
    real(sp),intent(in) :: mix
    mixing = mix
  end subroutine set_mixing_sp

  subroutine get_mixing_dp(mix)
    real(dp),intent(out) :: mix
    mix = mixing
  end subroutine get_mixing_dp
  !
  subroutine get_mixing_sp(mix)
    real(sp),intent(out) :: mix
    mix = mixing
  end subroutine get_mixing_sp

  subroutine set_error_dp(err)
    real(dp),intent(in) :: err
    error = err
  end subroutine set_error_dp
  !
  subroutine set_error_sp(err)
    real(sp),intent(in) :: err
    error = err
  end subroutine set_error_sp

  subroutine get_error_dp(err)
    real(dp),intent(out) :: err
    err = error
  end subroutine get_error_dp
  !
  subroutine get_error_sp(err)
    real(sp),intent(out) :: err
    err = error
  end subroutine get_error_sp

end module SCM_COMMON
