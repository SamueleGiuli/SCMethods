program test_scm_lmixing
  use SCM_COMMON, only: set_mixing, set_error
  use SCM_LMIXING, only: linear_mixing
  use iso_fortran_env, only: dp => real64
  implicit none

  integer  :: i
  real(dp) :: X_new(1), X_old(1), mixing, error
  logical  :: converged
  write(*,*) 'test_scm_lmixing : X = cos(X)'

  X_old=1.0_dp; mixing=0.9_dp; error=1.0d-5
  call linear_mixing(X_old,converged)
  write(*,*) "Using mixing=",mixing
  write(*,*) "Using error=",error
  !
  call set_mixing(mixing)
  call set_error(error)
  !
  do i=1,100
     write(*,*) " "
     X_new = cos(X_old)
     write(*,*) "istep=",i," - X=",X_new
     call linear_mixing(X_new,converged,X_old)
     write(*,*) "diff:",X_new-X_old, " - converged=",converged
     if(converged) exit
  end do
  

end program test_scm_lmixing
