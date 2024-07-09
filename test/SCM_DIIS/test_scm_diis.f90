program test_scm_diis
  use SCM_DIIS, only: diis
  use iso_fortran_env, only: dp => real64
  implicit none

  integer  :: i
  real(dp) :: X_new(1), X_old(1), E_new(1), mixing, error(1)
  logical  :: converged
  write(*,*) 'test_scm_diis : X = cos(X)'

  X_new=1.0_dp; mixing=0.9_dp; error=1.0d-6
  write(*,*) "Using mixing=",mixing
  write(*,*) "Using error=",error
  X_old = X_new
  X_new = cos(X_old)
  write(*,*) "istep=0 - X=",X_new
  X_old = X_new
  X_new = cos(X_old)
  !
  !
  do i=1,100
     write(*,*) " "
     X_old = X_new
     X_new = cos(X_old)
     E_new=X_new-X_old
     write(*,*) "istep=",i," - X=",X_new
     call diis(X_new,E_new,print_info_=.true.)
     if( all(abs(X_new-X_old) < error) ) converged=.true.
     write(*,*) "diff:",X_new-X_old, " - converged=",converged
     if(all(abs(E_new)<error)) exit
  end do
  

end program test_scm_diis
