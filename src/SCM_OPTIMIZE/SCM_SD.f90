module SCM_SD

  USE SCM_JACOBIAN
  
  abstract interface
     function sdfit_func(a)
       real(8),dimension(:)  ::  a
       real(8)               ::  sdfit_func
     end function sdfit_func
     function sdfit_fjac(a)
       real(8),dimension(:)       :: a
       real(8),dimension(size(a)) :: sdfit_fjac
     end function sdfit_fjac
  end interface

  integer                        :: ncom
  real(8), dimension(:), pointer :: pcom,xicom
  procedure(sdfit_func),pointer  :: func
  procedure(sdfit_fjac),pointer  :: fjac
  real(8), dimension(:),pointer   :: fmin_fvecp

contains

  ! Steepest Descent
  subroutine fmin_sd_f(p,f,iter,fret,ftol,itmax,istop,deps,iverbose)
    procedure(sdfit_func)                :: f
    real(8), dimension(:), intent(inout) :: p
    integer, intent(out)                 :: iter
    real(8), intent(out)                 :: fret
    real(8),optional                     :: ftol,deps
    integer, optional                    :: itmax,istop
    real(8)                              :: ftol_,deps_,a,b
    integer                              :: itmax_,istop_
    integer                              :: i,its
    real(8)                              :: dgg,fp,gam,gg,err_
    real(8), dimension(size(p))          :: g,h,xi,p_prev
    logical,optional                     :: iverbose
    logical                              :: iverbose_,converged
    !
    !this is just to ensure that routine needing dfunc allocated
    !and properly definted will continue to work.
    if(associated(func))nullify(func) ; func=>f
    if(associated(fjac))nullify(fjac) ; fjac=>df
    !
    iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
    ftol_=1.d-5
    if(present(ftol))then
       ftol_=ftol
       if(iverbose_)write(*,"(A,ES9.2)")"SD: ftol updated to:",ftol
    endif
    itmax_=500
    if(present(itmax))then
       itmax_=itmax
       if(iverbose_)write(*,"(A,I5)")"SD: itmax updated to:",itmax
    endif
    istop_=0
    if(present(istop))then
       istop_=istop
       if(iverbose_)write(*,"(A,I3)")"SD: istop update to:",istop
    endif
    deps_=1.d-8
    if(present(deps))then
       deps_=deps
       if(iverbose_)write(*,"(A,ES9.2)")"SD: deps update to:",deps
    endif
    df_eps = deps_
    !
    fp=func(p)
    xi=fjac(p)
    !
    do its=1,itmax_
       print*,"its",its," - xi:",xi
       fret = fp
       p_prev = p
       g    =-xi
       p    = p + g
       fp   = func(p)
       !
       a = abs(fret-fp)/(1d0+abs(fp))
       b = dot_product(p-p_prev,p-p_prev)/(1d0+dot_product(p,p))
       !
       if(iverbose_)then
          write(*,*)"Iter,F_n =",iter,fp
          write(*,"(A10)")"    gradF:"
          do i=1,size(xi)
             write(*,*)xi(i)
          enddo
          write(*,*)"A,B,ftol =",a,b,ftol_
       endif
       select case(istop_)
       case default
          converged = (a<ftol_) .AND. (b<ftol_)
          if( converged )print*,"Converged with (|F_n-F_n-1|/1+|F_n|), ||a_n - a_n-1||^2/1+||a_n||^2 < ftol_:",a,b
       case(1)
          converged = a<ftol_
          if( converged )print*,"Converged with (|F_n-F_n-1|/1+|F_n|)< ftol_:",a
       case(2)
          converged = b<ftol_
          if( converged )print*,"Converged with ||a_n - a_n-1||^2/1+||a_n||^2 < ftol_:",b
       end select
       if( converged )return
       if( iverbose_)write(*,*)""
       xi   = fjac(p)
    end do

    if(iverbose_)write(*,*)"SD: MaxIter",itmax_," exceeded."
    nullify(func)
    nullify(fjac)
    return
  end subroutine fmin_sd_f
  !
  function df(p)
    real(8),dimension(:)       :: p
    real(8),dimension(size(p)) :: df
    df=f_jac_1n_func(func,size(p),p)
  end function df

end module SCM_SD
