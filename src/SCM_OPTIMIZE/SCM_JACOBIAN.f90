module SCM_JACOBIAN
  
  real(8)                         :: df_eps=tiny(1d0)

contains
  
  !----------------------------------------------------------------------
  !           AUXILIARY JACOBIAN/GRADIENT CALCULATIONS
  !
  !          1 x N Jacobian (df_i/dx_j for i=1;j=1,...,N)
  !-----------------------------------------------------------------------
  subroutine fdjac_1n_func(funcv,x,fjac,epsfcn)
    implicit none
    interface 
       function funcv(x)
         implicit none
         real(8),dimension(:) :: x
         real(8)              :: funcv
       end function funcv
    end interface
    integer          ::  n
    real(8)          ::  x(:)
    real(8)          ::  fvec
    real(8)          ::  fjac(size(x))
    real(8),optional ::  epsfcn
    real(8)          ::  eps,eps_
    real(8)          ::  epsmch
    real(8)          ::  h,temp
    real(8)          ::  wa1
    real(8)          ::  wa2
    integer          :: i,j,k
    n=size(x)
    eps_= df_eps; if(present(epsfcn))eps_=epsfcn
    epsmch = epsilon(epsmch)
    eps  = sqrt(max(eps_,epsmch))
    !  Evaluate the function
    fvec = funcv(x)
    do j=1,n
       temp = x(j)
       h    = eps*abs(temp)
       if(h==0.d0) h = eps
       x(j) = temp + h
       wa1  = funcv(x)
       x(j) = temp
       fjac(j) = (wa1-fvec)/h
    enddo
  end subroutine fdjac_1n_func

  subroutine fdjac_1n_sub(funcv,x,fjac,epsfcn)
    implicit none
    interface 
       subroutine funcv(n,x,y)
         implicit none
         integer              :: n
         real(8),dimension(n) :: x
         real(8)              :: y
       end subroutine funcv
    end interface
    integer          ::  n
    real(8)          ::  x(:)
    real(8)          ::  fvec
    real(8)          ::  fjac(size(x))
    real(8),optional ::  epsfcn
    real(8)          ::  eps,eps_
    real(8)          ::  epsmch
    real(8)          ::  h,temp
    real(8)          ::  wa1
    real(8)          ::  wa2
    integer          :: i,j,k
    n=size(x)
    eps_= df_eps; if(present(epsfcn))eps_=epsfcn
    epsmch = epsilon(epsmch)
    eps  = sqrt(max(eps_,epsmch))
    !  Evaluate the function
    call funcv(n,x,fvec)
    !  Computation of dense approximate jacobian.
    do j=1,n
       temp = x(j)
       h    = eps*abs(temp)
       if(h==0.d0) h = eps
       x(j) = temp + h
       call funcv(n,x,wa1)
       x(j) = temp
       fjac(j) = (wa1-fvec)/h
    enddo
    return
  end subroutine fdjac_1n_sub

  function f_jac_1n_func(funcv,n,x,epsfcn) result(df)
    interface
       function funcv(x)
         implicit none
         real(8),dimension(:) :: x
         real(8)              :: funcv
       end function funcv
    end interface
    integer               :: n
    real(8), dimension(n) :: x
    real(8), dimension(n) :: df
    real(8),optional      :: epsfcn
    if(present(epsfcn))then
       call fdjac_1n_func(funcv,x,df,epsfcn=epsfcn)
    else
       call fdjac_1n_func(funcv,x,df)
    end if
  end function f_jac_1n_func

  function f_jac_1n_sub(funcv,n,x,epsfcn) result(df)
    interface
       subroutine funcv(n,x,y)
         implicit none
         integer               :: n
         real(8), dimension(n) :: x
         real(8)               :: y
       end subroutine funcv
    end interface
    integer               :: n
    real(8), dimension(n) :: x
    real(8), dimension(n) :: df
    real(8),optional      :: epsfcn
    if(present(epsfcn))then
       call fdjac_1n_sub(funcv,x,df,epsfcn=epsfcn)
    else
       call fdjac_1n_sub(funcv,x,df)
    end if
  end function f_jac_1n_sub












end module SCM_JACOBIAN
