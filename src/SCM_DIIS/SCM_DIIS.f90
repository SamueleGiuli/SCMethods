module SCM_DIIS
  !> Module containing the "direct inversion in the iterative subspace" (DIIS) algorithm
  !> The interface require at each call a vector of parameters "new_xval" and a vector of errors "new_eval"
  !> The module will store the last "memory" iterations and will perform a linear combination of the stored vectors
  !> following the DIIS algorithm
  !> The module will return the new_xval vector
  !> If old_xval and/or old_eval are present, the routine will return the old values inside those vectors
  !> The implementation has a maximum number of stored iterations "memory" and a maximum number of iterations
  !> after which it restart the algorithm. There are two types of restart, "best" and "last":
  !> - "best" will restart the algorithm using the best "Nsave" iterations
  !> - "last" will restart the algorithm using the last "Nsave" iterations
  !> default values are memory=10, Nsave=0, restart_type="best"
  !> If print_info is present and true, the routine will print the matrix B_err and C_lambda
  USE SCM_COMMON
  USE STDLIB_IO, only: open
  use iso_fortran_env, only: dp => real64, sp => real32
  implicit none
  
  public :: diis

  private
  
  real(sp),allocatable :: evec_rsp(:,:) , evec_old_rsp(:,:)
  real(sp),allocatable :: xvec_rsp(:,:) , xvec_old_rsp(:,:)
  real(sp),allocatable :: A_rsp(:,:) , B_rsp(:,:)
  !
  real(dp),allocatable :: evec_rdp(:,:) , evec_old_rdp(:,:)
  real(dp),allocatable :: xvec_rdp(:,:) , xvec_old_rdp(:,:)
  real(dp),allocatable :: A_rdp(:,:) , B_rdp(:,:)
  !
  interface diis
     !> Simple DIIS
     module procedure scm_diis_rdp, scm_diis_rsp
  end interface diis

contains

  subroutine scm_diis_rdp(new_xval,new_eval,old_xval,old_eval,print_info_,memory_,Nsave_,restart_type_)
    real(dp),intent(inout)                      :: new_xval(:), new_eval(:)
    real(dp),intent(inout),allocatable,optional :: old_xval(:,:), old_eval(:,:)
    !
    integer, optional                           :: memory_, Nsave_
    character(len=4),optional                   :: restart_type_ !type of restart. default is "best" using the best Nsave elements
    logical, optional                           :: print_info_
    !
    integer                                     :: memory, Nsave
    character(len=4)                            :: restart_type  ! "best" or "last"
    logical                                     :: print_info, desperation
    !
    integer                                     :: i,j,k
    integer                                     :: info, unit
    real(dp)                                    :: err
    !
    integer,allocatable :: ipiv(:)
    real(dp),allocatable :: wk(:),  err_min(:)
    
    print_info=.false.; if(present(print_info_)) print_info=print_info_
    memory=10; if(present(memory_)) memory=memory_
    Nsave=0; if(present(Nsave_)) Nsave=Nsave_
    restart_type="best"; if(present(restart_type_)) restart_type=restart_type_
    !
    if(present(old_xval).and.allocated(xvec_old_rdp)) old_xval = xvec_old_rdp
    if(present(old_eval).and.allocated(evec_old_rdp)) old_eval = evec_old_rdp
    !
    if( (.not.allocated(evec_rdp)) .or. (.not.allocated(xvec_rdp)) ) then
       if(allocated(evec_rdp)) deallocate(evec_rdp)
       if(allocated(xvec_rdp)) deallocate(xvec_rdp)
       allocate(evec_rdp(size(new_eval),1))
       evec_rdp(:,1) = new_eval
       allocate(xvec_rdp(size(new_xval),1))
       xvec_rdp(:,1) = new_xval
    else
       if(allocated(evec_old_rdp)) deallocate(evec_old_rdp)
       if(allocated(xvec_old_rdp)) deallocate(xvec_old_rdp)
       evec_old_rdp = evec_rdp;  xvec_old_rdp = xvec_rdp
       if(allocated(evec_rdp)) deallocate(evec_rdp)
       if(allocated(xvec_rdp)) deallocate(xvec_rdp)
       allocate(evec_rdp(size(new_eval),size(evec_old_rdp,2)+1))
       allocate(xvec_rdp(size(new_xval),size(xvec_old_rdp,2)+1))
       evec_rdp(:,1:size(evec_old_rdp,2)) = evec_old_rdp
       xvec_rdp(:,1:size(xvec_old_rdp,2)) = xvec_old_rdp
       evec_rdp(:,size(evec_rdp,2)) = new_eval
       xvec_rdp(:,size(xvec_rdp,2)) = new_xval
       !
    end if

    ! SIZE check
    if( (allocated(xvec_rdp).and.(size(xvec_rdp,2)>memory)) .or. &
         (allocated(evec_rdp).and.(size(evec_rdp,2)>memory)) )then
       if(Nsave>0)then
          ! Here restart
          if(allocated(xvec_old_rdp)) deallocate(xvec_old_rdp)
          if(allocated(evec_old_rdp)) deallocate(evec_old_rdp)
          if(allocated(err_min)) deallocate(err_min)
          select case(restart_type)
          case("best")
             allocate(xvec_old_rdp(size(xvec_rdp,1),Nsave))
             allocate(evec_old_rdp(size(evec_rdp,1),Nsave))
             allocate(err_min(Nsave))
             err_min(:)=huge(err)
             do i=1,size(evec_rdp,2)
                err=sum(evec_rdp(:,i)**2)
                k=1
                desperation=.true.
                do while(.not.(k>=Nsave) .and. desperation )
                   if( err<err_min(k))then
                      if( k<Nsave )then
                         xvec_old_rdp(:,k+1:) = xvec_old_rdp(:,k:Nsave-1)
                         evec_old_rdp(:,k+1:) = evec_old_rdp(:,k:Nsave-1)
                         err_min(k+1:) = err_min(k:Nsave-1)
                      end if
                      err_min(k)=err
                      xvec_old_rdp(:,k) = xvec_rdp(:,i)
                      evec_old_rdp(:,k) = evec_rdp(:,i)
                      desperation=.false.
                      !
                   end if
                   k=k+1
                end do
             end do
             if(allocated(xvec_rdp))deallocate(xvec_rdp);  if(allocated(evec_rdp))deallocate(evec_rdp)
             xvec_rdp = xvec_old_rdp(:,2:); evec_rdp = evec_old_rdp(:,:)
             if(print_info_) write(*,*) "DIIS_size > ",memory," - discarding all previous elements, restarting from the best ",Nsave
          case("last")
             xvec_old_rdp = xvec_rdp; evec_old_rdp = evec_rdp
             if(allocated(xvec_rdp))deallocate(xvec_rdp); if(allocated(evec_rdp))deallocate(evec_rdp)
             xvec_rdp = xvec_old_rdp; evec_rdp = evec_old_rdp;
             if(allocated(xvec_old_rdp))deallocate(xvec_old_rdp); if(allocated(evec_old_rdp))deallocate(evec_old_rdp)
             if(print_info_) write(*,*) "DIIS_size > ",memory," - discarding all previous elements, restarting from the last ",Nsave
          case default
             print*,"DIIS_size> Error in restart_type"
             stop
          end select
          
       else
          if(allocated(evec_rdp)) deallocate(evec_rdp)
          if(allocated(xvec_rdp)) deallocate(xvec_rdp)
          allocate(evec_rdp(size(new_eval),1))
          evec_rdp(:,1) = new_eval
          allocate(xvec_rdp(size(new_xval),1))
          xvec_rdp(:,1) = new_xval
       end if
    end if
    
    if(allocated(A_rdp)) deallocate(A_rdp)
    if(allocated(B_rdp)) deallocate(B_rdp)
    allocate(A_rdp(size(evec_rdp,2)+1,size(evec_rdp,2)+1))
    allocate(B_rdp(size(evec_rdp,2)+1,1))
    !
    A_rdp=0.d0
    do i=1,size(evec_rdp,2)
       do j=1,size(evec_rdp,2)
          A_rdp(i,j) = sum(evec_rdp(:,i)*evec_rdp(:,j))
       end do
    end do
    A_rdp(size(A_rdp,1),:size(A_rdp,2)-1) = 1.d0
    A_rdp(:size(A_rdp,1)-1,size(A_rdp,2)) = 1.d0
    B_rdp = 0.d0; B_rdp(size(B_rdp,1),1) = 1.d0
    !
    if(print_info)then
       write(*,*) 'B_err'
       unit = open("Berr.mat","w")
       do i=1,size(A_rdp,1)
          write(*,*) A_rdp(i,:)
          write(unit,*) A_rdp(i,:)
       end do
       close(unit)
    end if
    !
    allocate(ipiv(size(A_rdp,2)),wk(1))
    call dsysv('U',size(A_rdp,2),1,A_rdp,size(A_rdp,1),ipiv,B_rdp,size(B_rdp,1),wk,1,info)
    !
    if(print_info)then
       write(*,*) 'C_lambda'
       unit=open("C_lambda.mat","w")
       do i=1,size(B_rdp,1)
          write(*,*) B_rdp(i,1)
          write(unit,*) B_rdp(i,1)
       end do
       close(unit)
    end if
    !
    if(info<0) then
       print*,'Error in DIIS_rdp - info = ',info
       stop
    else if(info>0) then
       print*,'Warning in DIIS_rdp - info = ',info
       print*," discarding previous iterations"
       if(allocated(evec_rdp)) deallocate(evec_rdp)
       if(allocated(xvec_rdp)) deallocate(xvec_rdp)
    else
       new_xval = matmul(xvec_rdp,B_rdp(:size(B_rdp,1)-1,1))
    end if
    !
  end subroutine scm_diis_rdp
  
  subroutine scm_diis_rsp(new_xval,new_eval,old_xval,old_eval,print_info_,memory_,Nsave_,restart_type_)
    real(sp),intent(inout)                      :: new_xval(:), new_eval(:)
    real(sp),intent(inout),allocatable,optional :: old_xval(:,:), old_eval(:,:)
    !
    integer, optional                           :: memory_, Nsave_
    character(len=4),optional                   :: restart_type_ !type of restart. default is "best" using the best Nsave elements
    logical, optional                           :: print_info_
    !
    integer                                     :: memory, Nsave
    character(len=4)                            :: restart_type  ! "best" or "last"
    logical                                     :: print_info, desperation
    !
    integer                                     :: i,j,k
    integer                                     :: info, unit
    real(sp)                                    :: err
    !
    integer,allocatable :: ipiv(:)
    real(sp),allocatable :: wk(:),  err_min(:)
    
    print_info=.false.; if(present(print_info_)) print_info=print_info_
    memory=10; if(present(memory_)) memory=memory_
    Nsave=0; if(present(Nsave_)) Nsave=Nsave_
    restart_type="best"; if(present(restart_type_)) restart_type=restart_type_
    !
    if(present(old_xval).and.allocated(xvec_old_rsp)) old_xval = xvec_old_rsp
    if(present(old_eval).and.allocated(evec_old_rsp)) old_eval = evec_old_rsp
    !
    if( (.not.allocated(evec_rsp)) .or. (.not.allocated(xvec_rsp)) ) then
       if(allocated(evec_rsp)) deallocate(evec_rsp)
       if(allocated(xvec_rsp)) deallocate(xvec_rsp)
       allocate(evec_rsp(size(new_eval),1))
       evec_rsp(:,1) = new_eval
       allocate(xvec_rsp(size(new_xval),1))
       xvec_rsp(:,1) = new_xval
    else
       if(allocated(evec_old_rsp)) deallocate(evec_old_rsp)
       if(allocated(xvec_old_rsp)) deallocate(xvec_old_rsp)
       evec_old_rsp = evec_rsp;  xvec_old_rsp = xvec_rsp
       if(allocated(evec_rsp)) deallocate(evec_rsp)
       if(allocated(xvec_rsp)) deallocate(xvec_rsp)
       allocate(evec_rsp(size(new_eval),size(evec_old_rsp,2)+1))
       allocate(xvec_rsp(size(new_xval),size(xvec_old_rsp,2)+1))
       evec_rsp(:,1:size(evec_old_rsp,2)) = evec_old_rsp
       xvec_rsp(:,1:size(xvec_old_rsp,2)) = xvec_old_rsp
       evec_rsp(:,size(evec_rsp,2)) = new_eval
       xvec_rsp(:,size(xvec_rsp,2)) = new_xval
       !
    end if

    ! SIZE check
    if( (allocated(xvec_rsp).and.(size(xvec_rsp,2)>memory)) .or. &
         (allocated(evec_rsp).and.(size(evec_rsp,2)>memory)) )then
       if(Nsave>0)then
          ! Here restart
          if(allocated(xvec_old_rsp)) deallocate(xvec_old_rsp)
          if(allocated(evec_old_rsp)) deallocate(evec_old_rsp)
          if(allocated(err_min)) deallocate(err_min)
          select case(restart_type)
          case("best")
             allocate(xvec_old_rsp(size(xvec_rsp,1),Nsave))
             allocate(evec_old_rsp(size(evec_rsp,1),Nsave))
             allocate(err_min(Nsave))
             err_min(:)=huge(err)
             do i=1,size(evec_rsp,2)
                err=sum(evec_rsp(:,i)**2)
                k=1
                desperation=.true.
                do while(.not.(k>=Nsave) .and. desperation )
                   if( err<err_min(k))then
                      if( k<Nsave )then
                         xvec_old_rsp(:,k+1:) = xvec_old_rsp(:,k:Nsave-1)
                         evec_old_rsp(:,k+1:) = evec_old_rsp(:,k:Nsave-1)
                         err_min(k+1:) = err_min(k:Nsave-1)
                      end if
                      err_min(k)=err
                      xvec_old_rsp(:,k) = xvec_rsp(:,i)
                      evec_old_rsp(:,k) = evec_rsp(:,i)
                      desperation=.false.
                      !
                   end if
                   k=k+1
                end do
             end do
             if(allocated(xvec_rsp))deallocate(xvec_rsp);  if(allocated(evec_rsp))deallocate(evec_rsp)
             xvec_rsp = xvec_old_rsp(:,2:); evec_rsp = evec_old_rsp(:,:)
             if(print_info_) write(*,*) "DIIS_size > ",memory," - discarding all previous elements, restarting from the best ",Nsave
          case("last")
             xvec_old_rsp = xvec_rsp; evec_old_rsp = evec_rsp
             if(allocated(xvec_rsp))deallocate(xvec_rsp); if(allocated(evec_rsp))deallocate(evec_rsp)
             xvec_rsp = xvec_old_rsp; evec_rsp = evec_old_rsp;
             if(allocated(xvec_old_rsp))deallocate(xvec_old_rsp); if(allocated(evec_old_rsp))deallocate(evec_old_rsp)
             if(print_info_) write(*,*) "DIIS_size > ",memory," - discarding all previous elements, restarting from the last ",Nsave
          case default
             print*,"DIIS_size> Error in restart_type"
             stop
          end select
          
       else
          if(allocated(evec_rsp)) deallocate(evec_rsp)
          if(allocated(xvec_rsp)) deallocate(xvec_rsp)
          allocate(evec_rsp(size(new_eval),1))
          evec_rsp(:,1) = new_eval
          allocate(xvec_rsp(size(new_xval),1))
          xvec_rsp(:,1) = new_xval
       end if
    end if
    
    if(allocated(A_rsp)) deallocate(A_rsp)
    if(allocated(B_rsp)) deallocate(B_rsp)
    allocate(A_rsp(size(evec_rsp,2)+1,size(evec_rsp,2)+1))
    allocate(B_rsp(size(evec_rsp,2)+1,1))
    !
    A_rsp=0.d0
    do i=1,size(evec_rsp,2)
       do j=1,size(evec_rsp,2)
          A_rsp(i,j) = sum(evec_rsp(:,i)*evec_rsp(:,j))
       end do
    end do
    A_rsp(size(A_rsp,1),:size(A_rsp,2)-1) = 1.0
    A_rsp(:size(A_rsp,1)-1,size(A_rsp,2)) = 1.0
    B_rsp = 0.0; B_rsp(size(B_rsp,1),1) = 1.0
    !
    if(print_info)then
       write(*,*) 'B_err'
       unit = open("Berr.mat","w")
       do i=1,size(A_rsp,1)
          write(*,*) A_rsp(i,:)
          write(unit,*) A_rsp(i,:)
       end do
       close(unit)
    end if
    !
    allocate(ipiv(size(A_rsp,2)),wk(1))
    call ssysv('U',size(A_rsp,2),1,A_rsp,size(A_rsp,1),ipiv,B_rsp,size(B_rsp,1),wk,1,info)
    !
    if(print_info)then
       write(*,*) 'C_lambda'
       unit=open("C_lambda.mat","w")
       do i=1,size(B_rsp,1)
          write(*,*) B_rsp(i,1)
          write(unit,*) B_rsp(i,1)
       end do
       close(unit)
    end if
    !
    if(info<0) then
       print*,'Error in DIIS_rsp - info = ',info
       stop
    else if(info>0) then
       print*,'Warning in DIIS_rsp - info = ',info
       print*," discarding previous iterations"
       if(allocated(evec_rsp)) deallocate(evec_rsp)
       if(allocated(xvec_rsp)) deallocate(xvec_rsp)
    else
       new_xval = matmul(xvec_rsp,B_rsp(:size(B_rsp,1)-1,1))
    end if
    !
  end subroutine scm_diis_rsp


end module SCM_DIIS
