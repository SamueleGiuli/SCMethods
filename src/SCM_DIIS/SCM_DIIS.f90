module SCM_DIIS
  USE SCM_COMMON
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
     module procedure scm_diis_rdp, scm_diis_rsp
  end interface diis

contains

  subroutine scm_diis_rdp(new_xval,new_eval,old_xval,old_eval,print_info_)
    real(dp),intent(inout)                      :: new_xval(:), new_eval(:)
    real(dp),intent(inout),allocatable,optional :: old_xval(:,:),old_eval(:,:)
    logical, optional                           :: print_info_
    logical                                     :: print_info
    integer                                     :: i,j,k, info
    !
    integer,allocatable :: ipiv(:)
    real(dp),allocatable :: wk(:)
    ! size check
    if( (allocated(xvec_rdp).and.(size(xvec_rdp,2)>20)) .or. &
         (allocated(evec_rdp).and.(size(evec_rdp,2)>20)) )then
       if(allocated(xvec_rdp))deallocate(xvec_rdp)
       if(allocated(evec_rdp))deallocate(evec_rdp)
       if(print_info_) write(*,*) "DIIS_size > 20 - discarding all previous elements"
    end if
    

    !
    print_info=.false.; if(present(print_info_)) print_info=print_info_
    if(present(old_xval).and.allocated(xvec_old_rdp)) old_xval = xvec_old_rdp
    if(present(old_eval).and.allocated(evec_old_rdp)) old_eval = evec_old_rdp
    !
    ! if evec_rdp and eve_new_rdp not allocated -> allocate (evec_rdp) and put new_val
    ! else allocate evec_new_rdp and evec_rdp//put new_val -> evec_rdp
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
       do i=1,size(A_rdp,1)
          write(*,*) A_rdp(i,:)
       end do
    end if
    !
    allocate(ipiv(size(A_rdp,2)),wk(1))
    call dsysv('U',size(A_rdp,2),1,A_rdp,size(A_rdp,1),ipiv,B_rdp,size(B_rdp,1),wk,1,info)
    !
    if(print_info)then
       write(*,*) 'C_lambda'
       do i=1,size(B_rdp,1)
          write(*,*) B_rdp(i,1)
       end do
    end if
    !
    if(info/=0) then
       print*,'Error in DIIS_rdp'
       stop
    end if
    new_xval = matmul(xvec_rdp,B_rdp(:size(B_rdp,1)-1,1))
    !
  end subroutine scm_diis_rdp
  !
  subroutine scm_diis_rsp(new_xval,new_eval,old_xval,old_eval,print_info_)
    real(sp),intent(inout)                      :: new_xval(:), new_eval(:)
    real(sp),intent(inout),allocatable,optional :: old_xval(:,:),old_eval(:,:)
    logical,optional                            :: print_info_
    logical                                     :: print_info
    integer                                     :: i,j,k, info
    !
    integer,allocatable :: ipiv(:)
    real(sp),allocatable :: wk(:)
    !
    print_info=.false.; if(present(print_info_)) print_info=print_info_
    if(present(old_xval).and.allocated(xvec_old_rsp)) old_xval = xvec_old_rsp
    if(present(old_eval).and.allocated(evec_old_rsp)) old_eval = evec_old_rsp
    !
    ! if evec_rsp and eve_new_rsp not allocated -> allocate (evec_rsp) and put new_val
    ! else allocate evec_new_rsp and evec_rsp//put new_val -> evec_rsp
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
    A_rsp(size(A_rsp,1),:size(A_rsp,2)-1) = 1.d0
    A_rsp(:size(A_rsp,1)-1,size(A_rsp,2)) = 1.d0
    B_rsp = 0.d0; B_rsp(size(B_rsp,1),1) = 1.d0
    !
    !
    if(print_info)then
       write(*,*) 'B_err'
       do i=1,size(A_rsp,1)
          write(*,*) A_rsp(i,:)
       end do
    end if
    allocate(ipiv(size(A_rsp,2)),wk(1))
    call ssysv('U',size(A_rsp,2),1,A_rsp,size(A_rsp,1),ipiv,B_rsp,size(B_rsp,1),wk,1,info)
    !
    if(print_info)then
       write(*,*) 'C_lambda'
       do i=1,size(B_rsp,1)
          write(*,*) B_rsp(i,1)
       end do
    end if
    !
    if(info/=0) then
       print*,'Error in DIIS_rsp'
       stop
    end if
    new_xval = matmul(xvec_rsp,B_rsp(:size(B_rsp,1)-1,1))
    !
  end subroutine scm_diis_rsp


end module SCM_DIIS
