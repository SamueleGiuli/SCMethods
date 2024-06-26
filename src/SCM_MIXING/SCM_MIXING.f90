module SCM_LMIXING
  use SCM_COMMON
  use iso_fortran_env, only: dp => real64, sp => real32
  implicit none
  
  public :: linear_mixing

  private
  
  real(sp),allocatable :: prev_rsp_1(:), prev_rsp_2(:,:)
  real(dp),allocatable :: prev_rdp_1(:), prev_rdp_2(:,:)

  interface linear_mixing
     module procedure scm_lmixing_rsp_1, scm_lmixing_rdp_1
     module procedure scm_lmixing_rsp_2, scm_lmixing_rdp_2
  end interface linear_mixing

contains
  
  subroutine scm_lmixing_rdp_1(new_val,convg,old_val)
    real(dp),intent(inout)          :: new_val(:)
    logical,intent(out)             :: convg
    real(dp),intent(inout),optional :: old_val(:)
    !
    convg=.false.
    if(present(old_val)) old_val = new_val
    !
    if (.not.allocated(prev_rdp_1)) then
       allocate(prev_rdp_1(size(new_val)))
       prev_rdp_1 = new_val
    else
       new_val = mixing*new_val + (1.d0-mixing)*prev_rdp_1
       if(all(abs(new_val-prev_rdp_1)<error)) convg=.true.
       prev_rdp_1 = new_val
    end if
    !
  end subroutine scm_lmixing_rdp_1
  !
  subroutine scm_lmixing_rsp_1(new_val,convg,old_val)
    real(sp),intent(inout)          :: new_val(:)
    logical,intent(out)             :: convg
    real(sp),intent(inout),optional :: old_val(:)
    !
    convg=.false.
    if(present(old_val)) old_val = new_val
    !
    if (.not.allocated(prev_rsp_1)) then
       allocate(prev_rsp_1(size(new_val)))
       prev_rsp_1 = new_val
    else
       new_val = mixing*new_val + (1.d0-mixing)*prev_rsp_1
       if(all(abs(new_val-prev_rsp_1)<error)) convg=.true.
       prev_rsp_1 = new_val
    end if
    !
  end subroutine scm_lmixing_rsp_1

    subroutine scm_lmixing_rdp_2(new_val,convg,old_val)
    real(dp),intent(inout)          :: new_val(:,:)
    logical,intent(out)             :: convg
    real(dp),intent(inout),optional :: old_val(:,:)
    !
    convg=.false.
    if(present(old_val)) old_val = new_val
    !
    if (.not.allocated(prev_rdp_2)) then
       allocate(prev_rdp_2(size(new_val,1),size(new_val,2)))
       prev_rdp_2 = new_val
    else
       new_val = mixing*new_val + (1.d0-mixing)*prev_rdp_2
       if(all(abs(new_val-prev_rdp_2)<error)) convg=.true.
       prev_rdp_2 = new_val
    end if
    !
  end subroutine scm_lmixing_rdp_2
  !
  subroutine scm_lmixing_rsp_2(new_val,convg,old_val)
    real(sp),intent(inout)          :: new_val(:,:)
    logical,intent(out)             :: convg
    real(sp),intent(inout),optional :: old_val(:,:)
    !
    convg=.false.
    if(present(old_val)) old_val = new_val
    !
    if (.not.allocated(prev_rsp_2)) then
       allocate(prev_rsp_2(size(new_val,1),size(new_val,2)))
       prev_rsp_2 = new_val
    else
       new_val = mixing*new_val + (1.d0-mixing)*prev_rsp_2
       if(all(abs(new_val-prev_rsp_2)<error)) convg=.true.
       prev_rsp_2 = new_val
    end if
    !
  end subroutine scm_lmixing_rsp_2

end module SCM_LMIXING
