!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File lsqr_data_m.f90
!
! Defines real(dp) and a few constants for use in other modules.
!
! 24 Oct 2007: Allows floating-point precision dp to be defined
!              in exactly one place (here).  Note that we need
!              use lsqr_data_m at the beginning of modules 
!              AND inside interfaces. Zero and one are not 
!              currently used by LSQR, but this shows how 
!              they should be declared by a user routine that does 
!              need them.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module lsqr_data_m

    implicit none

    intrinsic :: selected_real_kind
    integer, parameter, public :: dp = selected_real_kind(15)
    real(dp), parameter, public :: zero = 0.0_dp, one = 1.0_dp

end module lsqr_data_m
