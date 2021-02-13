module mod_precision
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Setup of machine precision
! ------------------------------------------------------------------------------
   use, intrinsic :: iso_fortran_env
   implicit none
   integer, parameter :: sp = REAL32  ! single precision
   integer, parameter :: dp = REAL64  ! double precision
   integer, parameter :: qp = REAL128 ! quadruple precision
   integer, parameter :: wp = REAL64  ! working precision

end module mod_precision