module mod_time
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Setup and initialization for time variables
! ------------------------------------------------------------------------------
   use mod_mode
   implicit none
   integer  :: ntime
   integer  :: nsteps
   integer  :: nrk                   ! Type of Runge-Kutta scheme
   integer  :: istage                ! Actual and total number of RK stages
   logical  :: printstep

   real(wp) :: CFL, Fo             ! CFL and Viscous conditions

   real(wp) :: timemax
   real(wp), dimension(:,:), allocatable :: coeff_rk ! RK coefficients

   ! Work arrays
   real(wp), dimension(:,:,:), allocatable :: phitemp, dphidt

   contains

   subroutine init_time
   ! ---------------------------------------------------------------------------
   !  DESCRIPTION
   !  -----------
   !  Time variables initialization
   ! ---------------------------------------------------------------------------
      implicit none
      ! ------------------------------------------------------------------------

      !nrk = 4
      allocate( coeff_rk(nrk,1) )
      coeff_rk = 0.0_wp

      allocate( phitemp(nx1:nx2,ny1:ny2,4) &
              ,  dphidt(  1:nx ,  1:ny ,4) )

      phitemp = 0.0_wp
      dphidt  = 0.0_wp

   end subroutine init_time

end module mod_time