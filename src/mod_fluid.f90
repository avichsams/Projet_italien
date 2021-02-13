module mod_fluid
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Setup and initialization for fluid variables
! ------------------------------------------------------------------------------
   use mod_mode
   implicit none
   ! ---------------------------------------------------------------------------
   real(wp), parameter :: Rair  = 287.15_wp   ! Air gas constant [J/kg K]
   real(wp), parameter :: gamma = 1.4_wp      ! Specific Heat ratio
   real(wp), parameter :: Pr    = 0.7_wp      ! Prandtl number
   real(wp)            :: cv                  ! Isochoric Specific heat
   real(wp)            :: cp                  ! Isobaric Specific heat
   real(wp)            :: mu                  ! Dynamic viscosity
   real(wp)            :: lambda              ! Thermal conductivity
   ! ---------------------------------------------------------------------------

   contains

   subroutine init_fluid
   ! ---------------------------------------------------------------------------
   !  DESCRIPTION
   !  -----------
   !  Derivative variables initialization
   ! ---------------------------------------------------------------------------
      implicit none
      ! Local variables
      ! ------------------------------------------------------------------------

      cv = 1.0_wp/(gamma-1.0_wp)*Rair
      cp =  gamma/(gamma-1.0_wp)*Rair

      if ( is_visc ) then
         mu = 1e-5_wp
         lambda = mu*cp/Pr
      else
         mu = 0.0_wp
         lambda = 0.0_wp
      endif

   end subroutine init_fluid

end module mod_fluid