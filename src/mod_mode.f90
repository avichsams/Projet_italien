module mod_mode
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Setup and initialization for the simulation
! ------------------------------------------------------------------------------
   use mod_precision
   implicit none
   ! ---------------------------------------------------------------------------
   ! natural constants
   real(wp), parameter :: pi    = 3.14159265358979323846_wp ! Pi
   real(wp), parameter :: small = 1.e-10_wp                 ! constant for array init
   real(wp)            :: tiny  = epsilon( 1.0_wp )         ! machine error small
   real(wp)            :: vbig  = huge( 1.0_wp )            ! machine error large

   ! numbers and fractions
   real(wp), parameter :: ZERO = 0.0_wp
   real(wp), parameter :: ONE  = 1.0_wp

   real(wp), parameter :: ONE_THIRD = 1.0_wp/3.0_wp
   real(wp), parameter :: TWO_THIRD = 2.0_wp/3.0_wp
   real(wp), parameter :: root2 = sqrt(2.0_wp), root3 = sqrt(3.0_wp)

   ! flowtype: options for parameter(flowtype=)
   integer, parameter :: SHOCK_TUBE          = 1 &
                       , SHU_OSHER           = 2 &
                       , VORTEX_CONVECTION   = 3 &
                       , SHOCK_VORTEX_INT    = 4 &
                       , TAYLOR_GREEN_VORTEX = 5
   ! ---------------------------------------------------------------------------
   ! face bc types
   logical, dimension(2)   :: is_boundary         ! On/off directional boundary treatment
   integer, parameter      :: PERIODIC   = -1
   integer, parameter      :: OUTFLOW    =  1   ! outflow

   ! ---------------------------------------------------------------------------
   ! filter and derivative type: (filter_type, deriv_conv_type=)
   integer, parameter  :: EXPLICIT_STD  = 1
   integer, parameter  :: EXPLICIT_DRP  = 2
   integer, parameter  :: EXPLICIT_JAM  = 4

   ! time and time step bounds before files system IO breaks
   real(wp), parameter :: max_time = 99999.9990_wp

   ! Parameters of the simulation
   ! ----------------------------
   integer, parameter :: nvars = 4             ! Variables in state vector
   integer            :: ndim
   integer            :: nx, ny                ! local grid size

   integer            :: flowtype              ! Type of flow
   logical            :: is_visc               ! Viscous / inviscid simulation

   logical            :: is_shock              ! Flag to activate shock capturing
   character(len=1 )  :: fltshock_type = 'A'   ! Type of shock capturing

   integer            :: deriv_conv_type       ! Deriv. type for conv. fluxes
   integer            :: deriv_conv_order      ! Deriv. order for conv. fluxes
   integer            :: deriv_visc_type       ! Deriv. type for viscous fluxes
   integer            :: deriv_visc_order      ! Deriv. order for viscous fluxes
   integer            :: filter_order          ! Filter order

   real(wp)           :: dtprint

   integer            :: ngh                   ! Number of ghost cells
   integer            :: nx1, ny1, nz1         ! Lower array bounds
   integer            :: nx2, ny2, nz2         ! Upper array bounds

   logical            :: is_dt_var
   integer            :: printstdout, nprint
   real(wp)           :: time, tstar,  deltat
   real(wp)           :: cputime_ite, cputime_tot
   real(wp)           :: dtstar

   character(len=20)  :: filestamp

   logical            :: is_var_not_updated = .true. ! Check if all arrays are update

   real(wp)           :: Pinf, Tinf, Ainf, Minf, Rhoinf ! Reference p, T, a, M
   real(wp)           :: Uscale, Lscale, Tscale         ! Scale for U, L, T
   real(wp)           :: alpha_dir ! Direction of vortex
   contains

   subroutine init_mode
   ! ---------------------------------------------------------------------------
   !  DESCRIPTION
   !  -----------
   !  Simulation size initialization
   ! ---------------------------------------------------------------------------
      implicit none
      ! ------------------------------------------------------------------------

      ngh = 6

      ndim = 2

      if (ny==1) ndim=1

      nx1 =  1 - ngh
      nx2 = nx + ngh
      ny1 =  1 - ngh
      ny2 = ny + ngh

   end subroutine init_mode

   subroutine calcfilestamp_fromtime(time, filestamp)
   !-------------------------------------------------------------------------------
   !  AUTHOR
   !  ------
   !  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
   !
   !  DESCRIPTION
   !  -----------
   !  Subroutine for calculating filestamp, sequential string.
   !  format 00000.000 ~> 00000_000
   !-------------------------------------------------------------------------------
      use mod_precision
      implicit none
      ! Inputs/Outputs
      real(wp), intent(in)          :: time
      character(len=*), intent(out) :: filestamp
      ! Local variables
      integer :: residual
      ! ---------------------------------------------------------------------------

      residual = int(10000.0_wp*(time-int(time)))
      if (residual < 9998) then
         write( filestamp, '(I4.4,A,I4.4)' ) int(time), '_', residual
      else
         write( filestamp, '(I4.4,A)' ) int(time+1), '_0000'
      endif

   end subroutine calcfilestamp_fromtime

end module mod_mode