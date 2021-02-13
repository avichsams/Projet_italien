module mod_grid
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Setup and initialization for grid variables
! ------------------------------------------------------------------------------
   use mod_mode
   implicit none
   ! ---------------------------------------------------------------------------
   real(wp)              :: Xlength, xmin, xmax, dx
   real(wp)              :: Ylength, ymin, ymax, dy
   real(wp), allocatable :: gridx(:), gridy(:)
   ! ---------------------------------------------------------------------------

   contains

   subroutine init_grid
   ! ---------------------------------------------------------------------------
   !  DESCRIPTION
   !  -----------
   !  Grid variables initialization
   ! ---------------------------------------------------------------------------
      implicit none
      ! Local variables
      integer :: i, j
      ! ------------------------------------------------------------------------
      ! Init grid arrays
      allocate( gridx(1:nx),  gridy(1:ny) )

      gridx(:) = 0.0_wp
      gridy(:) = 0.0_wp

      ! Minimum/maximum extrema in the three directions, General case
      xmin = 0.0_wp
      xmax = xmin + Xlength
      ymin = 0.0_wp
      ymax = ymin + Ylength

      select case ( flowtype )
         case ( SHOCK_VORTEX_INT )
            xmin = -30.0_wp
            xmax = xmin + Xlength
            ymin = -45.0_wp
            ymax = ymin + Ylength

         case ( SHOCK_TUBE, SHU_OSHER )
            xmin = -Xlength/2.0_wp
            xmax = xmin + Xlength

         case ( TAYLOR_GREEN_VORTEX )
            Xlength = Xlength*pi
            xmax = xmin + Xlength
            Ylength = Ylength*pi
            ymax = ymin + Ylength
         case default
            ! Nothing to do
      end select

      ! Set correct boundary (periodic or not) according to test case
      select case ( flowtype )
         case ( SHOCK_TUBE, SHOCK_VORTEX_INT, SHU_OSHER )
            is_boundary(1) = .true.
            is_boundary(2) = .false.
         case default
         is_boundary(:) = .false.
      end select

      ! ------------------------------------------------------------------------
      ! generate uniform physical grids
      if (is_boundary(1)) then
         dx = Xlength/dble(nx-1)
      else
         dx = Xlength/dble(nx)
      endif
      do i=1, nx
         gridx(i) = xmin + dx*(i-1)
      enddo

      if ( is_boundary(2) ) then
         dy = Ylength/dble(ny-1)
      else
         dy = Ylength/dble(ny)
      endif
      do j = 1,ny
         gridy(j) = ymin + dy*(j-1)
      enddo

   end subroutine init_grid

end module mod_grid