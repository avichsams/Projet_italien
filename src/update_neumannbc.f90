subroutine update_neumannbc(dsxdx, dsydy, itag)
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Subroutine to update Neumann bc variables
!-------------------------------------------------------------------------------
   use mod_mode
   implicit none
   ! ---------------------------------------------------------------------------
   ! Input/Output arguments
   integer, intent(in) :: itag
   real(wp), dimension(1:nx,1:ny), intent(inout) :: dsxdx, dsydy
   ! ---------------------------------------------------------------------------

   select case ( itag )
      case ( 0 )
         ! Inviscid / viscous rhoE condition
         if ( is_boundary(1) ) then
            dsxdx( 1,1:ny) = 0.0_wp
            dsxdx(nx,1:ny) = 0.0_wp
         endif

         if ( is_boundary(2) ) then
            dsydy(1:nx, 1) = 0.0_wp
            dsydy(1:nx,ny) = 0.0_wp
         endif

      case ( 1 )
         ! viscous rho u1 condition
         if ( is_boundary(2) ) then
            dsydy(1:nx, 1) = 0.0_wp
            dsydy(1:nx,ny) = 0.0_wp
         endif

      case ( 2 )
         ! viscous rho u2 condition
         if ( is_boundary(1) ) then
            dsxdx( 1,1:ny) = 0.0_wp
            dsxdx(nx,1:ny) = 0.0_wp
         endif

      case default
         write(*,*) 'itag incorrect in update_neumannbc.f90'
         stop
   end select

end subroutine update_neumannbc

!-------------------------------------------------------------------------------

subroutine update_neumannbc_energy(dphidt, dvel, stress)
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Subroutine to update Neumann bc variables in energy equation
!-------------------------------------------------------------------------------
   use mod_mode
   implicit none
   ! ---------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2,2,2), intent(in) :: stress, dvel
   real(wp), dimension(  1:nx ,  1:ny ,4)  , intent(inout) :: dphidt
   ! ---------------------------------------------------------------------------

   ! Handling outflow
   ! ----------------
   if ( is_boundary(1) ) then
      dphidt(   1,1:ny,4) = dphidt( 1,1:ny,4) + dvel( 1,1:ny,2,1)*stress( 1,1:ny,1,2)
      dphidt(  nx,1:ny,4) = dphidt(nx,1:ny,4) + dvel(nx,1:ny,2,1)*stress(nx,1:ny,1,2)
   endif
   if ( is_boundary(2) ) then
      dphidt(1:nx,   1,4) = dphidt(1:nx, 1,4) + dvel(1:nx, 1,1,2)*stress(1:nx, 1,2,1)
      dphidt(1:nx,  ny,4) = dphidt(1:nx,ny,4) + dvel(1:nx,ny,1,2)*stress(1:nx,ny,2,1)
   endif

end subroutine update_neumannbc_energy