subroutine update_dirichbc( phivar )
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Subroutine to update Dirichlet bcs on conservative variables.
!-------------------------------------------------------------------------------
   use mod_mode
   implicit none
   ! ---------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2,nvars), intent(inout) :: phivar
   ! ---------------------------------------------------------------------------

   if ( is_boundary(1) ) then
      phivar( 1,1:ny,:) = (2.0_wp* phivar(   2,1:ny,:) + phivar(   3,1:ny,:))/3.0_wp
      phivar(nx,1:ny,:) = (2.0_wp* phivar(nx-1,1:ny,:) + phivar(nx-2,1:ny,:))/3.0_wp
   endif

   if ( is_boundary(2) ) then
      phivar(1:nx, 1,:) = (2.0_wp* phivar(1:nx,   2,:) + phivar(1:nx,   3,:))/3.0_wp
      phivar(1:nx,ny,:) = (2.0_wp* phivar(1:nx,ny-1,:) + phivar(1:nx,ny-2,:))/3.0_wp
   endif

end subroutine update_dirichbc