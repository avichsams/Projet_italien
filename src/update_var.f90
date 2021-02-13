subroutine update_var( phivar )
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Routine to update properties in RK.
!-------------------------------------------------------------------------------
   use mod_work
   implicit none
   ! ---------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2,4), intent(inout) :: phivar
   ! Local variables
   integer :: i, j, n
   ! ---------------------------------------------------------------------------

   ! First, apply boundary conditions on conservative variables
   call update_dirichbc( phivar )

   ! Compute thermodynamic properties over the whole domain
   call compute_thermo( phivar )

   call update_ghost( prs, 0 )
   do n=1,nvars
      call update_ghost( phivar(:,:,n), 0 )
   enddo

   ! Update gradients
   ! Update vel and yspv in ghost cells too
   ! --------------------------------------
   do j=ny1,ny2
   do i=nx1,nx2
      vel(i,j,1)  = phivar(i,j,2)/phivar(i,j,1)
      vel(i,j,2)  = phivar(i,j,3)/phivar(i,j,1)
   enddo
   enddo

   do n = 1,2
      call derivative_x( vel(:,:,n), dvel(:,:,n,1), deriv_visc_order )
      call derivative_y( vel(:,:,n), dvel(:,:,n,2), deriv_visc_order )
   enddo

   do j=1,ny
   do i=1,nx
      div(i,j) = dvel(i,j,1,1) + dvel(i,j,2,2)
   enddo
   enddo

   is_var_not_updated = .false.

end subroutine update_var