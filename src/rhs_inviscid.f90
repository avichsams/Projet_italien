subroutine rhs_inviscid(phivar, dphidt)
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  compute the R.H.S. of the NS equations: convective fluxes
!-------------------------------------------------------------------------------
   use mod_mode
   use mod_deriv
   use mod_work
   implicit none
   ! ---------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2,nvars), intent(in)  :: phivar
   real(wp), dimension(  1:nx ,  1:ny ,nvars), intent(out) :: dphidt
   ! Local variables
   integer :: i, j
   ! ---------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! m = 1 density equation
   !----------------------------------------------------------------------------
   call derivative_x( phivar(:,:,2), dsxdx, deriv_conv_order )
   call derivative_y( phivar(:,:,3), dsydy, deriv_conv_order )
   if (is_boundary(1)) call update_neumannbc( dsxdx, dsydy, 0 )

   ! The minus sign is to avoid temporary copy in the previous deriv calls
   dphidt(:,:,1) = -(dsxdx + dsydy)
   !----------------------------------------------------------------------------
   ! m = 2 rhou1 equation
   !----------------------------------------------------------------------------
   do j=ny1,ny2
   do i=nx1,nx2
      sx(i,j) = -(phivar(i,j,2)*vel(i,j,1) + prs(i,j) )
      sy(i,j) = - phivar(i,j,2)*vel(i,j,2)
   enddo
   enddo

   call derivative_x( sx, dsxdx, deriv_conv_order )
   call derivative_y( sy, dsydy, deriv_conv_order )
   if (is_boundary(1)) call update_neumannbc( dsxdx, dsydy, 0 )

   dphidt(:,:,2) = dsxdx + dsydy
   !----------------------------------------------------------------------------
   ! m = 3 rhou2  equation
   !----------------------------------------------------------------------------
   if (ndim==2) then
      do j=ny1,ny2
      do i=nx1,nx2
         sx(i,j) = - phivar(i,j,3)*vel(i,j,1)
         sy(i,j) = -(phivar(i,j,3)*vel(i,j,2) + prs(i,j) )
      enddo
      enddo

      call derivative_x( sx, dsxdx, deriv_conv_order )
      call derivative_y( sy, dsydy, deriv_conv_order )
      if (is_boundary(1)) call update_neumannbc( dsxdx, dsydy, 0 )

      dphidt(:,:,3) = dsxdx + dsydy
   endif
   !----------------------------------------------------------------------------
   ! m = 4 energy equation
   !----------------------------------------------------------------------------
   do j=ny1,ny2
   do i=nx1,nx2
      sx(i,j) = -( phivar(i,j,4) + prs(i,j) )*vel(i,j,1)
      sy(i,j) = -( phivar(i,j,4) + prs(i,j) )*vel(i,j,2)
   enddo
   enddo

   call derivative_x( sx, dsxdx, deriv_conv_order )
   call derivative_y( sy, dsydy, deriv_conv_order )
   if (is_boundary(1)) call update_neumannbc( dsxdx, dsydy, 0 )

   dphidt(:,:,4) = dsxdx + dsydy

end subroutine rhs_inviscid