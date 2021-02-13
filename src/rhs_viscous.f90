subroutine rhs_viscous(phivar, dphidt)
!------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  compute the R.H.S. of the NS equations: viscous fluxes
!-------------------------------------------------------------------------------
   use mod_mode
   use mod_deriv
   use mod_work
   use mod_fluid
   implicit none
   ! ---------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2,4), intent(in)    :: phivar
   real(wp), dimension(  1:nx ,  1:ny ,4), intent(inout) :: dphidt
   ! Local variables
   integer  :: i, j, idir
   real(wp) :: ss11, ss22, ss12
   ! ---------------------------------------------------------------------------

   do j=1,ny
   do i=1,nx
      ss11 =   dvel(i,j,1,1)
      ss22 =   dvel(i,j,2,2)
      ss12 = ( dvel(i,j,1,2) + dvel(i,j,2,1) )

      stress(i,j,1,1) = 2.0_wp*mu*ss11 + (bulk_lad(i,j) - TWO_THIRD*mu)*(ss11+ss22)
      stress(i,j,1,2) = mu*ss12
      stress(i,j,2,1) = mu*ss12
      stress(i,j,2,2) = 2.0_wp*mu*ss22 + (bulk_lad(i,j) - TWO_THIRD*mu)*(ss11+ss22)
   enddo
   enddo

   do j=1,2
   do i=1,2
      call update_ghost( stress(:,:,1,2), 0 )
   enddo
   enddo

   !----------------------------------------------------------------------------
   ! m = 2 rhou1 equation
   !----------------------------------------------------------------------------
   call derivative_x( stress(:,:,1,1), dsxdx, deriv_visc_order )
   call derivative_y( stress(:,:,1,2), dsydy, deriv_visc_order )
   if (is_boundary(1)) call update_neumannbc(dsxdx, dsydy, 1)

   dphidt(:,:,2) = dphidt(:,:,2) + dsxdx + dsydy
   !----------------------------------------------------------------------------
   ! m = 3 rhou2  equation
   !----------------------------------------------------------------------------
   if (ndim==2) then
      call derivative_x( stress(:,:,2,1), dsxdx, deriv_visc_order )
      call derivative_y( stress(:,:,2,2), dsydy, deriv_visc_order )
      if (is_boundary(1)) call update_neumannbc(dsxdx, dsydy, 2)

      dphidt(:,:,3) = dphidt(:,:,3) + dsxdx + dsydy
   endif
   !----------------------------------------------------------------------------
   ! m =4 rhoE equation
   !----------------------------------------------------------------------------
   call update_ghost( tmp, 0 )
   call derivative_x( tmp, qflux(:,:,1), deriv_visc_order )
   call derivative_y( tmp, qflux(:,:,2), deriv_visc_order )

   do idir = 1,2
      do j=1,ny
      do i=1,nx
         qflux(i,j,idir) = -(lambda+lambda_lad(i,j))*qflux(i,j,idir)
      enddo
      enddo
   enddo

   do j=1,ny
   do i=1,nx
      ! x direction
      sx(i,j) = vel(i,j,1)*stress(i,j,1,1) &
              + vel(i,j,2)*stress(i,j,1,2) - qflux(i,j,1)
      ! y direction
      sy(i,j) = vel(i,j,1)*stress(i,j,2,1) &
              + vel(i,j,2)*stress(i,j,2,2) - qflux(i,j,2)
   enddo
   enddo

   call update_ghost( sx, 1 )
   call update_ghost( sy, 2 )

   call derivative_x( sx, dsxdx, deriv_visc_order )
   call derivative_y( sy, dsydy, deriv_visc_order )

   if (is_boundary(1)) call update_neumannbc( dsxdx, dsydy, 0 )

   dphidt(:,:,4) = dphidt(:,:,4) + dsxdx + dsydy

   if (is_boundary(1)) call update_neumannbc_energy(dphidt, dvel, stress)

end subroutine rhs_viscous