subroutine deriv02_x_exp( var, dvar )
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  compute derivative of array var in the x direction
!-------------------------------------------------------------------------------
   use mod_deriv
   use mod_grid, only: dx
   implicit none
   ! ---------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2), intent(in)  :: var
   real(wp), dimension(  1:nx,  1:ny  ), intent(out) :: dvar
   ! Local variables
   integer  :: i, j
   ! ---------------------------------------------------------------------------
   
   ! X-direction, central differencing
   do j = 1,ny
   do i = ider1_02,ider2_02
      dvar(i,j)  = ( a02(1) * ( var(i+1,j)-var(i-1,j) ) ) / dx
   enddo
   enddo

   ! ---------------------------------------------------------------------------
   if ( is_boundary(1) ) then
      do j = 1,ny
         i = 1
         dvar(i,j)  = ( a02d(1)*var(i  ,j) + a02d(2)*var(i+1,j)  + a02d(3)*var(i+2,j) ) / dx

        

      enddo
   ! ------------------------------------------------------------------------
      do j = 1,ny
         i = nx
         dvar(i,j)  = -( a02d(1)*var(i  ,j) + a02d(2)*var(i-1,j) &
                       + a02d(3)*var(i-2,j)  ) / dx
      enddo
   endif

end subroutine deriv02_x_exp

! ------------------------------------------------------------------------------

subroutine deriv02_y_exp( var, dvar )
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  compute derivative of array var in the y direction
!-------------------------------------------------------------------------------
   use mod_deriv
   use mod_grid, only: dy
   implicit none
   ! ---------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2), intent(in)  :: var
   real(wp), dimension(  1:nx ,  1:ny ), intent(out) :: dvar
   ! Local variables
   integer  :: i, j
   ! ---------------------------------------------------------------------------

   dvar = 0.0_wp
   if (ndim==1) return
   ! Y-direction, central differencing
   do j = jder1_02,jder2_02
   do i = 1,nx
      dvar(i,j)  = ( a02(1) * ( var(i,j+1)-var(i,j-1) )) / dy
   enddo
   enddo

   ! ---------------------------------------------------------------------------
   if ( is_boundary(2) ) then
      do i = 1,nx
         j = 1
         dvar(i,j)  = ( a02d(1)*var(i,j ) + a02d(2)*var(i,j+1) &
                      + a02d(3)*var(i,j+2) ) / dy
      enddo
   ! ------------------------------------------------------------------------
      do i = 1,nx
         j = ny
         dvar(i,j)  = -( a02d(1)*var(i,j) + a02d(2)*var(i,j-1) &
                       + a02d(3)*var(i,j-2) ) / dy
      enddo
   endif

end subroutine deriv02_y_exp