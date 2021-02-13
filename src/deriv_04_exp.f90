subroutine deriv04_x_exp( var, dvar )
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
   do i = ider1_04,ider2_04
      dvar(i,j)  = ( a04(1) * ( var(i+1,j)-var(i-1,j) )   &
                   + a04(2) * ( var(i+2,j)-var(i-2,j) ) ) / dx
   enddo
   enddo

   ! ---------------------------------------------------------------------------
   if ( is_boundary(1) ) then
      do j = 1,ny
         i = 2
         dvar(i,j)  = ( a13d(1)*var(i-1,j) + a13d(2)*var(i  ,j) &
                      + a13d(3)*var(i+1,j) + a13d(4)*var(i+2,j) &
                      + a13d(5)*var(i+3,j) ) / dx

         i = 1
         dvar(i,j)  = ( a04d(1)*var(i  ,j) + a04d(2)*var(i+1,j) &
                      + a04d(3)*var(i+2,j) + a04d(4)*var(i+3,j) &
                      + a04d(5)*var(i+4,j) ) / dx
      enddo
   ! ------------------------------------------------------------------------
      do j = 1,ny
         i = nx-1
         dvar(i,j)  = -( a13d(1)*var(i+1,j) + a13d(2)*var(i  ,j) &
                       + a13d(3)*var(i-1,j) + a13d(4)*var(i-2,j) &
                       + a13d(5)*var(i-3,j) ) / dx

         i = nx
         dvar(i,j)  = -( a04d(1)*var(i  ,j) + a04d(2)*var(i-1,j) &
                       + a04d(3)*var(i-2,j) + a04d(4)*var(i-3,j) &
                       + a04d(5)*var(i-4,j) ) / dx
      enddo
   endif

end subroutine deriv04_x_exp

! ------------------------------------------------------------------------------

subroutine deriv04_y_exp( var, dvar )
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
   do j = jder1_04,jder2_04
   do i = 1,nx
      dvar(i,j)  = ( a04(1) * ( var(i,j+1)-var(i,j-1) )   &
                   + a04(2) * ( var(i,j+2)-var(i,j-2) ) ) / dy
   enddo
   enddo

   ! ---------------------------------------------------------------------------
   if ( is_boundary(2) ) then
      do i = 1,nx
         j = 2
         dvar(i,j)  = ( a13d(1)*var(i,j-1) + a13d(2)*var(i,j  ) &
                      + a13d(3)*var(i,j+1) + a13d(4)*var(i,j+2) &
                      + a13d(5)*var(i,j+3) ) / dy

         j = 1
         dvar(i,j)  = ( a04d(1)*var(i,j  ) + a04d(2)*var(i,j+1) &
                      + a04d(3)*var(i,j+2) + a04d(4)*var(i,j+3) &
                      + a04d(5)*var(i,j+4) ) / dy
      enddo
   ! ------------------------------------------------------------------------
      do i = 1,nx
         j = ny-1
         dvar(i,j)  = -( a13d(1)*var(i,j+1) + a13d(2)*var(i,j  ) &
                       + a13d(3)*var(i,j-1) + a13d(4)*var(i,j-2) &
                       + a13d(5)*var(i,j-3) ) / dy

         j = ny
         dvar(i,j)  = -( a04d(1)*var(i,j  ) + a04d(2)*var(i,j-1) &
                       + a04d(3)*var(i,j-2) + a04d(4)*var(i,j-3) &
                       + a04d(5)*var(i,j-4) ) / dy
      enddo
   endif

end subroutine deriv04_y_exp