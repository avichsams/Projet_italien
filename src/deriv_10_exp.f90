subroutine deriv10_x_exp(var, dvar)
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  compute derivative of var in the idir direction
!-------------------------------------------------------------------------------
   use mod_deriv
   use mod_grid, only: dx
   implicit none
   ! ---------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2), intent(in)  :: var
   real(wp), dimension(  1:nx ,  1:ny ), intent(out) :: dvar
   ! Local variables
   integer  :: i, j
   ! ---------------------------------------------------------------------------

   ! X-direction, central differencing
   do j = 1,ny
   do i = ider1_10,ider2_10
      dvar(i,j)  = ( a10(1) * ( var(i+1,j)-var(i-1,j) )   &
                   + a10(2) * ( var(i+2,j)-var(i-2,j) )   &
                   + a10(3) * ( var(i+3,j)-var(i-3,j) )   &
                   + a10(4) * ( var(i+4,j)-var(i-4,j) )   &
                   + a10(5) * ( var(i+5,j)-var(i-5,j) ) ) / dx
   enddo
   enddo

   ! ---------------------------------------------------------------------------
   if ( is_boundary(1) ) then
      do j = 1,ny
         i = 5
         dvar(i,j)  = ( a08(1) * ( var(i+1,j)-var(i-1,j) )   &
                      + a08(2) * ( var(i+2,j)-var(i-2,j) )   &
                      + a08(3) * ( var(i+3,j)-var(i-3,j) )   &
                      + a08(4) * ( var(i+4,j)-var(i-4,j) ) ) / dx

         i = 4
         dvar(i,j)  = ( a06(1) * ( var(i+1,j)-var(i-1,j) )   &
                      + a06(2) * ( var(i+2,j)-var(i-2,j) )   &
                      + a06(3) * ( var(i+3,j)-var(i-3,j) ) ) / dx

         i = 3
         dvar(i,j)  = ( a04(1) * ( var(i+1,j)-var(i-1,j) )   &
                      + a04(2) * ( var(i+2,j)-var(i-2,j) ) ) / dx

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
         i = nx-4
         dvar(i,j)  = ( a08(1) * ( var(i+1,j)-var(i-1,j) )   &
                      + a08(2) * ( var(i+2,j)-var(i-2,j) )   &
                      + a08(3) * ( var(i+3,j)-var(i-3,j) )   &
                      + a08(4) * ( var(i+4,j)-var(i-4,j) ) ) / dx

         i = nx-3
         dvar(i,j)  = ( a06(1) * ( var(i+1,j)-var(i-1,j) )   &
                      + a06(2) * ( var(i+2,j)-var(i-2,j) )   &
                      + a06(3) * ( var(i+3,j)-var(i-3,j) ) ) / dx

         i = nx-2
         dvar(i,j)  = ( a04(1) * ( var(i+1,j)-var(i-1,j) )   &
                      + a04(2) * ( var(i+2,j)-var(i-2,j) ) ) / dx

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

end subroutine deriv10_x_exp

! ------------------------------------------------------------------------------

subroutine deriv10_y_exp( var, dvar )
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  compute derivative of var in the idir direction
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
   do j = jder1_10,jder2_10
   do i = 1,nx
      dvar(i,j)  = ( a10(1) * ( var(i,j+1)-var(i,j-1) )   &
                   + a10(2) * ( var(i,j+2)-var(i,j-2) )   &
                   + a10(3) * ( var(i,j+3)-var(i,j-3) )   &
                   + a10(4) * ( var(i,j+4)-var(i,j-4) )   &
                   + a10(5) * ( var(i,j+5)-var(i,j-5) ) ) / dy
   enddo
   enddo

   ! ---------------------------------------------------------------------------
   if ( is_boundary(2) ) then
      do i = 1,nx
         j = 5
         dvar(i,j) = ( a08(1) * ( var(i,j+1)-var(i,j-1) )   &
                     + a08(2) * ( var(i,j+2)-var(i,j-2) )   &
                     + a08(3) * ( var(i,j+3)-var(i,j-3) )   &
                     + a08(4) * ( var(i,j+4)-var(i,j-4) ) ) / dy

         j = 4
         dvar(i,j)  = ( a06(1) * ( var(i,j+1)-var(i,j-1) )   &
                      + a06(2) * ( var(i,j+2)-var(i,j-2) )   &
                      + a06(3) * ( var(i,j+3)-var(i,j-3) ) ) / dy

         j = 3
         dvar(i,j)  = ( a04(1) * ( var(i,j+1)-var(i,j-1) )   &
                      + a04(2) * ( var(i,j+2)-var(i,j-2) ) ) / dy

         j = 2
         dvar(i,j)  = ( a13d(1)*var(i,j-1) + a13d(2)*var(i,j  ) &
                      + a13d(3)*var(i,j+1) + a13d(4)*var(i,j+2) &
                      + a13d(5)*var(i,j+3) ) / dy

         j = 1
         dvar(i,j)  = ( a04d(1)*var(i,j  ) + a04d(2)*var(i,j+1) &
                      + a04d(3)*var(i,j+2) + a04d(4)*var(i,j+3) &
                      + a04d(5)*var(i,j+4) ) / dy
         ! j = 2
         ! dvar(i,j)  = ( a02(1) * ( var(i,j+1)-var(i,j-1) ) ) / dy
         ! j = 1
         ! dvar(i,j)  = ( var(i,j+1)-var(i,j) ) / dy
      enddo
      ! ------------------------------------------------------------------------
      do i = 1,nx
         j = ny-4
         dvar(i,j)  = ( a08(1) * ( var(i,j+1)-var(i,j-1) )   &
                      + a08(2) * ( var(i,j+2)-var(i,j-2) )   &
                      + a08(3) * ( var(i,j+3)-var(i,j-3) )   &
                      + a08(4) * ( var(i,j+4)-var(i,j-4) ) ) / dy

         j = ny-3
         dvar(i,j)  = ( a06(1) * ( var(i,j+1)-var(i,j-1) )   &
                      + a06(2) * ( var(i,j+2)-var(i,j-2) )   &
                      + a06(3) * ( var(i,j+3)-var(i,j-3) ) ) / dy

         j = ny-2
         dvar(i,j)  = ( a04(1) * ( var(i,j+1)-var(i,j-1) )   &
                      + a04(2) * ( var(i,j+2)-var(i,j-2) ) ) / dy

         j = ny-1
         dvar(i,j)  = -( a13d(1)*var(i,j+1) + a13d(2)*var(i,j  ) &
                       + a13d(3)*var(i,j-1) + a13d(4)*var(i,j-2) &
                       + a13d(5)*var(i,j-3) ) / dy

         j = ny
         dvar(i,j)  = -( a04d(1)*var(i,j  ) + a04d(2)*var(i,j-1) &
                       + a04d(3)*var(i,j-2) + a04d(4)*var(i,j-3) &
                       + a04d(5)*var(i,j-4) ) / dy
         ! j = ny-1
         ! dvar(i,j)  = ( a04(1) * ( var(i,j+1)-var(i,j-1) ) ) / dy
         ! j = ny
         ! dvar(i,j)  = ( var(i,j)-var(i,j-1) ) / dy
      enddo
   endif

end subroutine deriv10_y_exp