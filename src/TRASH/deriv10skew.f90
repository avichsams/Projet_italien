subroutine deriv10_skew_x_exp( rho, vel, var, dvar )
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  L. Sciacovelli (luca.sciacovelli@ensam.eu)
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
   real(wp), dimension(nx1:nx2,ny1:ny2), intent(in)  :: rho, vel, var
   real(wp), dimension(  1:nx ,  1:ny ), intent(out) :: dvar
   ! Local variables
   integer  :: i, j
   ! ---------------------------------------------------------------------------

   ! X-direction, central differencing
   do j = 1,ny
   do i = ider1_10,ider2_10
      dvar(i,j) = ( 0.25_wp* ( a10(1) * ( rho(i+1,j)*vel(i+1,j)*var(i+1,j) - rho(i-1,j)*vel(i-1,j)*var(i-1,j) )     &
                               + a10(2) * ( rho(i+2,j)*vel(i+2,j)*var(i+2,j) - rho(i-2,j)*vel(i-2,j)*var(i-2,j) )     &
                               + a10(3) * ( rho(i+3,j)*vel(i+3,j)*var(i+3,j) - rho(i-3,j)*vel(i-3,j)*var(i-3,j) )     &
                               + a10(4) * ( rho(i+4,j)*vel(i+4,j)*var(i+4,j) - rho(i-4,j)*vel(i-4,j)*var(i-4,j) )     &
                               + a10(5) * ( rho(i+5,j)*vel(i+5,j)*var(i+5,j) - rho(i-5,j)*vel(i-5,j)*var(i-5,j) ) )   &
        + 0.25_wp*( vel(i,j)*( a10(1) * ( rho(i+1,j)*             var(i+1,j) - rho(i-1,j)*             var(i-1,j) )     &
                               + a10(2) * ( rho(i+2,j)*             var(i+2,j) - rho(i-2,j)*             var(i-2,j) )     &
                               + a10(3) * ( rho(i+3,j)*             var(i+3,j) - rho(i-3,j)*             var(i-3,j) )     &
                               + a10(4) * ( rho(i+4,j)*             var(i+4,j) - rho(i-4,j)*             var(i-4,j) )     &
                               + a10(5) * ( rho(i+5,j)*             var(i+5,j) - rho(i-5,j)*             var(i-5,j) ) )   &
                  + rho(i,j)*( a10(1) * (              vel(i+1,j)*var(i+1,j) -              vel(i-1,j)*var(i-1,j) )     &
                               + a10(2) * (              vel(i+2,j)*var(i+2,j) -              vel(i-2,j)*var(i-2,j) )     &
                               + a10(3) * (              vel(i+3,j)*var(i+3,j) -              vel(i-3,j)*var(i-3,j) )     &
                               + a10(4) * (              vel(i+4,j)*var(i+4,j) -              vel(i-4,j)*var(i-4,j) )     &
                               + a10(5) * (              vel(i+5,j)*var(i+5,j) -              vel(i-5,j)*var(i-5,j) ) )   &
                  + var(i,j)*( a10(1) * ( rho(i+1,j)*vel(i+1,j)              - rho(i-1,j)*vel(i-1,j)              )     &
                               + a10(2) * ( rho(i+2,j)*vel(i+2,j)              - rho(i-2,j)*vel(i-2,j)              )     &
                               + a10(3) * ( rho(i+3,j)*vel(i+3,j)              - rho(i-3,j)*vel(i-3,j)              )     &
                               + a10(4) * ( rho(i+4,j)*vel(i+4,j)              - rho(i-4,j)*vel(i-4,j)              )     &
                               + a10(5) * ( rho(i+5,j)*vel(i+5,j)              - rho(i-5,j)*vel(i-5,j)              ) ) ) &
        + 0.25_wp*( rho(i,j)*vel(i,j)*( a10(1)*(var(i+1,j) - var(i-1,j) )   &
                                          + a10(2)*(var(i+2,j) - var(i-2,j) )   &
                                          + a10(3)*(var(i+3,j) - var(i-3,j) )   &
                                          + a10(4)*(var(i+4,j) - var(i-4,j) )   &
                                          + a10(5)*(var(i+5,j) - var(i-5,j) ) ) &
                  + rho(i,j)*var(i,j)*( a10(1)*(vel(i+1,j) - vel(i-1,j) )   &
                                          + a10(2)*(vel(i+2,j) - vel(i-2,j) )   &
                                          + a10(3)*(vel(i+3,j) - vel(i-3,j) )   &
                                          + a10(4)*(vel(i+4,j) - vel(i-4,j) )   &
                                          + a10(5)*(vel(i+5,j) - vel(i-5,j) ) ) &
                  + vel(i,j)*var(i,j)*( a10(1)*(rho(i+1,j) - rho(i-1,j) )   &
                                          + a10(2)*(rho(i+2,j) - rho(i-2,j) )   &
                                          + a10(3)*(rho(i+3,j) - rho(i-3,j) )   &
                                          + a10(4)*(rho(i+4,j) - rho(i-4,j) )   &
                                          + a10(5)*(rho(i+5,j) - rho(i-5,j) ) ) ) )/dx
   enddo
   enddo

end subroutine deriv10_skew_x_exp

! ------------------------------------------------------------------------------

subroutine deriv10_skew_y_exp( rho, vel, var, dvar )
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  L. Sciacovelli (luca.sciacovelli@ensam.eu)
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
   real(wp), dimension(nx1:nx2,ny1:ny2), intent(in)  :: rho, vel, var
   real(wp), dimension(  1:nx ,  1:ny ), intent(out) :: dvar
   ! Local variables
   integer  :: i, j, k
   ! ---------------------------------------------------------------------------

   ! Y-direction, central differencing
   do j = jder1_10,jder2_10
   do i = 1,nx
      dvar(i,j) = ( 0.25_wp* ( a10(1) * ( rho(i,j+1)*vel(i,j+1)*var(i,j+1) - rho(i,j-1)*vel(i,j-1)*var(i,j-1) )     &
                               + a10(2) * ( rho(i,j+2)*vel(i,j+2)*var(i,j+2) - rho(i,j-2)*vel(i,j-2)*var(i,j-2) )     &
                               + a10(3) * ( rho(i,j+3)*vel(i,j+3)*var(i,j+3) - rho(i,j-3)*vel(i,j-3)*var(i,j-3) )     &
                               + a10(4) * ( rho(i,j+4)*vel(i,j+4)*var(i,j+4) - rho(i,j-4)*vel(i,j-4)*var(i,j-4) )     &
                               + a10(5) * ( rho(i,j+5)*vel(i,j+5)*var(i,j+5) - rho(i,j-5)*vel(i,j-5)*var(i,j-5) ) )   &
        + 0.25_wp*( vel(i,j)*( a10(1) * ( rho(i,j+1)*             var(i,j+1) - rho(i,j-1)*             var(i,j-1) )     &
                               + a10(2) * ( rho(i,j+2)*             var(i,j+2) - rho(i,j-2)*             var(i,j-2) )     &
                               + a10(3) * ( rho(i,j+3)*             var(i,j+3) - rho(i,j-3)*             var(i,j-3) )     &
                               + a10(4) * ( rho(i,j+4)*             var(i,j+4) - rho(i,j-4)*             var(i,j-4) )     &
                               + a10(5) * ( rho(i,j+5)*             var(i,j+5) - rho(i,j-5)*             var(i,j-5) ) )   &
                  + rho(i,j)*( a10(1) * (              vel(i,j+1)*var(i,j+1) -              vel(i,j-1)*var(i,j-1) )     &
                               + a10(2) * (              vel(i,j+2)*var(i,j+2) -              vel(i,j-2)*var(i,j-2) )     &
                               + a10(3) * (              vel(i,j+3)*var(i,j+3) -              vel(i,j-3)*var(i,j-3) )     &
                               + a10(4) * (              vel(i,j+4)*var(i,j+4) -              vel(i,j-4)*var(i,j-4) )     &
                               + a10(5) * (              vel(i,j+5)*var(i,j+5) -              vel(i,j-5)*var(i,j-5) ) )   &
                  + var(i,j)*( a10(1) * ( rho(i,j+1)*vel(i,j+1)              - rho(i,j-1)*vel(i,j-1)              )     &
                               + a10(2) * ( rho(i,j+2)*vel(i,j+2)              - rho(i,j-2)*vel(i,j-2)              )     &
                               + a10(3) * ( rho(i,j+3)*vel(i,j+3)              - rho(i,j-3)*vel(i,j-3)              )     &
                               + a10(4) * ( rho(i,j+4)*vel(i,j+4)              - rho(i,j-4)*vel(i,j-4)              )     &
                               + a10(5) * ( rho(i,j+5)*vel(i,j+5)              - rho(i,j-5)*vel(i,j-5)              ) ) ) &
        + 0.25_wp*( rho(i,j)*vel(i,j)*( a10(1)*(var(i,j+1) - var(i,j-1) )   &
                                          + a10(2)*(var(i,j+2) - var(i,j-2) )   &
                                          + a10(3)*(var(i,j+3) - var(i,j-3) )   &
                                          + a10(4)*(var(i,j+4) - var(i,j-4) )   &
                                          + a10(5)*(var(i,j+5) - var(i,j-5) ) ) &
                  + rho(i,j)*var(i,j)*( a10(1)*(vel(i,j+1) - vel(i,j-1) )   &
                                          + a10(2)*(vel(i,j+2) - vel(i,j-2) )   &
                                          + a10(3)*(vel(i,j+3) - vel(i,j-3) )   &
                                          + a10(4)*(vel(i,j+4) - vel(i,j-4) )   &
                                          + a10(5)*(vel(i,j+5) - vel(i,j-5) ) ) &
                  + vel(i,j)*var(i,j)*( a10(1)*(rho(i,j+1) - rho(i,j-1) )   &
                                          + a10(2)*(rho(i,j+2) - rho(i,j-2) )   &
                                          + a10(3)*(rho(i,j+3) - rho(i,j-3) )   &
                                          + a10(4)*(rho(i,j+4) - rho(i,j-4) )   &
                                          + a10(5)*(rho(i,j+5) - rho(i,j-5) ) ) ) )/ dy
   enddo
   enddo

end subroutine deriv10_skew_y_exp