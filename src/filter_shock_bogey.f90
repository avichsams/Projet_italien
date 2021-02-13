subroutine filter_shock_bogey( var )
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Apply adaptative spatial filtering - shock capturing methodology as in
!  Bogey, de Cacqueray, Bailly JCP 2009
!-------------------------------------------------------------------------------
   use mod_filter
   use mod_grid
   use mod_work
   implicit none
   ! ---------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2,nvars), intent(inout) :: var
   ! Local variables
   integer  :: i, j, k, m
   real(wp) :: xp12, xm12, dsnsm1, dsns0, dsnsp1, dsns_mod!, dsens
   real(wp), parameter :: cutoff = 1.e-16_wp
   integer , parameter :: nflt02 = 2
   ! ---------------------------------------------------------------------------

   call update_ghost( div, 0 )

   fltvar = 0.0_wp
   ! x direction
   ! -----------
   do j=1,ny
   do i=0,nx+1
      dsnsm1 = 0.25_wp*(-div(i  ,j) + 2.0_wp*div(i-1,j) - div(i-2,j))
      dsns0  = 0.25_wp*(-div(i+1,j) + 2.0_wp*div(i  ,j) - div(i-1,j))
      dsnsp1 = 0.25_wp*(-div(i+2,j) + 2.0_wp*div(i+1,j) - div(i  ,j))

      dsns_mod = (0.5_wp*((dsns0-dsnsp1)**2 + (dsns0-dsnsm1)**2))/as(i,j)*dx**2 + cutoff

      wrk(i,j) = 0.5_wp*( 1.0_wp-fltshock_amp/dsns_mod + abs(1.0_wp-fltshock_amp/dsns_mod))
   enddo
   enddo

   do m=1,nvars
      do j=1,ny
      do i=ibog1,ibog2
         xp12 = 0.5_wp*(wrk(i,j) + wrk(i+1,j))
         xm12 = 0.5_wp*(wrk(i,j) + wrk(i-1,j))
         fltvar(i,j) = (d12(-1)*var(i-1,j,m) + d12(0)*var(i  ,j,m) +      &
                        d12( 1)*var(i+1,j,m) + d12(2)*var(i+2,j,m))*xp12  &
                      -(d12(-1)*var(i-2,j,m) + d12(0)*var(i-1,j,m) +      &
                        d12( 1)*var(i  ,j,m) + d12(2)*var(i+1,j,m))*xm12
      enddo
      enddo
      var(1:nx,1:ny,m) = var(1:nx,1:ny,m) - fltvar(1:nx,1:ny)
   enddo

   if (ndim.eq.1) return
   ! ---------------------------------------------------------------------------
   ! y direction
   ! -----------
   fltvar = 0.0_wp
   do j=0,ny+1
   do i=1,nx
      dsnsm1 = 0.25_wp*(-div(i,j  ) + 2.0_wp*div(i,j-1) - div(i,j-2))
      dsns0  = 0.25_wp*(-div(i,j+1) + 2.0_wp*div(i,j  ) - div(i,j-1))
      dsnsp1 = 0.25_wp*(-div(i,j+2) + 2.0_wp*div(i,j+1) - div(i,j  ))

      dsns_mod = (0.5_wp*((dsns0-dsnsp1)**2 + (dsns0-dsnsm1)**2))/as(i,j)*dy**2 + cutoff

      wrk(i,j) = 0.5_wp*( 1.0_wp-fltshock_amp/dsns_mod + abs(1.0_wp-fltshock_amp/dsns_mod))
   enddo
   enddo

   do m=1,nvars
      do j=jbog1,jbog2
      do i=1,nx
         xp12 = 0.5_wp*(wrk(i,j) + wrk(i,j+1))
         xm12 = 0.5_wp*(wrk(i,j) + wrk(i,j-1))
         fltvar(i,j) = (d12(-1)*var(i,j-1,m) + d12(0)*var(i,j  ,m) +      &
                          d12( 1)*var(i,j+1,m) + d12(2)*var(i,j+2,m))*xp12  &
                        -(d12(-1)*var(i,j-2,m) + d12(0)*var(i,j-1,m) +      &
                          d12( 1)*var(i,j  ,m) + d12(2)*var(i,j+1,m))*xm12
      enddo
      enddo
      var(1:nx,1:ny,m) = var(1:nx,1:ny,m) - fltvar(1:nx,1:ny)
   enddo

end subroutine filter_shock_bogey