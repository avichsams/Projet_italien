subroutine filter_shock_jameson( var )
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Apply artificial viscosity (DNC-Jameson)
!-------------------------------------------------------------------------------
   use mod_filter
   use mod_grid, only: dx, dy
   use mod_work, only: vel, dvel, as, prs, div
   implicit none
   ! ---------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2,nvars), intent(inout) :: var
   ! Local variables
   real(wp), parameter :: cutoff = 1.e-15_wp
   real(wp), parameter :: coeff_corr10 = 12.0_wp/1260.0_wp &
                        , coeff_corr8  = 12.0_wp/ 280.0_wp &
                        , coeff_corr6  = 12.0_wp/  60.0_wp &
                        , coeff_corr4  = 12.0_wp/  12.0_wp
   integer  :: i, j, m
   real(wp) :: b3, xnu2, xnu4, xnu6, xnu8, xnu10
   real(wp) :: dh_0, x10p12, x10m12, lp12, lm12, xp12, xm12
   real(wp) :: dh_m5, dh_m4, dh_m3, dh_m2, dh_m1
   real(wp) :: dh_p5, dh_p4, dh_p3, dh_p2, dh_p1
   ! ---------------------------------------------------------------------------

   ! Compute Ducros sensor
   ! ---------------------
   do j=1,ny
   do i=1,nx
      if (div(i,j) < 0.0_wp) then
         b3 = (dvel(i,j,2,1) - dvel(i,j,1,2))**2
         ducros(i,j) = div(i,j)**2 / (div(i,j)**2 + b3 + cutoff)
      else
         ducros(i,j) = 0.0_wp
      endif
      rspec(i,j,1) = abs(vel(i,j,1)) + as(i,j)
      rspec(i,j,2) = abs(vel(i,j,2)) + as(i,j)
   enddo
   enddo

   ! Pressure is already known in ghostcells, no need to communicate
   do j=1,ny
   do i=1,nx
      psens(i,j,1) = abs( (prs(i-1,j)-2.0_wp*prs(i,j)+prs(i+1,j)) / &
                          (prs(i-1,j)+2.0_wp*prs(i,j)+prs(i+1,j)) ) * ducros(i,j)
   enddo
   enddo

   do j=1,ny
   do i=1,nx
      psens(i,j,2) = abs( (prs(i,j-1)-2.0_wp*prs(i,j)+prs(i,j+1)) / &
                          (prs(i,j-1)+2.0_wp*prs(i,j)+prs(i,j+1)) ) * ducros(i,j)
   enddo
   enddo

   do m=1,2
      call update_ghost(psens(:,:,m), m)
      call update_ghost(rspec(:,:,m), m)
   enddo

   xnu10 = fltamp/1260.0_wp
   xnu8  = fltamp/ 280.0_wp
   xnu6  = fltamp/  60.0_wp
   xnu4  = fltamp/  12.0_wp
   xnu2  = fltamp/   2.0_wp

   Varsloop: do m=1,nvars
      ! ------------------------------------------------------------------------
      ! x direction
      ! -----------
      ! no boundaries use ghost cell
      do j=1,ny
      do i=iflt1,iflt2
         lp12 = 0.5_wp*( rspec(i,j,1) + rspec(i+1,j,1) )/dx
         lm12 = 0.5_wp*( rspec(i,j,1) + rspec(i-1,j,1) )/dx

         xp12 = fltshock_amp*max( psens(i,j,1), psens(i+1,j,1) )
         xm12 = fltshock_amp*max( psens(i,j,1), psens(i-1,j,1) )
         x10p12 = lp12*max(0.0_wp,xnu10-coeff_corr10*xp12)
         x10m12 = lm12*max(0.0_wp,xnu10-coeff_corr10*xm12)
         xp12 = xp12*lp12
         xm12 = xm12*lm12

         dh_p5=          x10p12
         dh_p4=  -9.0_wp*x10p12         -x10m12
         dh_p3=  36.0_wp*x10p12  +9.0_wp*x10m12
         dh_p2= -84.0_wp*x10p12 -36.0_wp*x10m12
         dh_p1= 126.0_wp*x10p12 +84.0_wp*x10m12
         dh_0 =-126.0_wp*x10p12-126.0_wp*x10m12
         dh_m1=  84.0_wp*x10p12+126.0_wp*x10m12
         dh_m2= -36.0_wp*x10p12 -84.0_wp*x10m12
         dh_m3=   9.0_wp*x10p12 +36.0_wp*x10m12
         dh_m4=         -x10p12  -9.0_wp*x10m12
         dh_m5=                         +x10m12

         fltvar(i,j) = xp12*( var(i+1,j,m) - var(i  ,j,m) )     - &
                       xm12*( var(i  ,j,m) - var(i-1,j,m) )     + &
                        dh_0 *var(i  ,j,m)                      + &
                        dh_p1*var(i+1,j,m) + dh_m1*var(i-1,j,m) + &
                        dh_p2*var(i+2,j,m) + dh_m2*var(i-2,j,m) + &
                        dh_p3*var(i+3,j,m) + dh_m3*var(i-3,j,m) + &
                        dh_p4*var(i+4,j,m) + dh_m4*var(i-4,j,m) + &
                        dh_p5*var(i+5,j,m) + dh_m5*var(i-5,j,m)

      enddo
      enddo
      ! ------------------------------------------------------------------------
      if ( is_boundary(1) ) then
         ! left boundary
         do j=1,ny
            i=5
            lp12 = 0.5_wp*( rspec(i,j,1) + rspec(i+1,j,1) )/dx
            lm12 = 0.5_wp*( rspec(i,j,1) + rspec(i-1,j,1) )/dx

            xp12 = fltshock_amp*max( psens(i,j,1), psens(i+1,j,1) )
            xm12 = fltshock_amp*max( psens(i,j,1), psens(i-1,j,1) )
            x10p12 = lp12*max(0.0_wp,xnu10-coeff_corr10*xp12)
            x10m12 = lm12*max(0.0_wp,xnu8 -coeff_corr8 *xm12)
            xp12 = xp12*lp12
            xm12 = xm12*lm12

            dh_p5=          x10p12
            dh_p4=  -9.0_wp*x10p12
            dh_p3=  36.0_wp*x10p12        +x10m12
            dh_p2= -84.0_wp*x10p12 -7.0_wp*x10m12
            dh_p1= 126.0_wp*x10p12+21.0_wp*x10m12
            dh_0 =-126.0_wp*x10p12-35.0_wp*x10m12
            dh_m1=  84.0_wp*x10p12+35.0_wp*x10m12
            dh_m2= -36.0_wp*x10p12-21.0_wp*x10m12
            dh_m3=   9.0_wp*x10p12 +7.0_wp*x10m12
            dh_m4=         -x10p12        -x10m12

            fltvar(i,j) = xp12*( var(i+1,j,m) - var(i  ,j,m) )     - &
                            xm12*( var(i  ,j,m) - var(i-1,j,m) )     + &
                             dh_0 *var(i  ,j,m)                        + &
                             dh_p1*var(i+1,j,m) + dh_m1*var(i-1,j,m) + &
                             dh_p2*var(i+2,j,m) + dh_m2*var(i-2,j,m) + &
                             dh_p3*var(i+3,j,m) + dh_m3*var(i-3,j,m) + &
                             dh_p4*var(i+4,j,m) + dh_m4*var(i-4,j,m) + &
                             dh_p5*var(i+5,j,m)
            !----------------------------------------------------------------
            i=4
            lp12 = 0.5_wp*( rspec(i,j,1) + rspec(i+1,j,1) )/dx
            lm12 = 0.5_wp*( rspec(i,j,1) + rspec(i-1,j,1) )/dx

            xp12 = fltshock_amp*max( psens(i,j,1), psens(i+1,j,1) )
            xm12 = fltshock_amp*max( psens(i,j,1), psens(i-1,j,1) )
            x10p12 = lp12*max(0.0_wp,xnu8-coeff_corr8*xp12)
            x10m12 = lm12*max(0.0_wp,xnu6-coeff_corr6*xm12)
            xp12 = xp12*lp12
            xm12 = xm12*lm12

            dh_p4=        -x10p12
            dh_p3=  7.0_wp*x10p12
            dh_p2=-21.0_wp*x10p12        -x10m12
            dh_p1= 35.0_wp*x10p12 +5.0_wp*x10m12
            dh_0 =-35.0_wp*x10p12-10.0_wp*x10m12
            dh_m1= 21.0_wp*x10p12+10.0_wp*x10m12
            dh_m2= -7.0_wp*x10p12 -5.0_wp*x10m12
            dh_m3=         x10p12        +x10m12

            fltvar(i,j) = xp12*( var(i+1,j,m) - var(i  ,j,m) )     - &
                            xm12*( var(i  ,j,m) - var(i-1,j,m) )     + &
                             dh_0 *var(i  ,j,m)                        + &
                             dh_p1*var(i+1,j,m) + dh_m1*var(i-1,j,m) + &
                             dh_p2*var(i+2,j,m) + dh_m2*var(i-2,j,m) + &
                             dh_p3*var(i+3,j,m) + dh_m3*var(i-3,j,m) + &
                             dh_p4*var(i+4,j,m)
            !----------------------------------------------------------------
            i=3
            lp12 = 0.5_wp*( rspec(i,j,1) + rspec(i+1,j,1) )/dx
            lm12 = 0.5_wp*( rspec(i,j,1) + rspec(i-1,j,1) )/dx

            xp12 = fltshock_amp*max( psens(i,j,1), psens(i+1,j,1) )
            xm12 = fltshock_amp*max( psens(i,j,1), psens(i-1,j,1) )
            x10p12 = lp12*max(0.0_wp,xnu6-coeff_corr6*xp12)
            x10m12 = lm12*max(0.0_wp,xnu4-coeff_corr4*xm12)
            xp12 = xp12*lp12
            xm12 = xm12*lm12

            dh_p3=         x10p12
            dh_p2= -5.0_wp*x10p12
            dh_p1= 10.0_wp*x10p12       +x10m12
            dh_0 =-10.0_wp*x10p12-3.0_wp*x10m12
            dh_m1=  5.0_wp*x10p12+3.0_wp*x10m12
            dh_m2=        -x10p12       -x10m12

            fltvar(i,j) = xp12*( var(i+1,j,m) - var(i  ,j,m) ) - &
                            xm12*( var(i  ,j,m) - var(i-1,j,m) ) + &
                             dh_0 *var(i  ,j,m) + &
                             dh_p1*var(i+1,j,m) + dh_m1*var(i-1,j,m) + &
                             dh_p2*var(i+2,j,m) + dh_m2*var(i-2,j,m) + &
                             dh_p3*var(i+3,j,m)
            !----------------------------------------------------------------
            i=2
            lp12 = 0.5_wp*( rspec(i,j,1) + rspec(i+1,j,1) )/dx
            lm12 = 0.5_wp*( rspec(i,j,1) + rspec(i-1,j,1) )/dx

            xp12 = fltshock_amp*max( psens(i,j,1), psens(i+1,j,1) )
            xm12 = fltshock_amp*max( psens(i,j,1), psens(i-1,j,1) )
            x10p12 = lp12*max(0.0_wp,xnu4-coeff_corr4*xp12)
            x10m12 = lm12*max(0.0_wp,xnu2-coeff_corr4*xm12)
            xp12 = xp12*lp12
            xm12 = xm12*lm12

            dh_p2=       -x10p12
            dh_p1= 3.0_wp*x10p12
            dh_0 =-3.0_wp*x10p12-x10m12
            dh_m1=        x10p12+x10m12

            fltvar(i,j) = xp12*( var(i+1,j,m) - var(i  ,j,m) )     - &
                            xm12*( var(i  ,j,m) - var(i-1,j,m) )     + &
                             dh_0 *var(i  ,j,m)                        + &
                             dh_p1*var(i+1,j,m) + dh_m1*var(i-1,j,m) + &
                             dh_p2*var(i+2,j,m)
            !----------------------------------------------------------------
            i=1
            fltvar(i,j) = 0.0_wp
         enddo
         ! ---------------------------------------------------------------------
         ! right boundary
         do j=1,ny
            i=nx-4
            lp12 = 0.5_wp*( rspec(i,j,1) + rspec(i+1,j,1) )/dx
            lm12 = 0.5_wp*( rspec(i,j,1) + rspec(i-1,j,1) )/dx

            xp12 = fltshock_amp*max( psens(i,j,1), psens(i+1,j,1) )
            xm12 = fltshock_amp*max( psens(i,j,1), psens(i-1,j,1) )
            x10p12 = lp12*max(0.0_wp,xnu8 -coeff_corr8 *xp12)
            x10m12 = lm12*max(0.0_wp,xnu10-coeff_corr10*xm12)
            xp12 = xp12*lp12
            xm12 = xm12*lm12

            dh_p4=        -x10p12         -x10m12
            dh_p3=  7.0_wp*x10p12  +9.0_wp*x10m12
            dh_p2=-21.0_wp*x10p12 -36.0_wp*x10m12
            dh_p1= 35.0_wp*x10p12 +84.0_wp*x10m12
            dh_0 =-35.0_wp*x10p12-126.0_wp*x10m12
            dh_m1= 21.0_wp*x10p12+126.0_wp*x10m12
            dh_m2= -7.0_wp*x10p12 -84.0_wp*x10m12
            dh_m3=         x10p12 +36.0_wp*x10m12
            dh_m4=                 -9.0_wp*x10m12
            dh_m5=                        +x10m12

            fltvar(i,j) = xp12*( var(i+1,j,m) - var(i  ,j,m) )     - &
                            xm12*( var(i  ,j,m) - var(i-1,j,m) )     + &
                             dh_0 *var(i  ,j,m)                        + &
                             dh_p1*var(i+1,j,m) + dh_m1*var(i-1,j,m) + &
                             dh_p2*var(i+2,j,m) + dh_m2*var(i-2,j,m) + &
                             dh_p3*var(i+3,j,m) + dh_m3*var(i-3,j,m) + &
                             dh_p4*var(i+4,j,m) + dh_m4*var(i-4,j,m) + &
                                                    dh_m5*var(i-5,j,m)
            !----------------------------------------------------------------
            i=nx-3
            lp12 = 0.5_wp*( rspec(i,j,1) + rspec(i+1,j,1) )/dx
            lm12 = 0.5_wp*( rspec(i,j,1) + rspec(i-1,j,1) )/dx

            xp12 = fltshock_amp*max( psens(i,j,1), psens(i+1,j,1) )
            xm12 = fltshock_amp*max( psens(i,j,1), psens(i-1,j,1) )
            x10p12 = lp12*max(0.0_wp,xnu6-coeff_corr6*xp12)
            x10m12 = lm12*max(0.0_wp,xnu8-coeff_corr8*xm12)
            xp12 = xp12*lp12
            xm12 = xm12*lm12

            dh_p3=         x10p12        +x10m12
            dh_p2= -5.0_wp*x10p12 -7.0_wp*x10m12
            dh_p1= 10.0_wp*x10p12+21.0_wp*x10m12
            dh_0 =-10.0_wp*x10p12-35.0_wp*x10m12
            dh_m1=  5.0_wp*x10p12+35.0_wp*x10m12
            dh_m2=        -x10p12-21.0_wp*x10m12
            dh_m3=                +7.0_wp*x10m12
            dh_m4=                       -x10m12

            fltvar(i,j) = xp12*( var(i+1,j,m) - var(i  ,j,m) )     - &
                            xm12*( var(i  ,j,m) - var(i-1,j,m) )     + &
                             dh_0 *var(i  ,j,m)                        + &
                             dh_p1*var(i+1,j,m) + dh_m1*var(i-1,j,m) + &
                             dh_p2*var(i+2,j,m) + dh_m2*var(i-2,j,m) + &
                             dh_p3*var(i+3,j,m) + dh_m3*var(i-3,j,m) + &
                                                    dh_m4*var(i-4,j,m)
            !----------------------------------------------------------------
            i=nx-2
            lp12 = 0.5_wp*( rspec(i,j,1) + rspec(i+1,j,1) )/dx
            lm12 = 0.5_wp*( rspec(i,j,1) + rspec(i-1,j,1) )/dx

            xp12 = fltshock_amp*max( psens(i,j,1), psens(i+1,j,1) )
            xm12 = fltshock_amp*max( psens(i,j,1), psens(i-1,j,1) )
            x10p12 = lp12*max(0.0_wp,xnu4-coeff_corr4*xp12)
            x10m12 = lm12*max(0.0_wp,xnu6-coeff_corr6*xm12)
            xp12 = xp12*lp12
            xm12 = xm12*lm12

            dh_p2=       -x10p12        -x10m12
            dh_p1= 3.0_wp*x10p12 +5.0_wp*x10m12
            dh_0 =-3.0_wp*x10p12-10.0_wp*x10m12
            dh_m1=        x10p12+10.0_wp*x10m12
            dh_m2=               -5.0_wp*x10m12
            dh_m3=                      +x10m12

            fltvar(i,j) = xp12*( var(i+1,j,m) - var(i  ,j,m) )     - &
                          xm12*( var(i  ,j,m) - var(i-1,j,m) )     + &
                           dh_0 *var(i  ,j,m)                      + &
                           dh_p1*var(i+1,j,m) + dh_m1*var(i-1,j,m) + &
                           dh_p2*var(i+2,j,m) + dh_m2*var(i-2,j,m) + &
                                                dh_m3*var(i-3,j,m)
            !-------------------------------------------------------------------
            i=nx-1
            lp12 = 0.5_wp*( rspec(i,j,1) + rspec(i+1,j,1) )/dx
            lm12 = 0.5_wp*( rspec(i,j,1) + rspec(i-1,j,1) )/dx

            xp12 = fltshock_amp*max( psens(i,j,1), psens(i+1,j,1) )
            xm12 = fltshock_amp*max( psens(i,j,1), psens(i-1,j,1) )
            x10p12 = lp12*max(0.0_wp,xnu2-coeff_corr4*xp12)
            x10m12 = lm12*max(0.0_wp,xnu4-coeff_corr4*xm12)
            xp12 = xp12*lp12
            xm12 = xm12*lm12

            dh_p1= x10p12       +x10m12
            dh_0 =-x10p12-3.0_wp*x10m12
            dh_m1=       +3.0_wp*x10m12
            dh_m2=              -x10m12

            fltvar(i,j) = xp12*( var(i+1,j,m) - var(i  ,j,m) ) - &
                          xm12*( var(i  ,j,m) - var(i-1,j,m) ) + &
                           dh_0 *var(i  ,j,m) + &
                           dh_p1*var(i+1,j,m) + dh_m1*var(i-1,j,m) + &
                                                dh_m2*var(i-2,j,m)
            !-------------------------------------------------------------------
            i=nx
            fltvar(i,j) = 0.0_wp
         enddo
      endif
      ! ------------------------------------------------------------------------
      ! y direction
      ! -----------
      ! no boundaries use ghost cell
      if (ndim==2) then
         do j=jflt1,jflt2
         do i=1,nx
            lp12 = 0.5_wp*(rspec(i,j,2)+rspec(i,j+1,2))/dy
            lm12 = 0.5_wp*(rspec(i,j,2)+rspec(i,j-1,2))/dy

            xp12 = fltshock_amp*max(psens(i,j,2),psens(i,j+1,2))
            xm12 = fltshock_amp*max(psens(i,j,2),psens(i,j-1,2))
            x10p12 = lp12*max(0.0_wp,xnu10-coeff_corr10*xp12)
            x10m12 = lm12*max(0.0_wp,xnu10-coeff_corr10*xm12)
            xp12 = xp12*lp12
            xm12 = xm12*lm12

            dh_p5=          x10p12
            dh_p4=  -9.0_wp*x10p12         -x10m12
            dh_p3=  36.0_wp*x10p12  +9.0_wp*x10m12
            dh_p2= -84.0_wp*x10p12 -36.0_wp*x10m12
            dh_p1= 126.0_wp*x10p12 +84.0_wp*x10m12
            dh_0 =-126.0_wp*x10p12-126.0_wp*x10m12
            dh_m1=  84.0_wp*x10p12+126.0_wp*x10m12
            dh_m2= -36.0_wp*x10p12 -84.0_wp*x10m12
            dh_m3=   9.0_wp*x10p12 +36.0_wp*x10m12
            dh_m4=         -x10p12  -9.0_wp*x10m12
            dh_m5=                         +x10m12

            fltvar(i,j)= xp12*( var(i,j+1,m) - var(i,j  ,m) )     - &
                         xm12*( var(i,j  ,m) - var(i,j-1,m) )     + &
                          dh_0 *var(i,j  ,m)                      + &
                          dh_p1*var(i,j+1,m) + dh_m1*var(i,j-1,m) + &
                          dh_p2*var(i,j+2,m) + dh_m2*var(i,j-2,m) + &
                          dh_p3*var(i,j+3,m) + dh_m3*var(i,j-3,m) + &
                          dh_p4*var(i,j+4,m) + dh_m4*var(i,j-4,m) + &
                          dh_p5*var(i,j+5,m) + dh_m5*var(i,j-5,m) + fltvar(i,j)
         enddo
         enddo
         ! ---------------------------------------------------------------------
         if ( is_boundary(2) ) then
            ! bottom boundary
            do i=1,nx
               j=5
               lp12 = 0.5_wp*( rspec(i,j,2) + rspec(i,j+1,2) )/dy
               lm12 = 0.5_wp*( rspec(i,j,2) + rspec(i,j-1,2) )/dy

               xp12 = fltshock_amp*max( psens(i,j,2), psens(i,j+1,2) )
               xm12 = fltshock_amp*max( psens(i,j,2), psens(i,j-1,2) )
               x10p12 = lp12*max(0.0_wp,xnu10-coeff_corr10*xp12)
               x10m12 = lm12*max(0.0_wp,xnu8 -coeff_corr8 *xm12)
               xp12 = xp12*lp12
               xm12 = xm12*lm12

               dh_p5=          x10p12
               dh_p4=  -9.0_wp*x10p12
               dh_p3=  36.0_wp*x10p12        +x10m12
               dh_p2= -84.0_wp*x10p12 -7.0_wp*x10m12
               dh_p1= 126.0_wp*x10p12+21.0_wp*x10m12
               dh_0 =-126.0_wp*x10p12-35.0_wp*x10m12
               dh_m1=  84.0_wp*x10p12+35.0_wp*x10m12
               dh_m2= -36.0_wp*x10p12-21.0_wp*x10m12
               dh_m3=   9.0_wp*x10p12 +7.0_wp*x10m12
               dh_m4=         -x10p12        -x10m12

               fltvar(i,j) = xp12*( var(i,j+1,m) - var(i,j  ,m) )     - &
                               xm12*( var(i,j  ,m) - var(i,j-1,m) )     + &
                                dh_0 *var(i,j  ,m)                        + &
                                dh_p1*var(i,j+1,m) + dh_m1*var(i,j-1,m) + &
                                dh_p2*var(i,j+2,m) + dh_m2*var(i,j-2,m) + &
                                dh_p3*var(i,j+3,m) + dh_m3*var(i,j-3,m) + &
                                dh_p4*var(i,j+4,m) + dh_m4*var(i,j-4,m) + &
                                dh_p5*var(i,j+5,m) + fltvar(i,j)
               !----------------------------------------------------------------
               j=4
               lp12 = 0.5_wp*( rspec(i,j,2) + rspec(i,j+1,2) )/dy
               lm12 = 0.5_wp*( rspec(i,j,2) + rspec(i,j-1,2) )/dy

               xp12 = fltshock_amp*max( psens(i,j,2), psens(i,j+1,2) )
               xm12 = fltshock_amp*max( psens(i,j,2), psens(i,j-1,2) )
               x10p12 = lp12*max(0.0_wp,xnu8-coeff_corr8*xp12)
               x10m12 = lm12*max(0.0_wp,xnu6-coeff_corr6*xm12)
               xp12 = xp12*lp12
               xm12 = xm12*lm12

               dh_p4=        -x10p12
               dh_p3=  7.0_wp*x10p12
               dh_p2=-21.0_wp*x10p12        -x10m12
               dh_p1= 35.0_wp*x10p12 +5.0_wp*x10m12
               dh_0 =-35.0_wp*x10p12-10.0_wp*x10m12
               dh_m1= 21.0_wp*x10p12+10.0_wp*x10m12
               dh_m2= -7.0_wp*x10p12 -5.0_wp*x10m12
               dh_m3=         x10p12        +x10m12

               fltvar(i,j) = xp12*( var(i,j+1,m) - var(i,j  ,m) )     - &
                               xm12*( var(i,j  ,m) - var(i,j-1,m) )     + &
                                dh_0 *var(i,j  ,m)                        + &
                                dh_p1*var(i,j+1,m) + dh_m1*var(i,j-1,m) + &
                                dh_p2*var(i,j+2,m) + dh_m2*var(i,j-2,m) + &
                                dh_p3*var(i,j+3,m) + dh_m3*var(i,j-3,m) + &
                                dh_p4*var(i,j+4,m) + fltvar(i,j)
               !----------------------------------------------------------------
               j=3
               lp12 = 0.5_wp*( rspec(i,j,2) + rspec(i,j+1,2) )/dy
               lm12 = 0.5_wp*( rspec(i,j,2) + rspec(i,j-1,2) )/dy

               xp12 = fltshock_amp*max( psens(i,j,2), psens(i,j+1,2) )
               xm12 = fltshock_amp*max( psens(i,j,2), psens(i,j-1,2) )
               x10p12 = lp12*max(0.0_wp,xnu6-coeff_corr6*xp12)
               x10m12 = lm12*max(0.0_wp,xnu4-coeff_corr4*xm12)
               xp12 = xp12*lp12
               xm12 = xm12*lm12

               dh_p3=         x10p12
               dh_p2= -5.0_wp*x10p12
               dh_p1= 10.0_wp*x10p12       +x10m12
               dh_0 =-10.0_wp*x10p12-3.0_wp*x10m12
               dh_m1=  5.0_wp*x10p12+3.0_wp*x10m12
               dh_m2=        -x10p12       -x10m12

               fltvar(i,j) = xp12*( var(i,j+1,m) - var(i,j  ,m) )     - &
                               xm12*( var(i,j  ,m) - var(i,j-1,m) )     + &
                                dh_0 *var(i,j  ,m)                        + &
                                dh_p1*var(i,j+1,m) + dh_m1*var(i,j-1,m) + &
                                dh_p2*var(i,j+2,m) + dh_m2*var(i,j-2,m) + &
                                dh_p3*var(i,j+3,m) + fltvar(i,j)
               !----------------------------------------------------------------
               j=2
               lp12 = 0.5_wp*( rspec(i,j,2) + rspec(i,j+1,2) )/dy
               lm12 = 0.5_wp*( rspec(i,j,2) + rspec(i,j-1,2) )/dy

               xp12 = fltshock_amp*max( psens(i,j,2), psens(i,j+1,2) )
               xm12 = fltshock_amp*max( psens(i,j,2), psens(i,j-1,2) )
               x10p12 = lp12*max(0.0_wp,xnu4-coeff_corr4*xp12)
               x10m12 = lm12*max(0.0_wp,xnu2-coeff_corr4*xm12)
               xp12 = xp12*lp12
               xm12 = xm12*lm12

               dh_p2=       -x10p12
               dh_p1= 3.0_wp*x10p12
               dh_0 =-3.0_wp*x10p12-x10m12
               dh_m1=        x10p12+x10m12

               fltvar(i,j) = xp12*( var(i,j+1,m) - var(i,j  ,m) )     - &
                               xm12*( var(i,j  ,m) - var(i,j-1,m) )     + &
                                dh_0 *var(i,j  ,m)                        + &
                                dh_p1*var(i,j+1,m) + dh_m1*var(i,j-1,m) + &
                                dh_p2*var(i,j+2,m) + fltvar(i,j)
            enddo
         ! ---------------------------------------------------------------------
            ! top boundary
            do i=1,nx
               j=ny-4
               lp12 = 0.5_wp*( rspec(i,j,2) + rspec(i,j+1,2) )/dy
               lm12 = 0.5_wp*( rspec(i,j,2) + rspec(i,j-1,2) )/dy

               xp12 = fltshock_amp*max( psens(i,j,2), psens(i,j+1,2) )
               xm12 = fltshock_amp*max( psens(i,j,2), psens(i,j-1,2) )
               x10p12 = lp12*max(0.0_wp,xnu8 -coeff_corr8 *xp12)
               x10m12 = lm12*max(0.0_wp,xnu10-coeff_corr10*xm12)
               xp12 = xp12*lp12
               xm12 = xm12*lm12

               dh_p4=        -x10p12         -x10m12
               dh_p3=  7.0_wp*x10p12  +9.0_wp*x10m12
               dh_p2=-21.0_wp*x10p12 -36.0_wp*x10m12
               dh_p1= 35.0_wp*x10p12 +84.0_wp*x10m12
               dh_0 =-35.0_wp*x10p12-126.0_wp*x10m12
               dh_m1= 21.0_wp*x10p12+126.0_wp*x10m12
               dh_m2= -7.0_wp*x10p12 -84.0_wp*x10m12
               dh_m3=         x10p12 +36.0_wp*x10m12
               dh_m4=                 -9.0_wp*x10m12
               dh_m5=                        +x10m12

               fltvar(i,j) = xp12*( var(i,j+1,m) - var(i,j  ,m) )     - &
                               xm12*( var(i,j  ,m) - var(i,j-1,m) )     + &
                                dh_0 *var(i,j  ,m)                        + &
                                dh_p1*var(i,j+1,m) + dh_m1*var(i,j-1,m) + &
                                dh_p2*var(i,j+2,m) + dh_m2*var(i,j-2,m) + &
                                dh_p3*var(i,j+3,m) + dh_m3*var(i,j-3,m) + &
                                dh_p4*var(i,j+4,m) + dh_m4*var(i,j-4,m) + &
                                                       dh_m5*var(i,j-5,m) + fltvar(i,j)
               !----------------------------------------------------------------
               j=ny-3
               lp12 = 0.5_wp*( rspec(i,j,2) + rspec(i,j+1,2) )/dy
               lm12 = 0.5_wp*( rspec(i,j,2) + rspec(i,j-1,2) )/dy

               xp12 = fltshock_amp*max( psens(i,j,2), psens(i,j+1,2) )
               xm12 = fltshock_amp*max( psens(i,j,2), psens(i,j-1,2) )
               x10p12 = lp12*max(0.0_wp,xnu6-coeff_corr6*xp12)
               x10m12 = lm12*max(0.0_wp,xnu8-coeff_corr8*xm12)
               xp12 = xp12*lp12
               xm12 = xm12*lm12

               dh_p3=         x10p12        +x10m12
               dh_p2= -5.0_wp*x10p12 -7.0_wp*x10m12
               dh_p1= 10.0_wp*x10p12+21.0_wp*x10m12
               dh_0 =-10.0_wp*x10p12-35.0_wp*x10m12
               dh_m1=  5.0_wp*x10p12+35.0_wp*x10m12
               dh_m2=        -x10p12-21.0_wp*x10m12
               dh_m3=                +7.0_wp*x10m12
               dh_m4=                       -x10m12

               fltvar(i,j) = xp12*( var(i,j+1,m) - var(i,j  ,m) )     - &
                               xm12*( var(i,j  ,m) - var(i,j-1,m) )     + &
                                dh_0 *var(i,j  ,m)                        + &
                                dh_p1*var(i,j+1,m) + dh_m1*var(i,j-1,m) + &
                                dh_p2*var(i,j+2,m) + dh_m2*var(i,j-2,m) + &
                                dh_p3*var(i,j+3,m) + dh_m3*var(i,j-3,m) + &
                                                       dh_m4*var(i,j-4,m) + fltvar(i,j)
               !----------------------------------------------------------------
               j=ny-2
               lp12 = 0.5_wp*( rspec(i,j,2) + rspec(i,j+1,2) )/dy
               lm12 = 0.5_wp*( rspec(i,j,2) + rspec(i,j-1,2) )/dy

               xp12 = fltshock_amp*max( psens(i,j,2), psens(i,j+1,2) )
               xm12 = fltshock_amp*max( psens(i,j,2), psens(i,j-1,2) )
               x10p12 = lp12*max(0.0_wp,xnu4-coeff_corr4*xp12)
               x10m12 = lm12*max(0.0_wp,xnu6-coeff_corr6*xm12)
               xp12 = xp12*lp12
               xm12 = xm12*lm12

               dh_p2=       -x10p12        -x10m12
               dh_p1= 3.0_wp*x10p12 +5.0_wp*x10m12
               dh_0 =-3.0_wp*x10p12-10.0_wp*x10m12
               dh_m1=        x10p12+10.0_wp*x10m12
               dh_m2=               -5.0_wp*x10m12
               dh_m3=                      +x10m12

               fltvar(i,j) = xp12*( var(i,j+1,m) - var(i,j  ,m) )     - &
                               xm12*( var(i,j  ,m) - var(i,j-1,m) )     + &
                                dh_0 *var(i,j  ,m)                        + &
                                dh_p1*var(i,j+1,m) + dh_m1*var(i,j-1,m) + &
                                dh_p2*var(i,j+2,m) + dh_m2*var(i,j-2,m) + &
                                                       dh_m3*var(i,j-3,m) + fltvar(i,j)
               !----------------------------------------------------------------
               j=ny-1
               lp12 = 0.5_wp*( rspec(i,j,2) + rspec(i,j+1,2) )/dy
               lm12 = 0.5_wp*( rspec(i,j,2) + rspec(i,j-1,2) )/dy

               xp12 = fltshock_amp*max( psens(i,j,2), psens(i,j+1,2) )
               xm12 = fltshock_amp*max( psens(i,j,2), psens(i,j-1,2) )
               x10p12 = lp12*max(0.0_wp,xnu2-coeff_corr4*xp12)
               x10m12 = lm12*max(0.0_wp,xnu4-coeff_corr4*xm12)
               xp12 = xp12*lp12
               xm12 = xm12*lm12

               dh_p1= x10p12       +x10m12
               dh_0 =-x10p12-3.0_wp*x10m12
               dh_m1=       +3.0_wp*x10m12
               dh_m2=              -x10m12

               fltvar(i,j) = xp12*( var(i,j+1,m) - var(i,j  ,m) )     - &
                               xm12*( var(i,j  ,m) - var(i,j-1,m) )     + &
                                dh_0 *var(i,j  ,m)                        + &
                                dh_p1*var(i,j+1,m) + dh_m1*var(i,j-1,m) + &
                                                       dh_m2*var(i,j-2,m) + fltvar(i,j)
            enddo
         endif
      endif
      var(1:nx,1:ny,m) = var(1:nx,1:ny,m) + deltat*fltvar(1:nx,1:ny)
   enddo Varsloop

end subroutine filter_shock_jameson