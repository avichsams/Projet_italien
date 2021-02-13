subroutine compute_lad( phivar )
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Compute Localized Artificial Diffusivity.
!  Standard values for the coefficients are c_k  = 0.01, c_beta = 1.5
!-------------------------------------------------------------------------------
   use mod_mode
   use mod_grid
   use mod_work
   use mod_deriv
   use mod_filter
   implicit none
   !----------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2,nvars), intent(in) :: phivar
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i, j
   real(wp), parameter :: cutoff = 1.e-30_wp
   real(wp) :: fsw, b3, b4
   real(wp) :: gradx , grady , norm_i
   real(wp) :: gradx2, grady2, norm2_i
   !----------------------------------------------------------------------------

   ! -------------------------
   ! Artificial bulk viscosity
   ! -------------------------
   call update_ghost( div, 0 )
   call deriv4th_x_exp(div, dsxdx)
   call deriv4th_y_exp(div, dsydy)

   call derivative_x( phivar(:,:,1), drho(:,:,1), deriv_visc_order )
   call derivative_y( phivar(:,:,1), drho(:,:,2), deriv_visc_order )

   do j=1,ny
   do i=1,nx

      if ( div(i,j) >=0.0_wp ) then
         wrk(i,j) = 0.0_wp
      else
         b3 = dvel(i,j,2,1) - dvel(i,j,1,2)
         b4 = div(i,j)**2
         fsw = b4 / (b4 + b3**2 + cutoff)

         gradx2 = drho(i,j,1)**2
         grady2 = drho(i,j,2)**2

         norm2_i = 1.0_wp/(gradx2 + grady2 + cutoff)

         wrk(i,j) = c_beta*phivar(i,j,1)*fsw      &
                  * abs( dsxdx(i,j)*gradx2*dx**2 + dsydy(i,j)*grady2*dy**2 )*norm2_i
      endif
   enddo
   enddo

   call update_ghost( wrk, 0 )
   call filter4_gaussian( wrk, bulk_lad )

   !----------------------------------------------------------------------------
   ! -------------------------------
   ! Artificial thermal conductivity
   ! -------------------------------
   if (c_k >small) then
      ! Internal energy, no need to pass ghostcells
      do j=1,ny
      do i=1,nx
         wrk(i,j) = (phivar(i,j,4) - 0.5_wp*( phivar(i,j,2)**2 &
                                            + phivar(i,j,3)**2 )/phivar(i,j,1)) / phivar(i,j,1)
      enddo
      enddo
      ! For grad(e) norm, the qflux array is reused!
      call derivative_x( wrk, qflux(:,:,1), deriv_visc_order)
      call derivative_y( wrk, qflux(:,:,2), deriv_visc_order)

      call update_ghost( wrk, 0 )
      call deriv4th_x_exp( wrk, dsxdx )
      call deriv4th_y_exp( wrk, dsydy )

      do j=1,ny
      do i=1,nx
         gradx = qflux(i,j,1)
         grady = qflux(i,j,2)
         norm_i = 1.0_wp/sqrt(gradx**2 + grady**2 + cutoff)

         wrk(i,j) = c_k*phivar(i,j,1)*as(i,j)/Tmp(i,j) &
                  * abs( dsxdx(i,j)*gradx*dx + dsydy(i,j)*grady*dy )*norm_i
      enddo
      enddo

      call update_ghost( wrk, 0 )
      call filter4_gaussian( wrk, lambda_lad )
   endif

   contains

   subroutine deriv4th_x_exp( var, dvar )
   !----------------------------------------------------------------------------
   !  AUTHOR
   !  ------
   !  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
   !
   !  DESCRIPTION
   !  -----------
   !  Compute 4th derivative using explicit fourth order central difference
   !  scheme - Kawai 2010 JCP
   !----------------------------------------------------------------------------
      use mod_deriv
      implicit none
      ! ------------------------------------------------------------------------
      ! Input/Output arguments
      real(wp), dimension(nx1:nx2,ny1:ny2), intent(in)  :: var
      real(wp), dimension(  1:nx ,  1:ny ), intent(out) :: dvar
      ! ------------------------------------------------------------------------
      ! Local variables
      real(wp), parameter :: c0 = 28.0_wp/3.0_wp, c1 = -13.0_wp/2.0_wp, &
                             c2 = 2.0_wp        , c3 = -1.0_wp /6.0_wp
      integer :: i, j
      ! ------------------------------------------------------------------------

      dvar(:,:) = 0.0_wp
      ! X-direction, central differencing
      do j = 1,ny
      do i = ider1_06,ider2_06
         dvar(i,j) =  c0* var(i  ,j)                &
                   +  c1*(var(i+1,j) + var(i-1,j))  &
                   +  c2*(var(i+2,j) + var(i-2,j))  &
                   +  c3*(var(i+3,j) + var(i-3,j))
      enddo
      enddo

   end subroutine deriv4th_x_exp

   !----------------------------------------------------------------------------

   subroutine deriv4th_y_exp( var, dvar )
   !----------------------------------------------------------------------------
   !  AUTHOR
   !  ------
   !  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
   !
   !  DESCRIPTION
   !  -----------
   !  Compute 4th derivative using explicit fourth order central difference
   !  scheme - Kawai 2010 JCP
   !----------------------------------------------------------------------------
      use mod_deriv
      implicit none
      ! ------------------------------------------------------------------------
      ! Input/Output arguments
      real(wp), dimension(nx1:nx2,ny1:ny2), intent(in)  :: var
      real(wp), dimension(  1:nx ,  1:ny ), intent(out) :: dvar
      ! ------------------------------------------------------------------------
      ! Local variables
      real(wp), parameter :: c0 = 28.0_wp/3.0_wp, c1 = -13.0_wp/2.0_wp, &
                             c2 = 2.0_wp        , c3 = -1.0_wp /6.0_wp
      integer :: i, j
      ! ------------------------------------------------------------------------

      dvar(:,:) = 0.0_wp
      if (ndim==1) return
      ! Y-direction, central differencing
      do j = jder1_06,jder2_06
      do i = 1,nx
         dvar(i,j) =  c0* var(i,j  )                &
                   +  c1*(var(i,j+1) + var(i,j-1))  &
                   +  c2*(var(i,j+2) + var(i,j-2))  &
                   +  c3*(var(i,j+3) + var(i,j-3))
      enddo
      enddo

   end subroutine deriv4th_y_exp

   !----------------------------------------------------------------------------
   subroutine filter4_gaussian( var, varf )
   !----------------------------------------------------------------------------
   !  AUTHOR
   !  ------
   !  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
   !
   !  DESCRIPTION
   !  -----------
   !  Truncated Gaussian filter, Kawai 2010 JCP
   !  For non-periodic boundaries, var is mirrored across the boundary (Kawai 2008)
   !----------------------------------------------------------------------------
      use mod_filter
      implicit none
      ! ------------------------------------------------------------------------
      ! Input/Output arguments
      real(wp), dimension(nx1:nx2,ny1:ny2), intent(in)  :: var
      real(wp), dimension(  1:nx ,  1:ny ), intent(out) :: varf
      ! ------------------------------------------------------------------------
      ! Local variables
      integer :: i, j
      real(wp), parameter :: a1 = 3565.0_wp/ 10368.0_wp, &
                             a2 = 3091.0_wp/ 12960.0_wp, &
                             a3 = 1997.0_wp/ 25920.0_wp, &
                             a4 =  149.0_wp/ 12960.0_wp, &
                             a5 =  107.0_wp/103680.0_wp
      ! ------------------------------------------------------------------------

      varf = 0.0_wp
      ! X-direction, central filtering
      do j = 1,ny
      do i = igaus1,igaus2
         varf(i,j) = a1*var(i,j) &
                   + a2*(var(i+1,j) + var(i-1,j)) &
                   + a3*(var(i+2,j) + var(i-2,j)) &
                   + a4*(var(i+3,j) + var(i-3,j)) &
                   + a5*(var(i+4,j) + var(i-4,j))
      enddo
      enddo

      if (ndim==1) return
      ! ------------------------------------------------------------------------
      ! Y-direction, central filtering
      do j = jgaus1,jgaus2
      do i = 1,nx
         varf(i,j) = a1*var(i,j) &
                     + a2*(var(i,j+1) + var(i,j-1)) &
                     + a3*(var(i,j+2) + var(i,j-2)) &
                     + a4*(var(i,j+3) + var(i,j-3)) &
                     + a5*(var(i,j+4) + var(i,j-4)) + varf(i,j)
      enddo
      enddo

   end subroutine filter4_gaussian
end subroutine compute_lad