subroutine runge_kutta4
!------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Perform time integration with classical 4-step Runge-Kutta
!-------------------------------------------------------------------------------
   use mod_mode
   use mod_time
   use mod_work
   implicit none
   ! ---------------------------------------------------------------------------
   ! Local variables
   integer :: m, i, j
   real(8) :: g1, g2, g3, g4
   ! ---------------------------------------------------------------------------

   g1=1.0_wp; g2=0.5_wp; g3=1.0_wp/6.0_wp; g4 = 1.0_wp/24.0_wp

   coeff_rk(4,1) = g1                                ! 1
   coeff_rk(3,1) = g2                                ! 1/2
   coeff_rk(2,1) = g3/ coeff_rk(3,1)                 ! 1/3
   coeff_rk(1,1) = g4/(coeff_rk(3,1)*coeff_rk(2,1))  ! 1/4

   coeff_rk = deltat*coeff_rk

   ! integrate governing equations
   ! =============================
   ! save current solution
   phitemp(:,:,:) = phi(:,:,:)
   rungeloop: do istage = 1, nrk

      if ( is_var_not_updated ) call update_var( phitemp )

      !-------------------------------------------------------------------------
      ! calculate inviscid and viscous fluxes
      ! -------------------------------------
      call rhs_inviscid( phitemp, dphidt )

      if ( is_visc .or. fltshock_type=='L' ) call rhs_viscous( phitemp, dphidt )

      if (ndim==1) dphidt(:,:,3) = 0.0_wp
      ! update final solution
      do m=1,nvars
         do j=1,ny
         do i=1,nx
            phitemp(i,j,m) = phi(i,j,m) + coeff_rk(istage,1)*dphidt(i,j,m)
         enddo
         enddo
      enddo

      is_var_not_updated = .true.
   enddo rungeloop

   phi(1:nx,1:ny,1:nvars) = phitemp(1:nx,1:ny,1:nvars)
end subroutine runge_kutta4