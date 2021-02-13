subroutine setup_init_cond
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Set up intial conditions for various problems
!-------------------------------------------------------------------------------
   use mod_grid
   use mod_fluid
   use mod_work
   implicit none
   !----------------------------------------------------------------------------
   ! Local variables
   real(wp), dimension(nx1:nx2,ny1:ny2,2) :: U_init
   integer  :: i, j
   real(wp) :: rho_left , p_left , u_left , T_left
   real(wp) :: rho_right, p_right, u_right, T_right
   real(wp) :: xloc, yloc, rloc, xc_vort, yc_vort
   real(wp) :: u, v, rho, p, t
   !----------------------------------------------------------------------------
   ! Variables for vortex convection init
   real(wp) :: beta, sigma, floc, dT, Omega
   !----------------------------------------------------------------------------
   ! Variables for shock tube init
   real(wp) :: x_dia
   !----------------------------------------------------------------------------
   ! Variables for shock vortex interaction init
   real(wp) :: theta, Mach_v, kn
   !----------------------------------------------------------------------------

   write(*,*) 'Setting up initial conditions'

   select case ( flowtype )
   !----------------------------------------------------------------------------
   case ( SHOCK_TUBE )

      Pinf   = 1e5_wp
      Rhoinf = 1.0_wp

      Uscale = sqrt(Pinf/rhoinf)
      Lscale = 1.0_wp
      Tscale = Lscale/Uscale

      rho_left  = rhoinf                   ! 1     [kg/m^3]
      rho_right = rhoinf  /8.0_wp          ! 0.125 [kg/m^3]
      p_left    = Pinf                     ! 1     [Pa]
      p_right   = Pinf    /10.0_wp         ! 0.1   [Pa]
      u_left    = Uscale  *0.0_wp          ! 0     [m/s]
      u_right   = Uscale  *0.0_wp          ! 0     [m/s]
      T_left    = p_left /(Rair*rho_left)  ! [K]
      T_right   = p_right/(Rair*rho_right) ! [K]

      x_dia = xmin + Xlength/2.0_wp
      do j = 1, ny
      do i = 1, nx
         if ( gridx(i)  <= x_dia) then
            prs(i,j)      = p_left
            tmp(i,j)      = T_left
            U_init(i,j,1) = u_left
            U_init(i,j,2) = 0.0_wp
         else
            prs(i,j)      = p_right
            tmp(i,j)      = T_right
            U_init(i,j,1) = 0.0_wp
            U_init(i,j,2) = 0.0_wp
         endif
      enddo
      enddo

   !----------------------------------------------------------------------------
   case( SHU_OSHER )

      Pinf   = 1e5_wp
      Rhoinf = 1.0_wp

      Uscale = sqrt(Pinf/rhoinf)
      Lscale = 1.0_wp
      Tscale = Lscale/Uscale
      do j=1,ny
      do i=1,nx
         if ( gridx(i) < -4.0_wp*Lscale ) then
            rho           = Rhoinf*3.857143_wp
            prs(i,j)      = Pinf*10.33333_wp
            tmp(i,j)      = prs(i,j) / (Rair*rho)
            U_init(i,j,1) = Uscale*2.629369_wp
            U_init(i,j,2) = 0.0_wp
         else
            rho           = Rhoinf*( 1.0_wp + 0.2_wp*sin(5.0_wp*gridx(i)) )
            prs(i,j)      = Pinf
            tmp(i,j)      = prs(i,j) / (Rair*rho)
            U_init(i,j,1) = 0.0_wp
            U_init(i,j,2) = 0.0_wp
         endif
      enddo
      enddo

   !----------------------------------------------------------------------------
   case ( VORTEX_CONVECTION )

      Tinf = 300.0_wp
      Pinf = 1e5_wp
      Ainf = sqrt(gamma*Rair*Tinf)
      Minf = sqrt(2.0_wp/gamma)

      alpha_dir = 45.0_wp*(pi/180.0_wp)
      Uscale = Minf*Ainf
      Lscale = 1.0_wp
      Tscale = Lscale*(cos(alpha_dir)+sin(alpha_dir))/Uscale

      Rloc  = 1.0_wp
      sigma = 1.0_wp
      beta  = Minf*5.0_wp*sqrt(2.0_wp)/(4.0_wp*pi)*exp(0.5_wp)

      xc_vort = Xlength*0.5_wp
      yc_vort = Ylength*0.5_wp

      do j=1,ny
      do i=1,nx
         ! vortex location
         xloc = (gridx(i) - xc_vort)/Rloc
         yloc = (gridy(j) - yc_vort)/Rloc
         floc = -1.0_wp/(2.0_wp*sigma**2)*( xloc**2 + yloc**2 )
         Omega = beta*exp(floc)
         ! velocity perturbation
         u = -yloc*Omega
         v =  xloc*Omega
         U_init(i,j,1) = ( Minf*cos(alpha_dir) + u )*Ainf
         U_init(i,j,2) = ( Minf*sin(alpha_dir) + v )*Ainf
         ! Temperature + temperature perturbation
         dT = - (gamma-1.0_wp)/2.0_wp*Omega**2
         Tmp(i,j) = Tinf*(1.0_wp + dT)
         ! pressure + pressure perturbation
         prs(i,j) = Pinf*(1.0_wp + dT)**( gamma/(gamma-1.0_wp) )
      enddo
      enddo

   !----------------------------------------------------------------------------
   case ( SHOCK_VORTEX_INT )

      Tinf = 300.0_wp
      Pinf = 1e5_wp
      Minf = 1.2_wp
      Ainf = sqrt(gamma*Rair*Tinf)

      Mach_v = 0.25_wp ! Vortex Mach number

      Uscale = Minf*Ainf
      Lscale = 1.0_wp
      Tscale = Lscale/Uscale !(Ainf*(Minf-Mach_v))

      ! Upstream mach number Ms is Mach0
      T_left   = Tinf
      p_left   = Pinf
      u_left   = Uscale
      rho_left = p_left/(Rair*T_left)

      ! Postshock state found by classical Rankine-Hugoniot relations for perfect gases
      p_right   = p_left*( 1.0_wp + 2.0_wp*gamma/(gamma+1.0_wp)*(Minf**2 - 1.0_wp) )
      rho_right = rho_left*(gamma+1.0_wp)*Minf**2 / ( 2.0_wp + (gamma-1.0_wp)*Minf**2 )
      u_right   = u_left*rho_left/rho_right
      T_right   = p_right/(rho_right*Rair)

      ! Lscale is the Radius of the vortex
      xc_vort = -20.0_wp*Lscale
      yc_vort =   0.0_wp*Lscale
      do j=1,ny
      do i=1,nx
         ! vortex location
         xloc = gridx(i) - xc_vort
         yloc = gridy(j) - yc_vort
         rloc = sqrt( xloc**2 + yloc**2 )
         theta = atan2(yloc, xloc)
         if (gridx(i) >0.0_wp) then
            U_init(i,j,1) = u_right
            U_init(i,j,2) = 0.0_wp
            prs(i,j) = p_right
            tmp(i,j) = T_right
         else
            ! velocity perturbation
            u = -Mach_v*rloc*exp(0.5_wp*(1.0_wp-rloc**2))*sin(theta)*Ainf
            v =  Mach_v*rloc*exp(0.5_wp*(1.0_wp-rloc**2))*cos(theta)*Ainf
            U_init(i,j,1) = u_left + u
            U_init(i,j,2) = v
            ! density + density perturbation
            rho = rho_left*(1.0_wp - (gamma-1.0_wp)/2.0_wp*Mach_v**2*exp(1-rloc**2))**(1.0_wp/(gamma-1.0_wp))
            ! pressure, no perturbation
            prs(i,j) = p_left/rho_left**gamma*rho**gamma
            ! temperature + perturbation
            tmp(i,j) = prs(i,j)/(Rair*rho)
         endif
      enddo
      enddo

   !----------------------------------------------------------------------------
   case ( TAYLOR_GREEN_VORTEX )

      Pinf   = 1e5_wp
      Rhoinf = 1.0_wp
      Minf   = 0.3_wp
      Tinf   = Pinf/(Rair*Rhoinf)
      Ainf   = sqrt(gamma*Rair*Tinf)

      Uscale = Minf*Ainf
      Lscale = 1.0_wp
      Tscale = Lscale/Uscale

      mu = Rhoinf*Uscale*Lscale/2000.0_wp
      lambda = mu*cp/Pr

      kn = 4.0_wp
      do j = 1,ny
      do i = 1,nx
         xloc = gridx(i)
         yloc = gridy(j)
         ! velocity
         U_init(i,j,1) =  Uscale*sin(kn*xloc)*cos(kn*yloc)
         U_init(i,j,2) = -Uscale*cos(kn*xloc)*sin(kn*yloc)
         ! pressure (assume ideal gas for initial condition)
         prs(i,j) = Pinf + Rhoinf*Uscale**2 &
                  * ( cos(2.0_wp*kn*xloc) + cos(2.0_wp*kn*yloc) )/4.0_wp
         ! temperature (assume ideal gas for initial condition)
         tmp(i,j) = prs(i,j) / (Rair*Rhoinf)
      enddo
      enddo

   !----------------------------------------------------------------------------
   case default
      write(*,*) 'Case not known in setup_init_cond.f90. ABORT'
      stop
   end select
   !----------------------------------------------------------------------------

   ! ---------------------------------------------------------------------------
   ! Now fill the state vector
   ! Compute density,internal energy using equation of state
   do j = 1,ny
   do i = 1,nx

      u = U_init(i,j,1)
      v = U_init(i,j,2)

      P = prs(i,j)
      T = tmp(i,j)

      rho = P/(Rair*T)

      phi(i,j,1) = rho
      phi(i,j,2) = rho*u
      phi(i,j,3) = rho*v

      phi(i,j,4) = rho*( cv*T + 0.5_wp*(u**2+v**2) )

      as(i,j) = sqrt(gamma*Rair*T)
   enddo
   enddo

end subroutine setup_init_cond