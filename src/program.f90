program TP
! ------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Program for testing high-order methods in 1D and 2D canonical configurations.
! ------------------------------------------------------------------------------
   use mod_mode
   use mod_grid
   use mod_filter
   use mod_fluid
   use mod_deriv
   use mod_time
   use mod_work
   implicit none
   ! ---------------------------------------------------------------------------
   ! Local variables
   integer, parameter :: io=30
   real(wp) :: time_after, time_before
   ! ---------------------------------------------------------------------------
   time = 0.0_wp
   tstar = 0.0_wp
   cputime_tot = 0.0_wp
   nprint = 1

   open(UNIT=io, FILE='param.inp', STATUS='old', ACTION='read')
   read(io,*) flowtype
   read(io,*) nx, Xlength
   read(io,*) ny, Ylength
   read(io,*) is_visc
   read(io,*) timemax
   read(io,*) dtprint
   read(io,*) printstdout
   read(io,*) CFL, Fo
   read(io,*) deriv_conv_type, deriv_visc_type
   read(io,*) deriv_conv_order, deriv_visc_order
   read(io,*) filter_type
   read(io,*) nrk
   read(io,*) filter_order
   read(io,*) fltamp
   read(io,*) is_shock
   if ( is_shock ) then
      read(io,*) fltshock_type
      read(io,*) fltshock_amp
      read(io,*) c_beta, c_k
   endif
   close(UNIT=io)
   ! ---------------------------------------------------------------------------
   nsteps  = 999999
   ! ---------------------------------------------------------------------------
   call init_mode
   call init_work
   call init_grid
   call init_time
   call init_deriv
   call init_filter
   call init_fluid

   call setup_init_cond
   call update_var( phi )
   call write_solution( phi )

   call deltat_fromCFL
   ! ---------------------------------------------------------------------------
   ! START TIME INTEGRATION LOOP
   ! ---------------------------------------------------------------------------
   print*,'Rung Kutta',nrk
   looptime: do ntime = 1,nsteps

      call cpu_time(time_before)

      ! Evaluate time step size
      if ( fltshock_type=='L' ) call compute_lad( phi )
      call deltat_fromCFL

      ! Adjust timestep to print field at the desired frequency
      if ( (time+deltat)/tscale > nprint*dtprint ) then
         deltat = tscale*nprint*dtprint - time
         printstep = .true.
         nprint = nprint+1
      endif

      time   = time + deltat
      tstar  = time / tscale
      dtstar = deltat / tscale

      if ( mod(ntime,printstdout)==0 ) then
         write(*,'(1x, a, 2(i0, a), 1pE15.7)') &
              'Iteration: ', ntime, '/', ntime+int((timemax-tstar)/dtstar), ', Tstar: ', tstar
         write(*,"(1x,'dt*, dt: ', 2(1X,1pe12.5))") dtstar, deltat
      endif

      ! ------------------------------------------------------------------------
      ! Runge-Kutta integration
      if ( nrk==4) then
         call runge_kutta4
      elseif (nrk==2) then 
         call runge_kutta2
      elseif (nrk==6) then 
         call runge_kutta6
      endif
      
      ! ------------------------------------------------------------------------
      ! Filtering
      if ( is_shock ) then
         if ( fltshock_type=='J' ) then
            call filter_shock_jameson( phi )
         elseif ( fltshock_type=='L' ) then
            call filter_exp( phi )
         elseif ( fltshock_type=='B' ) then
            call filter_exp( phi )
            call filter_shock_bogey( phi )
         endif
      else
         call filter_exp( phi )
      endif
      ! ------------------------------------------------------------------------
      ! Check if NaN is issued
      if (isnan( sum(phi(:,:,1)) ) ) then
         write(*,*) 'NaN detected at it', ntime
         stop
      endif
      ! ------------------------------------------------------------------------
      ! Handle Output
      if ( printstep ) then
         call write_solution( phi )
         printstep = .false.
      endif
      ! ------------------------------------------------------------------------
      ! Determine elapsed time
      call cpu_time(time_after)
      cputime_ite = time_after - time_before
      cputime_tot = cputime_tot + cputime_ite
      if ( mod(ntime,printstdout)==0 ) then
         write(*,20) cputime_ite, int(cputime_tot)
         write(*,'(a)') repeat("-", 80)
      endif
 20 format('CPU time/it:',f12.5,' s, Run time: ',i0,' s')
      ! ------------------------------------------------------------------------
      ! Handle exit conditions
      if ( tstar >=timemax ) then
         write(*,*) 'Max Tstar reached = ', timemax
         exit looptime
      endif
   enddo looptime
   ! ---------------------------------------------------------------------------
   ! END TIME INTEGRATION LOOP
   ! ---------------------------------------------------------------------------

   write(*,*) 'The end!'
end program TP