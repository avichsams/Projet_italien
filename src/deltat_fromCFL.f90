subroutine deltat_fromCFL
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  evaluate timestep from CFL
!-------------------------------------------------------------------------------
   use mod_grid, only: dx, dy
   use mod_mode
   use mod_work
   use mod_time
   use mod_fluid, only: mu
   implicit none
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i, j
   real(wp) :: testc, testv, deltat_c, deltat_v
   !----------------------------------------------------------------------------

   deltat_c = 9.e+99_wp
   deltat_v = 9.e+99_wp

   ! ---------------------------------------------------------------------------

   if (is_visc .or. fltshock_type=='L') then
      ! Convective and viscous conditions
      do j = 1,ny
      do i = 1,nx
         testc = CFL*min( dx/( abs(vel(i,j,1)) + as(i,j) ) &
                        , dy/( abs(vel(i,j,2)) + as(i,j) ) )
         deltat_c = min(deltat_c, testc)

         testv = Fo*phi(i,j,1)*min(dx,dy)**2 / max((mu + bulk_lad(i,j)), 1e-10_wp)
         deltat_v = min(deltat_v, testv)
      enddo
      enddo
   else
      ! Convective condition only
      do j = 1,ny
      do i = 1,nx
         testc = CFL*min( dx/( abs(vel(i,j,1)) + as(i,j) ) &
                        , dy/( abs(vel(i,j,2)) + as(i,j) ) )
         deltat_c = min(deltat_c, testc)
      enddo
      enddo
   endif

   deltat = min(deltat_c, deltat_v)

end subroutine deltat_fromCFL