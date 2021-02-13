subroutine compute_thermo( phivar )
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Compute thermodynamic properties from state vector
!-------------------------------------------------------------------------------
   use mod_work
   use mod_fluid
   implicit none
   ! ---------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2,4), intent(in) :: phivar
   ! Local Variables
   integer  :: i, j
   real(wp) :: T, etot, u1, u2, rom1, rho
   ! ---------------------------------------------------------------------------

   do j = 1,ny
   do i = 1,nx
      rho  = phivar(i,j,1)
      rom1 = 1.0_wp/rho
      u1   = phivar(i,j,2)*rom1
      u2   = phivar(i,j,3)*rom1
      etot = phivar(i,j,4)*rom1

      T = abs( etot - 0.5_wp*( u1**2 + u2**2 ) )/cv

      Tmp(i,j) = T
      prs(i,j) = rho*Rair*T
      as(i,j)  = sqrt( gamma*Rair*T )
   enddo
   enddo

end subroutine compute_thermo