subroutine rhs_inviscid(phivar, dphidt)
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  compute the R.H.S. of the NS equations: convective fluxes
!-------------------------------------------------------------------------------
   use mod_mode
   use mod_deriv
   use mod_work
   implicit none
   ! ---------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2,nvars), intent(in)  :: phivar
   real(wp), dimension(  1:nx ,  1:ny ,nvars), intent(out) :: dphidt
   ! Local variables
   integer :: i, j
   ! ---------------------------------------------------------------------------

   wrk = 1.0_wp

   !----------------------------------------------------------------------------
   ! m = 1 density equation
   !----------------------------------------------------------------------------
   call deriv10_skew_x_exp( phivar(:,:,1), vel(:,:,1), wrk, dsxdx )
   call deriv10_skew_y_exp( phivar(:,:,1), vel(:,:,2), wrk, dsydy )

   ! The minus sign is to avoid temporary copy in the previous deriv calls
   dphidt(:,:,1) = -( dsxdx + dsydy )

   !----------------------------------------------------------------------------
   ! m = 2 rhou1 equation
   !----------------------------------------------------------------------------
   call deriv10_x_exp(prs, dprs(:,:,1))
   call deriv10_skew_x_exp( phivar(:,:,1), vel(:,:,1), vel(:,:,1), dsxdx )
   call deriv10_skew_y_exp( phivar(:,:,1), vel(:,:,2), vel(:,:,1), dsydy )

   dphidt(:,:,2) = -( dsxdx + dsydy + dprs(:,:,1) )
   !----------------------------------------------------------------------------
   ! m = 3 rhou2  equation
   !----------------------------------------------------------------------------
   call deriv10_y_exp(prs, dprs(:,:,2))
   call deriv10_skew_x_exp( phivar(:,:,1), vel(:,:,1), vel(:,:,2), dsxdx )
   call deriv10_skew_y_exp( phivar(:,:,1), vel(:,:,2), vel(:,:,2), dsydy )

   !----------------------------------------------------------------------------
   ! m = 5 energy equation
   !----------------------------------------------------------------------------
   wrk = ( phivar(:,:,4) + prs(:,:) ) / phivar(:,:,1)
   call deriv10_skew_x_exp( phivar(:,:,1), vel(:,:,1), wrk, dsxdx )
   call deriv10_skew_y_exp( phivar(:,:,1), vel(:,:,2), wrk, dsydy )

   ! dphidt(:,:,4) = -( dsxdx + dsydy + dprs(:,:,1) + dprs(:,:,2) + dprs(:,:,3) )
   dphidt(:,:,4) = -( dsxdx + dsydy )

end subroutine rhs_inviscid