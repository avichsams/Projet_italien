subroutine deriv_coeff_explicit
! ------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  List of coefficients for explicit derivatives
!  N.B. for centered schemes only half of the coefficients is listed because
!  they are symmetric. e.g. for 4th order a(-2) = - a(2), a(-1) = -a(1), a(0) = 0
!
!  Ref. for DRP centered schemes:
!  Bogey, Bailly: A family of low dispersive and low dissipative explicit schemes
!  for flow and noise computations, JCP 2004
! ------------------------------------------------------------------------------
   use mod_deriv
   implicit none
   ! ---------------------------------------------------------------------------
   ! Local variables
   ! 3 and 5 points
   real(wp) :: a02c_std(1), a02d_std(3), a04c_std(2), a13d_std(5), a04d_std(5)
   ! 7 points
   real(wp) :: a06c_std(3)
   real(wp) :: a06c_drp(3)
   ! 9 points
   real(wp) :: a08c_std(4)
   real(wp) :: a08c_drp(4)
   ! 11 points
   real(wp) :: a10c_std(5)
   real(wp) :: a10c_drp(5)
   ! 13 points
   real(wp) :: a12c_std(6)
   real(wp) :: a12c_drp(6)

   ! ---------------------------------------------------------------------------

   ! ---------------------------------------------------------------------------
   !                                  3 POINTS
   ! ---------------------------------------------------------------------------
   ! 2nd order, centered, standard
   a02c_std(1) =  0.5_wp
   ! 2nd order, decentered 1-3, standard
   a02d_std(1) = -3.0_wp
   a02d_std(2) =  4.0_wp
   a02d_std(3) = -1.0_wp
   a02d_std = a02d_std/2.0_wp

   ! ---------------------------------------------------------------------------
   !                                 5 POINTS
   ! ---------------------------------------------------------------------------
   ! 4th order, centered, standard
   a04c_std(1) =  8.0_wp
   a04c_std(2) = -1.0_wp
   a04c_std = a04c_std/12.0_wp
   ! 4th order, decentered 1-3, standard
   a13d_std(1) = - 3.0_wp
   a13d_std(2) = -10.0_wp
   a13d_std(3) =  18.0_wp
   a13d_std(4) = - 6.0_wp
   a13d_std(5) =   1.0_wp
   a13d_std = a13d_std/12.0_wp
   ! 4th order, decentered 0-4, standard
   a04d_std(1) = -25.0_wp
   a04d_std(2) =  48.0_wp
   a04d_std(3) = -36.0_wp
   a04d_std(4) =  16.0_wp
   a04d_std(5) = - 3.0_wp
   a04d_std = a04d_std/12.0_wp

   ! ---------------------------------------------------------------------------
   !                                 7 POINTS
   ! ---------------------------------------------------------------------------
   ! 6th order, centered, standard
   a06c_std(1) =  45.0_wp
   a06c_std(2) = - 9.0_wp
   a06c_std(3) =   1.0_wp
   a06c_std = a06c_std/60.0_wp
   ! ---------------------------------------------------------------------------
   ! optimised, centered, DRP
   a06c_drp(1) =  0.790803914666667_wp
   a06c_drp(2) = -0.182643131733333_wp
   a06c_drp(3) =  0.0248274496_wp

   ! ---------------------------------------------------------------------------
   !                                 9 POINTS
   ! ---------------------------------------------------------------------------
   ! 8th order, centered, standard
   a08c_std(1) =  672.0_wp
   a08c_std(2) = -168.0_wp
   a08c_std(3) =   32.0_wp
   a08c_std(4) = -  3.0_wp
   a08c_std = a08c_std/840.0_wp
   ! ---------------------------------------------------------------------------
   ! optimised, centered, DRP
   a08c_drp(1) =  0.841570125482_wp
   a08c_drp(2) = -0.244678631765_wp
   a08c_drp(3) =  0.059463584768_wp
   a08c_drp(4) = -0.007650904064_wp

   ! ---------------------------------------------------------------------------
   !                                 11 POINTS
   ! ---------------------------------------------------------------------------
   ! 10th order, centered, standard
   a10c_std(1) =  2100.0_wp
   a10c_std(2) = - 600.0_wp
   a10c_std(3) =   150.0_wp
   a10c_std(4) = -  25.0_wp
   a10c_std(5) =     2.0_wp
   a10c_std = a10c_std/2520.0_wp
   ! ---------------------------------------------------------------------------
   ! optimised, centered, DRP
   a10c_drp(1) =  0.872756993962_wp
   a10c_drp(2) = -0.286511173973_wp
   a10c_drp(3) =  0.090320001280_wp
   a10c_drp(4) = -0.020779405824_wp
   a10c_drp(5) =  0.002484594688_wp

   ! ---------------------------------------------------------------------------
   !                                13 POINTS
   ! ---------------------------------------------------------------------------
   ! 12th order, centered, standard
   a12c_std(1) =  23760.0_wp
   a12c_std(2) = - 7425.0_wp
   a12c_std(3) =   2200.0_wp
   a12c_std(4) = -  495.0_wp
   a12c_std(5) =     72.0_wp
   a12c_std(6) = -    5.0_wp
   a12c_std = a12c_std/27720._wp
   ! ---------------------------------------------------------------------------
   ! optimised, centered, DRP
   a12c_drp(1) =  0.907646591371_wp
   a12c_drp(2) = -0.337048393268_wp
   a12c_drp(3) =  0.133442885327_wp
   a12c_drp(4) = -0.045246480208_wp
   a12c_drp(5) =  0.011169294114_wp
   a12c_drp(6) = -0.001456501759_wp

   ! ---------------------------------------------------------------------------
   ! Assign the correct coefficients
   ! ---------------------------------------------------------------------------
   ! Coefficients for convective (high-order) derivatives
   ! ---------------------------------------------------------------------------

   ! Coefficients for low-order derivatives are only standards
   a02 = a02c_std
   a04 = a04c_std
   ! 4 order
   a13d = a13d_std
   a04d = a04d_std
   !
   if ( deriv_conv_type==EXPLICIT_STD ) then
      ! Centered
      a06 = a06c_std
      a08 = a08c_std
      a10 = a10c_std
      a12 = a12c_std
   elseif ( deriv_conv_type==EXPLICIT_DRP ) then
      a06 = a06c_drp
      a08 = a08c_drp
      a10 = a10c_drp
      a12 = a12c_drp
   endif

   ! ! ---------------------------------------------------------------------------
   ! ! Boundary treatment
   ! ! ---------------------------------------------------------------------------
   ! do idfc = 1, nfc

   !    face(idfc)%a13d = 0.0_wp
   !    face(idfc)%a04d = 0.0_wp

   !    ! Lower-order centered bnd scheme
   !    face(idfc)%a13d(1) = -0.5_wp
   !    face(idfc)%a13d(3) =  0.5_wp
   !    face(idfc)%a04d(1) = -1.0_wp
   !    face(idfc)%a04d(2) =  1.0_wp
   ! enddo

end subroutine deriv_coeff_explicit