subroutine filter_coeff_explicit
! ------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  List of coefficients for explicit filters.
!  Reference for standard filters: Kennedy & Carpenter 1994
!  Reference for DRP filter: Bogey & Bailly 2004
! ------------------------------------------------------------------------------
   use mod_filter
   implicit none
   ! ---------------------------------------------------------------------------
   ! Local variables
   ! integer  :: i
   ! Centered formulas coefficients
   real(wp) :: dc2(0:1), dc4(0:2) , dc6(0:3), dc8(0:4), dc10(0:5), dc12(0:6)
   ! Boundary points coefficients
   real(wp) ::  dd2_std(1, 2),  dd4_std(2, 4) &
             ,  dd6_std(3, 6),  dd6_drp(3, 7) &
             ,  dd8_std(4, 8),  dd8_drp(4, 9) &
             , dd10_std(5,10), dd10_drp(5,11) &
             , dd12_std(6,12), dd12_drp(6,13)
   ! ---------------------------------------------------------------------------

   if ( filter_type==explicit_std .or. filter_type==explicit_jam ) then
      ! Pay attetion: for filter_order=2,4,6 the sign is changed wrt to
      ! Kennedy & Carpenter's paper
      ! ------------------------------------------------------------------------
      !                                 3 POINTS
      ! ------------------------------------------------------------------------
      ! Interior filter 2th order
      dc2(0) = + 2.0_wp/4.0_wp
      dc2(1) = - 1.0_wp/4.0_wp
      ! Boundary filter 2th order
      dd2_std = 0.0_wp
      ! Coefficients at point 1
      dd2_std(1,1) = + 1.0_wp/4.0_wp
      dd2_std(1,2) = - 1.0_wp/4.0_wp

      ! ------------------------------------------------------------------------
      !                                 5 POINTS
      ! ------------------------------------------------------------------------
      ! Interior filter 4th order
      dc4(0) = - 6.0_wp/16.0_wp
      dc4(1) = + 4.0_wp/16.0_wp
      dc4(2) = - 1.0_wp/16.0_wp
      ! Boundary filter 4th order
      dd4_std = 0.0_wp
      ! Coefficients at point 1
      dd4_std(1,1) = - 1.0_wp/16.0_wp
      dd4_std(1,2) = + 2.0_wp/16.0_wp
      dd4_std(1,3) = - 1.0_wp/16.0_wp
      ! Coefficients at point 2
      dd4_std(2,1) = + 2.0_wp/16.0_wp
      dd4_std(2,2) = - 5.0_wp/16.0_wp
      dd4_std(2,3) = + 4.0_wp/16.0_wp
      dd4_std(2,4) = - 1.0_wp/16.0_wp
      dc4 = -dc4
      dd4_std = -dd4_std

      ! ------------------------------------------------------------------------
      !                                 7 POINTS
      ! ------------------------------------------------------------------------
      ! Interior filter 6th order
      dc6(0) = + 20.0_wp/64.0_wp
      dc6(1) = - 15.0_wp/64.0_wp
      dc6(2) = +  6.0_wp/64.0_wp
      dc6(3) = -  1.0_wp/64.0_wp
      ! ------------------------------------------------------------------------
      !                                 9 POINTS
      ! ------------------------------------------------------------------------
      ! Interior filter 8th order
      dc8(0) = - 70.0_wp/256.0_wp
      dc8(1) = + 56.0_wp/256.0_wp
      dc8(2) = - 28.0_wp/256.0_wp
      dc8(3) = +  8.0_wp/256.0_wp
      dc8(4) = -  1.0_wp/256.0_wp
      dc8 = -dc8
      ! ------------------------------------------------------------------------
      !                                11 POINTS
      ! ------------------------------------------------------------------------
      ! Interior filter 10th order
      dc10(0) = + 252.0_wp/1024.0_wp
      dc10(1) = - 210.0_wp/1024.0_wp
      dc10(2) = + 120.0_wp/1024.0_wp
      dc10(3) = -  45.0_wp/1024.0_wp
      dc10(4) = +  10.0_wp/1024.0_wp
      dc10(5) = -   1.0_wp/1024.0_wp
      ! ------------------------------------------------------------------------
      !                                13 POINTS
      ! ------------------------------------------------------------------------
      ! Interior filter 12th order
      dc12(0) = - 924.0_wp/4096.0_wp
      dc12(1) = + 792.0_wp/4096.0_wp
      dc12(2) = - 495.0_wp/4096.0_wp
      dc12(3) = + 220.0_wp/4096.0_wp
      dc12(4) = -  66.0_wp/4096.0_wp
      dc12(5) = +  12.0_wp/4096.0_wp
      dc12(6) = -   1.0_wp/4096.0_wp
      dc12 = -dc12

   ! ---------------------------------------------------------------------------
   elseif (filter_type==explicit_drp) then
      ! ------------------------------------------------------------------------
      !                                 7 POINTS
      ! ------------------------------------------------------------------------
      ! Interior filter 7 points (o(4) optimized ) Tam
      dc6(0) =  0.287392842460_wp
      dc6(1) = -0.226146951809_wp
      dc6(2) =  0.106303578770_wp
      dc6(3) = -0.023853048191_wp
      ! Boundary optimized filters
      dd6_drp = 0.0_wp
      ! Coefficients at point 1 (decentre 0-3), 4 points
      dd6_drp(1,1) =  0.320882352941_wp
      dd6_drp(1,2) = -0.465000000000_wp
      dd6_drp(1,3) =  0.179117647059_wp
      dd6_drp(1,4) = -0.035000000000_wp
      ! Filter at point 1 is too much dissipative. Berland suggests to use
      ! 1/10 of the filter strength
      dd6_drp(1,:) =  dd6_drp(1,:)/10.0_wp
      ! Coefficients at point 2
      ! Filtre de Berland modifie (decentre 1-5), 7 points
      dd6_drp(2,1) = -0.085777408970_wp + 0.000000000001_wp
      dd6_drp(2,2) =  0.277628171524_wp
      dd6_drp(2,3) = -0.356848072173_wp
      dd6_drp(2,4) =  0.223119093072_wp
      dd6_drp(2,5) = -0.057347064865_wp
      dd6_drp(2,6) = -0.000747264596_wp
      dd6_drp(2,7) = -0.000027453993_wp
      ! Coefficients at point 3
      ! Filtre de Berland (decentre 2-4), 7 points
      ! Attention: correction on dd6_drp(3,4)
      dd6_drp(3,1) =  0.032649010764_wp
      dd6_drp(3,2) = -0.143339502575_wp
      dd6_drp(3,3) =  0.273321177980_wp
      dd6_drp(3,4) = -0.294622121167_wp - 0.000000000002_wp
      dd6_drp(3,5) =  0.186711738069_wp
      dd6_drp(3,6) = -0.062038376258_wp
      dd6_drp(3,7) =  0.007318073189_wp
      ! ------------------------------------------------------------------------
      !                                 9 POINTS
      ! ------------------------------------------------------------------------
      ! Interior filter 9 points (o(4) optimized )
      dc8(0) =  0.243527493120_wp
      dc8(1) = -0.204788880640_wp
      dc8(2) =  0.120007591680_wp
      dc8(3) = -0.045211119360_wp
      dc8(4) =  0.008228661760_wp
      ! ------------------------------------------------------------------------
      !                                11 POINTS
      ! ------------------------------------------------------------------------
      ! Interior filter 11 points (o(6) optimized ) (Bogey Bailly JCP 2009)
      dc10(0) =  0.234810479761700_wp
      dc10(1) = -0.199250131285813_wp
      dc10(2) =  0.120198310245186_wp
      dc10(3) = -0.049303775636020_wp
      dc10(4) =  0.012396449873964_wp
      dc10(5) = -0.001446093078167_wp
      ! ------------------------------------------------------------------------
      !                                13 POINTS
      ! ------------------------------------------------------------------------
      ! Interior filter 13 points (o(4) optimized )
      dc12(0) =  0.190899511506_wp
      dc12(1) = -0.171503832236_wp
      dc12(2) =  0.123632891797_wp
      dc12(3) = -0.069975429105_wp
      dc12(4) =  0.029662754736_wp
      dc12(5) = -0.008520738659_wp
      dc12(6) =  0.001254597714_wp
      ! ------------------------------------------------------------------------
   endif
   ! ---------------------------------------------------------------------------
   ! Selection of the chosen filter
   dc = 0.0_wp
   dd = 0.0_wp
   if (filter_order==2) then
      dc = dc2
      dd = dd2_std
   elseif (filter_order==4) then
      dc = dc4
      dd = dd4_std
   elseif (filter_order==6) then
      dc = dc6
      if (filter_type==explicit_std) then
         dd = dd6_std
      else
         dd = dd6_drp
      endif
   elseif (filter_order==8) then
      dc = dc8
      if (filter_type==explicit_std) then
         dd = dd8_std
      else
         dd = dd8_drp
      endif
   elseif (filter_order==10) then
      dc = dc10
      if (filter_type==explicit_std) then
         ! dd = dd10_std
         dd(5,5:9) = dc8(0:4)
         dd(5,4) = dc8(1); dd(5,3) = dc8(2); dd(5,2) = dc8(3); dd(5,1) = dc8(4)
         dd(4,4:7) = dc6(0:3)
         dd(4,3) = dc6(1); dd(4,2) = dc6(2); dd(4,1) = dc6(3)
         dd(3,3:5) = dc4(0:2)
         dd(3,2) = dc4(1); dd(3,1) = dc4(2)
         dd(2,1:4) = dd4_std(2,1:4)
         ! dd(2,2:3) = dc2(0:1)*2.0_wp
         ! dd(2,1) = dc2(1)*2.0_wp
      else
         dd = dd10_drp
      endif
   elseif (filter_order==12) then
      dc = dc12
      if (filter_type==explicit_std) then
         dd = dd12_std
      else
         dd = dd12_drp
      endif
   else
      write(*,*) 'incorrect filter_order'
      stop
   endif

   ! Divide by two the central point coefficient
   dc(0) = dc(0)/2.0_wp
   ! Multiply by the amplitude
   dc = dc*fltamp
   dd = dd*fltamp

end subroutine filter_coeff_explicit

! ----

subroutine filter_coeff_shock
! ------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  List of coefficients for explicit filters.
!  Reference for standard filters: Kennedy & Carpenter 1994
!  Reference for DRP filter: Bogey & Bailly 2004
! ------------------------------------------------------------------------------
   use mod_filter
   implicit none
   ! ---------------------------------------------------------------------------
   ! Local variables
   ! ---------------------------------------------------------------------------
   ! coef filtre conservative opt. shock (Bogey JCP 2009)
   d12( 1) =-0.210383_wp
   d12( 2) = 0.039617_wp
   d12(-1) =-d12(2)
   d12( 0) =-d12(1)

   ! ! coef filtre conservative std. shock
   ! d12( 1) =-0.25_wp
   ! d12( 2) = 0.0_wp
   ! d12(-1) =-d12(2)
   ! d12( 0) =-d12(1)

end subroutine filter_coeff_shock