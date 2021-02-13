module mod_deriv
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Setup and initialization for derivative variables
! ------------------------------------------------------------------------------
   use mod_mode
   implicit none
   ! ---------------------------------------------------------------------------
   character(3) :: deriv_type_string

   ! work arrays
   real(wp), dimension(:,:), allocatable :: sx, sy, dsxdx, dsydy

   ! Derivative coefficients
   ! Centered scheme
   real(wp) :: a02(1), a04(2), a06(3), a08(4), a10(5), a12(6)
   ! Non-centered scheme
   real(wp) :: a04d(5) , a13d(5)
   real(wp) :: a24d(7) , a15d(7) , a06d(7)
   real(wp) :: a35d(9) , a26d(9) , a17d(9) , a08d(9)
   real(wp) :: a46d(11), a37d(11), a28d(11), a19d(11) , a010d(11)
   real(wp) :: a57d(13), a48d(13), a39d(13), a210d(13), a111d(13), a012d(13)

   integer, parameter :: nder02=1, nder04=2, nder06=3, nder08=4, nder10=5, nder12=6
   integer :: ider1_02, ider2_02, jder1_02, jder2_02
   integer :: ider1_04, ider2_04, jder1_04, jder2_04
   integer :: ider1_06, ider2_06, jder1_06, jder2_06
   integer :: ider1_08, ider2_08, jder1_08, jder2_08
   integer :: ider1_10, ider2_10, jder1_10, jder2_10
   integer :: ider1_12, ider2_12, jder1_12, jder2_12

   ! ---------------------------------------------------------------------------

   contains

   subroutine init_deriv
   ! ---------------------------------------------------------------------------
   !  DESCRIPTION
   !  -----------
   !  Derivative variables initialization
   ! ---------------------------------------------------------------------------
      implicit none
      ! ------------------------------------------------------------------------

      allocate( sx(nx1:nx2,ny1:ny2), dsxdx(1:nx,1:ny) &
              , sy(nx1:nx2,ny1:ny2), dsydy(1:nx,1:ny) )

      sx = 0.0_wp; dsxdx = 0.0_wp
      sy = 0.0_wp; dsydy = 0.0_wp

      ! ------------------------------------------------------------------------

      call deriv_coeff_explicit
      ! ------------------------------------------------------------------------
      ! Initialize derivative indexes
      ider1_02 =  1; jder1_02 =  1
      ider2_02 = nx; jder2_02 = ny
      if (is_boundary(1)) then
         ider1_02 =  1 + nder02
         ider2_02 = nx - nder02
      endif
      if (is_boundary(2)) then
         jder1_02 =  1 + nder02
         jder2_02 = ny - nder02
      endif
      ! ------------------------------------------------------------------------
      ider1_04 =  1; jder1_04 =  1
      ider2_04 = nx; jder2_04 = ny
      if (is_boundary(1)) then
         ider1_04 =  1 + nder04
         ider2_04 = nx - nder04
      endif
      if (is_boundary(2)) then
         jder1_04 =  1 + nder04
         jder2_04 = ny - nder04
      endif
      ! ------------------------------------------------------------------------
      ider1_06 =  1; jder1_06 =  1
      ider2_06 = nx; jder2_06 = ny
      if (is_boundary(1)) then
         ider1_06 =  1 + nder06
         ider2_06 = nx - nder06
      endif
      if (is_boundary(2)) then
         jder1_06 =  1 + nder06
         jder2_06 = ny - nder06
      endif
      ! ------------------------------------------------------------------------
      ider1_08 =  1; jder1_08 =  1
      ider2_08 = nx; jder2_08 = ny
      if (is_boundary(1)) then
         ider1_08 =  1 + nder08
         ider2_08 = nx - nder08
      endif
      if (is_boundary(2)) then
         jder1_08 =  1 + nder08
         jder2_08 = ny - nder08
      endif
      ! ------------------------------------------------------------------------
      ider1_10 =  1; jder1_10 =  1
      ider2_10 = nx; jder2_10 = ny
      if (is_boundary(1)) then
         ider1_10 =  1 + nder10
         ider2_10 = nx - nder10
      endif
      if (is_boundary(2)) then
         jder1_10 =  1 + nder10
         jder2_10 = ny - nder10
      endif
      ! ------------------------------------------------------------------------
      ider1_12 =  1; jder1_12 =  1
      ider2_12 = nx; jder2_12 = ny
      if (is_boundary(1)) then
         ider1_12 =  1 + nder12
         ider2_12 = nx - nder12
      endif
      if (is_boundary(2)) then
         jder1_12 =  1 + nder12
         jder2_12 = ny - nder12
      endif

   end subroutine init_deriv

end module mod_deriv