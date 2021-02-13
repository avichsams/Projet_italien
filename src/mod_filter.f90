module mod_filter
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Setup and initialization for filter variables
!-------------------------------------------------------------------------------
use mod_mode
implicit none

integer   :: filter_type         ! Type of filter
integer   :: filter_freq         ! Frequency of filter
real(wp)  :: fltamp              ! Filter Amplitude
real(wp)  :: fltshock_amp        ! Shock capturing Amplitude / threshold (J/B)
real(wp)  :: c_k, c_beta   ! Coefficients for LAD

! Filter coefficients
real(wp), allocatable :: dc(:), dd(:,:)
! Centered scheme
real(wp) :: f02(2)
! Non-centered scheme
real(wp) :: f02d(2)

integer  :: nflt, nflt2
integer  :: iflt1, iflt2, jflt1, jflt2
integer  :: igaus1, igaus2, jgaus1, jgaus2
integer  :: ibog1, ibog2, jbog1, jbog2
! Coefficient for shock capturing
real(wp) :: d12(-1:2)

! Work arrays
real(wp), dimension(:,:), allocatable :: fltvar, ducros
! Arrays for Jameson
real(wp), dimension(:,:,:), allocatable :: psens, rspec

contains

   subroutine init_filter
   ! ---------------------------------------------------------------------------
   !  DESCRIPTION
   !  -----------
   !  Filter variables initialization
   ! ---------------------------------------------------------------------------
      implicit none
      ! ------------------------------------------------------------------------

      ! nflt: nmb of points to consider at each side for interior filter
      nflt  = filter_order/2
      ! nflt2: maximum nmb of points to consider for one-sided boundary filter
      if (filter_type==explicit_std) then
         nflt2 = filter_order
      else
         nflt2 = filter_order+1
      endif

      allocate(dc(0:nflt), dd(nflt,nflt2))

      dc = 0.0_wp
      dd = 0.0_wp

      call filter_coeff_explicit

      ! Initialize filter indexes
      iflt1 =  1; jflt1 =  1
      iflt2 = nx; jflt2 = ny
      if (is_boundary(1)) then
         iflt1 =  1 + nflt
         iflt2 = nx - nflt
      endif
      if (is_boundary(2)) then
         jflt1 =  1 + nflt
         jflt2 = ny - nflt
      endif

      if (is_shock) then

         select case (fltshock_type)

            case ('B')
               call filter_coeff_shock

               ibog1 =  1; jbog1 =  1
               ibog2 = nx; jbog2 = ny
               if (is_boundary(1)) then
                  ibog1 =  1 + 2
                  ibog2 = nx - 2
               endif
               if (is_boundary(2)) then
                  jbog1 =  1 + 2
                  jbog2 = ny - 2
               endif

            case ('J')

               allocate(ducros(1:nx,1:ny)         &
                       , psens(nx1:nx2,ny1:ny2,2) &
                       , rspec(nx1:nx2,ny1:ny2,2) )
               ducros = 0.0_wp
               psens  = 0.0_wp
               rspec  = 0.0_wp

            case ('L')
               ! Index boundaries for truncated gaussian filter
               igaus1 =  1; jgaus1 =  1
               igaus2 = nx; jgaus2 = ny
               if (is_boundary(1)) then
                  igaus1 =  1 + 4
                  igaus2 = nx - 4
               endif
               if (is_boundary(2)) then
                  jgaus1 =  1 + 4
                  jgaus2 = ny - 4
               endif

            case default
               write(*,*) 'mod_filter.f90: Shock type not known'
               stop
         end select

      endif

      ! Allocate work arrays
      allocate( fltvar(1:nx,1:ny) )
      fltvar = 0.0_wp

   end subroutine init_filter

end module mod_filter