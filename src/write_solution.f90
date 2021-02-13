subroutine write_solution( phivar )
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Write solution to file.
!-------------------------------------------------------------------------------
   use mod_mode
   use mod_grid
   use mod_time
   implicit none
   ! ---------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2,4), intent(in) :: phivar
   ! Local variables
   integer, parameter :: io = 20
   integer          :: i, j, n
   logical          :: iexist
   character(len=30):: solname  = 'restart' &
                     , gridname = 'grid'
   ! ---------------------------------------------------------------------------
   ! Write grid file if it does not exist
   inquire(FILE=trim(gridname), EXIST=iexist)

   if ( .not.iexist ) then
      open( UNIT=io, FILE=trim(gridname)//'.bin', FORM='unformatted', ACCESS='stream' )
      write(io) nx
      write(io) ny
      write(io) ( gridx(i), i=1,nx )
      write(io) ( gridy(j), j=1,ny )
      close(io)
   endif

   call calcfilestamp_fromtime( tstar, filestamp )

   select case ( flowtype )

      case ( VORTEX_CONVECTION, SHOCK_VORTEX_INT, TAYLOR_GREEN_VORTEX )
         write(*,*) 'Write solution to file ~> ', trim(solname)//trim(filestamp)//'.bin'
         ! Write solution file at the actual time
         open( UNIT=io, FILE=trim(solname)//trim(filestamp)//'.bin', FORM='unformatted', ACCESS='stream' )
         write(io) nx
         write(io) ny
         write(io) ((( phivar(i,j,n), i=1,nx ), j=1,ny), n=1,nvars)
         close(io)

      case ( SHOCK_TUBE, SHU_OSHER )
         write(*,*) 'Write solution to file ~> ', trim(solname)//trim(filestamp)//'.dat'
         open( UNIT=io, FILE=trim(solname)//trim(filestamp)//'.dat', FORM='formatted' )
         do i=1,nx
            write(io,'(5(1X,G0))') gridx(i), (phivar(i,1,n), n=1,nvars)
         enddo
      case default

   end select

end subroutine write_solution