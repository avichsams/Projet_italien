subroutine derivative_x(var, dvar, ider_order)
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  compute derivative of var in the x direction
!-------------------------------------------------------------------------------
   use mod_mode
   implicit none
   ! ---------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2), intent(in)    :: var
   real(wp), dimension(  1:nx ,  1:ny ), intent(inout) :: dvar
   integer, intent(in) :: ider_order
   ! ---------------------------------------------------------------------------

   select case ( ider_order )
      case (  2 )
         ! call deriv02_x_exp( var, dvar )
         !write(*,*) 'derivative.f90: derivative order 2 to be implemented!'
         !stop
         call deriv02_x_exp( var, dvar )
      case (  4 )
         call deriv04_x_exp( var, dvar )
      case (  6 )
         call deriv06_x_exp( var, dvar )
         !write(*,*) 'derivative.f90: derivative order 6 to be implemented!'
         !stop
      case (  8 )
         ! call deriv08_x_exp( var, dvar )
         write(*,*) 'derivative.f90: derivative order 8 to be implemented!'
         stop
      case ( 10 )
         call deriv10_x_exp( var, dvar )
      case ( 12 )
         ! call deriv12_x_exp( var, dvar )
         write(*,*) 'derivative.f90: derivative order 12 to be implemented!'
         stop
   end select
end subroutine derivative_x

!-------------------------------------------------------------------------------

subroutine derivative_y(var, dvar, ider_order)
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  compute derivative of var in the y direction
!-------------------------------------------------------------------------------
   use mod_mode
   implicit none
   ! ---------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2), intent(in)    :: var
   real(wp), dimension(  1:nx ,  1:ny ), intent(inout) :: dvar
   integer, intent(in) :: ider_order
   ! ---------------------------------------------------------------------------

   select case ( ider_order )
      case (  2 )
         ! call deriv02_y_exp( var, dvar )
         !write(*,*) 'derivative.f90: derivative order 2 to be implemented!'
         !stop
         call deriv02_y_exp( var, dvar )
      case (  4 )
         call deriv04_y_exp( var, dvar )
      case (  6 )
         call deriv06_y_exp( var, dvar )
         !write(*,*) 'derivative.f90: derivative order 6 to be implemented!'
         !stop
      case (  8 )
         ! call deriv08_y_exp( var, dvar )
         write(*,*) 'derivative.f90: derivative order 8 to be implemented!'
         stop
      case ( 10 )
         call deriv10_y_exp( var, dvar )
      case ( 12 )
         ! call deriv12_y_exp( var, dvar )
         write(*,*) 'derivative.f90: derivative order 12 to be implemented!'
         stop
   end select
end subroutine derivative_y