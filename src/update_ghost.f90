subroutine update_ghost( var, idir )
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli
!
!  DESCRIPTION
!  -----------
!  subroutine to update ghost cells for taking derivatives
!-------------------------------------------------------------------------------
   use mod_mode
   implicit none
   !----------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2), intent(inout) :: var
   integer, intent(in) :: idir
   !----------------------------------------------------------------------------

   if (idir==0) then
      if ( .not.is_boundary(1) ) then
         var(  -ngh+1:0 ,1:ny) = var(nx-ngh+1:nx,1:ny) ! send in positive x-direction
         var(nx+1:nx+ngh,1:ny) = var(   1:ngh   ,1:ny) ! send in negative x-direction
      endif
      if ( .not.is_boundary(2) ) then
         var(1:nx,  -ngh+1:0 ) = var(1:nx,ny-ngh+1:ny) ! send in positive y-direction
         var(1:nx,ny+1:ny+ngh) = var(1:nx,   1:ngh   ) ! send in negative y-direction
      endif
   ! ---------------------------------------------------------------------------
   elseif (idir==1) then
      if ( .not.is_boundary(1) ) then
         var(  -ngh+1:0 ,1:ny) = var(nx-ngh+1:nx,1:ny) ! send in positive x-direction
         var(nx+1:nx+ngh,1:ny) = var(   1:ngh   ,1:ny) ! send in negative x-direction
      endif
   ! ---------------------------------------------------------------------------
   elseif (idir==2) then
      if ( .not.is_boundary(2) ) then
         var(1:nx,  -ngh+1:0 ) = var(1:nx,ny-ngh+1:ny) ! send in positive y-direction
         var(1:nx,ny+1:ny+ngh) = var(1:nx,   1:ngh   ) ! send in negative y-direction
      endif
   else
      write(*,*) 'invalid idir in update_ghost.f90'
      stop
   endif

end subroutine update_ghost