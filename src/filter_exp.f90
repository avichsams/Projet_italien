subroutine filter_exp( var )
!-------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Apply explicit n-th order filter to conservative variables vector
!-------------------------------------------------------------------------------
   use mod_filter
   implicit none
   ! ---------------------------------------------------------------------------
   ! Input/Output arguments
   real(wp), dimension(nx1:nx2,ny1:ny2,nvars), intent(inout) :: var
   ! Local variables
   integer, parameter :: nstart = 1
   integer :: i, j, m, n, n2
   ! ---------------------------------------------------------------------------

   ! ---------------------------------------------------------------------------
   ! filtering along x direction
   ! ===========================
   do m = 1,nvars
      call update_ghost( var(:,:,m), 1 )
      ! Inner points
      do j = 1,ny
      do i = iflt1,iflt2
         fltvar(i,j) = 0.0_wp
         do n = 0, nflt
            fltvar(i,j) = fltvar(i,j) + dc(n)*( var(i+n,j,m)+var(i-n,j,m) )
         enddo
      enddo
      enddo
      ! ------------------------------------------------------------------------
      if ( is_boundary(1) ) then
         ! left boundary
         do j = 1,ny
            ! Loop over number of boundary points (n)
            do n = nstart, nflt
               fltvar(n,j) = 0.0_wp
               ! Loop over number of points considered
               do n2 = 1,nflt2
                  fltvar(n,j) = fltvar(n,j) + dd(n,n2)*var(n2,j,m)
               enddo
            enddo
         enddo
         ! ---------------------------------------------------------------------
         ! right boundary
         do j = 1,ny
            ! Loop over number of boundary points (n)
            do n = nstart, nflt
               fltvar(nx-n+1,j) = 0.0_wp
               ! Loop over number of points considered
               do n2 = 1,nflt2
                  fltvar(nx-n+1,j) =  fltvar(nx-n+1,j) + dd(n,n2)*var(nx-n2+1,j,m)
               enddo
            enddo
         enddo
      endif
      ! update var
      var(1:nx,1:ny,m) = var(1:nx,1:ny,m) - fltvar(1:nx,1:ny)
   enddo

   if (ndim==1) return
   ! ------------------------------------------------------------------------
   ! filtering along y direction
   ! ===========================
   do m = 1,nvars
      call update_ghost( var(:,:,m), 2 )
      ! Inner points
      do j = jflt1,jflt2
      do i = 1,nx
         fltvar(i,j) = 0.0_wp
         do n = 0, nflt
            fltvar(i,j) = fltvar(i,j) + dc(n)*( var(i,j+n,m)+var(i,j-n,m) )
         enddo
      enddo
      enddo
      ! ------------------------------------------------------------------------
      if ( is_boundary(2) ) then
         ! bottom boundary
         do i = 1,nx
            ! Loop over number of boundary points (n)
            do n = nstart, nflt
               fltvar(i,n) = 0.0_wp
               ! Loop over number of points considered
               do n2 = 1,nflt2
                  fltvar(i,n) = fltvar(i,n) + dd(n,n2)*var(i,n2,m)
               enddo
            enddo
         enddo
         ! ---------------------------------------------------------------------
         ! top boundary
         do i = 1,nx
            ! Loop over number of boundary points (n)
            do n = nstart, nflt
               fltvar(i,ny-n+1) = 0.0_wp
               ! Loop over number of points considered
               do n2 = 1,nflt2
                  fltvar(i,ny-n+1) = fltvar(i,ny-n+1) + dd(n,n2)*var(i,ny-n2+1,m)
               enddo
            enddo
         enddo
      endif
      ! update var
      var(1:nx,1:ny,m) = var(1:nx,1:ny,m) - fltvar(1:nx,1:ny)
   enddo

end subroutine filter_exp