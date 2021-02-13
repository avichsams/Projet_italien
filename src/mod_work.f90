module mod_work
! ------------------------------------------------------------------------------
!  AUTHOR
!  ------
!  Luca Sciacovelli (luca.sciacovelli@ensam.eu)
!
!  DESCRIPTION
!  -----------
!  Setup and initialization of common working arrays
! ------------------------------------------------------------------------------
   use mod_mode
   implicit none
   real(wp), allocatable :: phi(:,:,:)   &   ! Conservative variables
                    ,       vel(:,:,:)   &   ! Velocities
                    ,        as(:,:)     &   ! Speed of sound
                    ,  bulk_lad(:,:)     &   ! Bulk viscosity
                    ,lambda_lad(:,:)     &   ! Artificial conductivity
                    ,       prs(:,:)     &   ! Pressure
                    ,       tmp(:,:)     &   ! Temperature
                    ,       div(:,:)     &   ! Velocity Divergence
                    ,      drho(:,:,:)   &   ! Density gradient
                    ,      dprs(:,:,:)   &   ! Pressure gradient
                    ,      dvel(:,:,:,:) &   ! Velocity gradient
                    ,    stress(:,:,:,:) &   ! Stress tensor
                    ,     qflux(:,:,:)   &   ! Heat flux
                    ,       wrk(:,:)         ! 2 dim work array
   contains

   subroutine init_work
   ! ---------------------------------------------------------------------------
   !  DESCRIPTION
   !  -----------
   !  General variables initialization
   ! ---------------------------------------------------------------------------
      implicit none

      ! Arrays with ghost-cells
      allocate(   wrk(nx1:nx2,ny1:ny2)     &
              ,   prs(nx1:nx2,ny1:ny2)     &
              ,   tmp(nx1:nx2,ny1:ny2)     &
              ,   div(nx1:nx2,ny1:ny2)     &
              ,   vel(nx1:nx2,ny1:ny2,2)   &
              ,stress(nx1:nx2,ny1:ny2,2,2) &
              ,   phi(nx1:nx2,ny1:ny2,4) )

      ! Arrays without ghost-cells
      allocate(         as(1:nx,1:ny)     &
              ,   bulk_lad(1:nx,1:ny)     &
              , lambda_lad(1:nx,1:ny)     &
              ,       drho(1:nx,1:ny,2)   &
              ,      qflux(1:nx,1:ny,2)   &
              ,       dvel(1:nx,1:ny,2,2) )

      prs    = 1.0_wp
      tmp    = 0.0_wp
      div    = 0.0_wp
      vel    = 0.0_wp
      phi    = 1.0_wp
      drho   = 0.0_wp
      dvel   = 0.0_wp
      stress = 0.0_wp
      qflux  = 0.0_wp
      bulk_lad   = 0.0_wp
      lambda_lad = 0.0_wp

   end subroutine init_work

end module mod_work