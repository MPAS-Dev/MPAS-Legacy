!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_advection_std_vadv4
!
!> \brief MPAS ocean 4th order vertical tracer advection driver (non-monotonic/fct)
!> \author Doug Jacobsen
!> \date   03/09/12
!> \version SVN:$Id:$
!> \details
!>  This module contains routines for 4th order vertical tracer advection.
!
!-----------------------------------------------------------------------
module mpas_ocn_tracer_advection_std_vadv4

   use mpas_kind_types
   use mpas_grid_types
   use mpas_configure
   use mpas_dmpar

   use mpas_ocn_tracer_advection_helpers

   implicit none
   private
   save

   public :: mpas_ocn_tracer_advection_std_vadv4_tend

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_ocn_tracer_advection_std_vadv4_tend
!
!> \brief MPAS ocean 4th order standard vertical tracer advection tendency
!> \author Doug Jacobsen
!> \date   03/09/12
!> \version SVN:$Id:$
!> \details
!>  This routine computes the 4th order vertical tracer advection tendencies.
!
!-----------------------------------------------------------------------
   subroutine mpas_ocn_tracer_advection_std_vadv4_tend(tracers, w, verticalCellSize, grid, tend)!{{{
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tend !< Input/Output: Tracer tendency
      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers !< Input: Tracer Values
      real (kind=RKIND), dimension(:,:), intent(in) :: w !< Input: Vertical Velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: verticalCellSize !< Input: Distance between vertical interfaces of cell
      type (mesh_type), intent(in) :: grid !< Input: Grid information

      integer :: i, iCell, iEdge, k, iTracer, cell1, cell2

      real (kind=RKIND), dimension(:,:), allocatable :: vert_flux
      real (kind=RKIND) :: weightK, weightKm1
      integer :: nVertLevels, num_tracers
      integer, dimension(:), pointer :: maxLevelCell

      nVertLevels = grid % nVertLevels
      num_tracers = size(tracers, dim=1)
      maxLevelCell => grid % maxLevelCell % array

      allocate(vert_flux(num_tracers, nVertLevels+1))

      !  vertical flux divergence
      !

      ! zero fluxes at top and bottom

      vert_flux(:,1) = 0.

      do iCell=1,grid % nCellsSolve

        k = 2
        do iTracer=1,num_tracers
          weightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k-1, iCell) + verticalCellSize(k, iCell))
          weightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k-1, iCell) + verticalCellSize(k, iCell))
          vert_flux(iTracer,k) = w(k,iCell)*(weightK*tracers(iTracer,k,iCell)+weightKm1*tracers(iTracer,k-1,iCell))
        enddo
        do k=3,nVertLevels-1
          do iTracer=1,num_tracers
            vert_flux(iTracer,k) = mpas_ocn_tracer_advection_vflux4( tracers(iTracer,k-2,iCell),tracers(iTracer,k-1,iCell),  &
                                     tracers(iTracer,k  ,iCell),tracers(iTracer,k+1,iCell), w(k,iCell) )
          end do
        end do

        k = maxLevelCell(iCell)
        do iTracer=1,num_tracers
          weightK = verticalCellSize(k-1, iCell) / (verticalCellSize(k-1, iCell) + verticalCellSize(k, iCell))
          weightKm1 = verticalCellSize(k, iCell) / (verticalCellSize(k-1, iCell) + verticalCellSize(k, iCell))
          vert_flux(iTracer,k) = w(k,iCell)*(weightK*tracers(iTracer,k,iCell)+weightKm1*tracers(iTracer,k-1,iCell))
        enddo

        vert_flux(:,maxLevelCell(iCell)+1) = 0.0

        do k=1,maxLevelCell(iCell)
           do iTracer=1,num_tracers
             tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + (vert_flux(iTracer, k+1) - vert_flux(iTracer, k))
           end do
        end do

      end do

      deallocate(vert_flux)
                                                                                        
   end subroutine mpas_ocn_tracer_advection_std_vadv4_tend!}}}

end module mpas_ocn_tracer_advection_std_vadv4
