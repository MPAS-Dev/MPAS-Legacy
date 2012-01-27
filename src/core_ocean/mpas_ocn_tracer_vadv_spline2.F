!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_vadv_spline2
!
!> \brief MPAS ocean vertical tracer advection driver
!> \author Doug Jacobsen
!> \date   16 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the main driver routine for computing 
!>  vertical advection tendencies.  
!
!-----------------------------------------------------------------------

module ocn_tracer_vadv_spline2

   use mpas_grid_types
   use mpas_configure

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_tracer_vadv_spline2_tend, &
             ocn_tracer_vadv_spline2_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: spline2On


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_vadv_spline2_tend
!
!> \brief   Computes tendency term for vertical tracer advection 2nd order spline
!> \author  Doug Jacobsen
!> \date    16 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the vertical advection tendency for tracer
!>  based on current state using a 2nd order spline.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vadv_spline2_tend(grid, h, wTop, tracers, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         h, &    !< Input: layer thickness
         wTop    !< Input: vertical tracer in top layer

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers     !< Input: tracers

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend          !< Input/Output: tracer tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, nCells, nCellsSolve, k, iTracer, num_tracers, nVertLevels

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND), dimension(:,:,:), allocatable :: tracerTop

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.spline2On) return
      ! Compute tracerTop using linear interpolation.

      nCells = grid % nCells
      nCellsSolve = grid % nCellsSolve
      nVertLevels = grid % nVertLevels
      num_tracers = size(tracers, 1)
      maxLevelCell => grid % maxLevelCell % array

      allocate(tracerTop(num_tracers,nVertLevels+1,nCells))

      do iCell=1,nCellsSolve 
         tracerTop(:,1,iCell) = tracers(:,1,iCell) 
         do k=2,maxLevelCell(iCell)
            do iTracer=1,num_tracers
               ! Note h on the k side is multiplied by tracer at k-1
               ! and h on the Km1 (k-1) side is mult. by tracer at k.
               tracerTop(iTracer,k,iCell) = &
                  (  h(k  ,iCell)*tracers(iTracer,k-1,iCell) &
                   + h(k-1,iCell)*tracers(iTracer,k  ,iCell) ) &
                  / (h(k-1,iCell) + h(k,iCell))   
            end do
         end do
         tracerTop(:,maxLevelCell(iCell)+1,iCell) = tracers(:,maxLevelCell(iCell),iCell)
      end do

      do iCell=1,nCellsSolve 
         do k=1,maxLevelCell(iCell)  
            do iTracer=1,num_tracers
               tend(iTracer,k,iCell) = tend(iTracer,k,iCell) &
                  - (   wTop(k  ,iCell)*tracerTop(iTracer,k  ,iCell) &
                      - wTop(k+1,iCell)*tracerTop(iTracer,k+1,iCell))
            end do
         end do
      end do

      deallocate(tracerTop)
   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vadv_spline2_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_vadv_spline2_init
!
!> \brief   Initializes ocean tracer vertical advection quantities
!> \author  Doug Jacobsen
!> \date    16 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  2nd order spline based vertical tracer advection in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vadv_spline2_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      spline2On = .false.

      if(config_vert_tracer_adv_order.eq.2) then
        spline2On = .true.
      endif


   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vadv_spline2_init!}}}

!***********************************************************************

end module ocn_tracer_vadv_spline2

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
