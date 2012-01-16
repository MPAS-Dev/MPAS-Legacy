!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_vadv_stencil2
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

module ocn_tracer_vadv_stencil2

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

   public :: ocn_tracer_vadv_stencil2_tend, &
             ocn_tracer_vadv_stencil2_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: stencil2On


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_vadv_stencil2_tend
!
!> \brief   Computes tendency term for vertical tracer advection 2nd order stencil
!> \author  Phil Jones, Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the vertical advection tendency for tracer
!>  based on current state using a 2nd order stencil.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vadv_stencil2_tend(grid, wTop, tracers, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
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

      integer :: nCellsSolve, iCell, k, iTracer, num_tracers, nVertLevels
      integer :: nCells

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

      if(.not.stencil2On) return

      nCells = grid % nCells
      nCellsSolve = grid % nCellsSolve
      num_tracers = size(tracers, 1)
      nVertLevels = grid % nVertLevels
      maxLevelCell => grid % maxLevelCell % array

      allocate(tracerTop(num_tracers,nVertLevels+1,nCells))

      ! Compute tracerTop using centered stencil, a simple average.

      do iCell=1,nCellsSolve 
         tracerTop(:,1,iCell) = tracers(:,1,iCell)
         do k=2,maxLevelCell(iCell)
            do iTracer=1,num_tracers
               tracerTop(iTracer,k,iCell) = &
                  ( tracers(iTracer,k-1,iCell) &
                   +tracers(iTracer,k  ,iCell))/2.0
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

   end subroutine ocn_tracer_vadv_stencil2_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_vadv_stencil2_init
!
!> \brief   Initializes ocean tracer vertical advection quantities
!> \author  Phil Jones, Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  a 2nd order stencil based vertical tracer advection in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vadv_stencil2_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      integer :: err1, err2, err3

      err = 0
      stencil2On = .false.

      if(config_vert_tracer_adv_order.eq.2) then
          stencil2On = .true.
      endif


   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vadv_stencil2_init!}}}

!***********************************************************************

end module ocn_tracer_vadv_stencil2

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
