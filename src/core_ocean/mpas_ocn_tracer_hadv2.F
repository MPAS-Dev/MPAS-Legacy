!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_hadv2
!
!> \brief MPAS ocean horizontal tracer advection 2nd order
!> \author Doug Jacobsen
!> \date   16 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the main driver routine for computing 
!>  horizontal advection tendencies.  
!
!-----------------------------------------------------------------------

module ocn_tracer_hadv2

   use grid_types
   use configure
   use timer

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

   public :: ocn_tracer_hadv2_tend, &
             ocn_tracer_hadv2_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: hadv2On !< Flag to turn on/off 2nd order hadv

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_hadv2_tend
!
!> \brief   Computes tendency term for horizontal tracer advection 2nd order
!> \author  Phil Jones, Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the horizontal advection tendency for tracer
!>  based on current state using a 2nd order formulation.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hadv2_tend(grid, u, h_edge, tracers , tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         u    !< Input: tracer

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         h_edge     !< Input: thickness at edge

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

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, nEdges, cell1, cell2, iTracer, num_tracers, k

      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND) :: flux, tracer_edge

      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.hadv2On) return

      call timer_start("compute_scalar_tend-horiz adv 2")

      nEdges = grid % nEdges
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      cellsOnEdge => grid % cellsOnEdge % array
      dvEdge => grid % dvEdge % array
      areaCell => grid % areaCell % array
      num_tracers = size(tracers, 1)

      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,maxLevelEdgeTop(iEdge)
            do iTracer=1,num_tracers
               tracer_edge = 0.5 * (tracers(iTracer,k,cell1) + tracers(iTracer,k,cell2))
               flux = u(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge) * tracer_edge
               tend(iTracer,k,cell1) = tend(iTracer,k,cell1) - flux/areaCell(cell1)
               tend(iTracer,k,cell2) = tend(iTracer,k,cell2) + flux/areaCell(cell2)
            end do
         end do
      end do

      call timer_stop("compute_scalar_tend-horiz adv 2")
   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hadv2_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_hadv2_init
!
!> \brief   Initializes ocean tracer horizontal advection quantities
!> \author  Phil Jones, Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  2nd order horizontal tracer advection in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hadv2_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      err = 0
      hadv2On = .false.

      if (config_tracer_adv_order == 2) then
          hadv2On = .true.
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hadv2_init!}}}

!***********************************************************************

end module ocn_tracer_hadv2

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
