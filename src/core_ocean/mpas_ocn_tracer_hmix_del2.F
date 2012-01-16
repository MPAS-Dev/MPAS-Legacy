!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_hmix_del2
!
!> \brief MPAS ocean horizontal tracer mixing driver
!> \author Doug Jacobsen
!> \date   19 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the main driver routine for computing 
!>  horizontal mixing tendencies.  
!>
!>  It provides an init and a tend function. Each are described below.
!
!-----------------------------------------------------------------------

module ocn_tracer_hmix_del2

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

   public :: ocn_tracer_hmix_del2_tend, &
             ocn_tracer_hmix_del2_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: del2On

   real (kind=RKIND) :: eddyDiff2


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_hmix_del2_tend
!
!> \brief   Computes laplacian tendency term for horizontal tracer mixing
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the horizontal mixing tendency for tracers
!>  based on current state using a laplacian parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hmix_del2_tend(grid, h_edge, tracers, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         h_edge    !< Input: thickness at edge

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
        tracers !< Input: tracer quantities

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

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

      integer :: iEdge, nEdges, nVertLevels, cell1, cell2
      integer :: k, iTracer, num_tracers

      integer, dimension(:,:), allocatable :: boundaryMask

      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge, edgeMask

      real (kind=RKIND) :: invAreaCell1, invAreaCell2
      real (kind=RKIND) :: tracer_turb_flux, flux, r_tmp

      real (kind=RKIND), dimension(:), pointer :: areaCell, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: meshScalingDel2

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if (.not.del2On) return

      nEdges = grid % nEdges
      nVertLevels = grid % nVertLevels
      num_tracers = size(tracers, dim=1)

      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      cellsOnEdge => grid % cellsOnEdge % array
      edgeMask => grid % edgeMask % array
      areaCell => grid % areaCell % array
      dvEdge => grid % dvEdge % array
      dcEdge => grid % dcEdge % array
      meshScalingDel2 => grid % meshScalingDel2 % array

      !
      ! compute a boundary mask to enforce insulating boundary conditions in the horizontal
      !
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         invAreaCell1 = 1.0/areaCell(cell1)
         invAreaCell2 = 1.0/areaCell(cell2)

         r_tmp = meshScalingDel2(iEdge) * eddyDiff2 * dvEdge(iEdge) / dcEdge(iEdge)

         do k=1,maxLevelEdgeTop(iEdge)
           do iTracer=1,num_tracers
              ! \kappa_2 \nabla \phi on edge
              tracer_turb_flux = tracers(iTracer,k,cell2) - tracers(iTracer,k,cell1)

              ! div(h \kappa_2 \nabla \phi) at cell center
              flux = h_edge(k,iEdge) * tracer_turb_flux * edgeMask(k, iEdge) * r_tmp

              tend(iTracer,k,cell1) = tend(iTracer,k,cell1) + flux * invAreaCell1
              tend(iTracer,k,cell2) = tend(iTracer,k,cell2) - flux * invAreaCell2
           end do
         end do

      end do
   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_del2_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_hmix_del2_init
!
!> \brief   Initializes ocean tracer horizontal mixing quantities
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  laplacian horizontal velocity mixing in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hmix_del2_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      del2on = .false.

      if ( config_h_tracer_eddy_diff2 > 0.0 ) then
          del2On = .true.
          eddyDiff2 = config_h_tracer_eddy_diff2
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_del2_init!}}}

!***********************************************************************

end module ocn_tracer_hmix_del2

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
