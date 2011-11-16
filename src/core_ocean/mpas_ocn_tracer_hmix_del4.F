!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_hmix_del4
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

module ocn_tracer_hmix_del4

   use mpas_grid_types
   use mpas_configure
   use mpas_timer

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

   public :: ocn_tracer_hmix_del4_tend, &
             ocn_tracer_hmix_del4_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: Del4On

   real (kind=RKIND) :: eddyDiff4


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_hmix_del4_tend
!
!> \brief   Computes biharmonic tendency term for horizontal tracer mixing
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the horizontal mixing tendency for tracers
!>  based on current state using a biharmonic parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hmix_del4_tend(grid, h_edge, tracers, tend, err)!{{{

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

      integer :: iEdge, nEdges, num_tracers, nVertLevels, nCells
      integer :: iTracer, k, iCell, cell1, cell2

      integer, dimension(:,:), allocatable :: boundaryMask

      integer, dimension(:), pointer :: maxLevelEdgeTop, maxLevelCell
      integer, dimension(:,:), pointer :: boundaryEdge, cellsOnEdge

      real (kind=RKIND) :: invAreaCell1, invAreaCell2, r, tracer_turb_flux, flux

      real (kind=RKIND), dimension(:,:,:), allocatable :: delsq_tracer

      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaCell, meshScalingDel4


      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if (.not.Del4On) return

      call mpas_timer_start("compute_scalar_tend-horiz diff 4")

      nEdges = grid % nEdges
      nCells = grid % nCells
      num_tracers = size(tracers, dim=1)
      nVertLevels = grid % nVertLevels

      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      maxLevelCell => grid % maxLevelCell % array
      boundaryEdge => grid % boundaryEdge % array
      cellsOnEdge => grid % cellsOnEdge % array

      dcEdge => grid % dcEdge % array
      dvEdge => grid % dvEdge % array
      areaCell => grid % areaCell % array
      meshScalingDel4 => grid % meshScalingDel4 % array

      allocate(boundaryMask(nVertLevels, nEdges+1))
      boundaryMask = 1.0
      where(boundaryEdge.eq.1) boundaryMask=0.0

      allocate(delsq_tracer(num_tracers,nVertLevels, nCells+1))

      delsq_tracer(:,:,:) = 0.

      ! first del2: div(h \nabla \phi) at cell center
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k=1,maxLevelEdgeTop(iEdge)
           do iTracer=1,num_tracers
              delsq_tracer(iTracer,k,cell1) = delsq_tracer(iTracer,k,cell1) &
                 + dvEdge(iEdge)*h_edge(k,iEdge) &
                   *(tracers(iTracer,k,cell2) - tracers(iTracer,k,cell1)) &
                   /dcEdge(iEdge) * boundaryMask(k,iEdge)
              delsq_tracer(iTracer,k,cell2) = delsq_tracer(iTracer,k,cell2) &
                 - dvEdge(iEdge)*h_edge(k,iEdge) &
                 *(tracers(iTracer,k,cell2) - tracers(iTracer,k,cell1)) &
                 /dcEdge(iEdge) * boundaryMask(k,iEdge)
           end do
         end do
      end do

      do iCell = 1,nCells
         r = 1.0 / areaCell(iCell)
         do k=1,maxLevelCell(iCell)
            do iTracer=1,num_tracers
               delsq_tracer(iTracer,k,iCell) = delsq_tracer(iTracer,k,iCell) * r
            end do
         end do
      end do

      ! second del2: div(h \nabla [delsq_tracer]) at cell center
      do iEdge=1,grid % nEdges
         cell1 = grid % cellsOnEdge % array(1,iEdge)
         cell2 = grid % cellsOnEdge % array(2,iEdge)
         invAreaCell1 = 1.0 / areaCell(cell1)
         invAreaCell2 = 1.0 / areaCell(cell2)

         do k=1,maxLevelEdgeTop(iEdge)
            do iTracer=1,num_tracers
               tracer_turb_flux = meshScalingDel4(iEdge) * eddyDiff4 &
                  *(  delsq_tracer(iTracer,k,cell2)  &
                    - delsq_tracer(iTracer,k,cell1))/dcEdge(iEdge)
               flux = dvEdge (iEdge) * tracer_turb_flux

               tend(iTracer,k,cell1) = tend(iTracer,k,cell1) & 
                  - flux * invAreaCell1 * boundaryMask(k,iEdge)
               tend(iTracer,k,cell2) = tend(iTracer,k,cell2) &
                  + flux * invAreaCell2 * boundaryMask(k,iEdge)

            enddo
         enddo
      end do

      deallocate(delsq_tracer)
      deallocate(boundaryMask)
      call mpas_timer_stop("compute_scalar_tend-horiz diff 4")
   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_del4_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_hmix_del4_init
!
!> \brief   Initializes ocean tracer horizontal mixing quantities
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  biharmonic horizontal velocity mixing in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hmix_del4_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0
      Del4on = .false.

      if ( config_h_tracer_eddy_diff4 > 0.0 ) then
          Del4On = .true.
          eddyDiff4 = config_h_tracer_eddy_diff4
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_del4_init!}}}

!***********************************************************************

end module ocn_tracer_hmix_del4

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
