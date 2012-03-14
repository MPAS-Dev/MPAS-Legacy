!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_hadv3
!
!> \brief MPAS ocean horizontal tracer advection 3rd order
!> \author Doug Jacobsen
!> \date   16 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the main driver routine for computing 
!>  horizontal advection tendencies.  
!
!-----------------------------------------------------------------------

module ocn_tracer_hadv3

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

   public :: ocn_tracer_hadv3_tend, &
             ocn_tracer_hadv3_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: hadv3On !< Flag to turn on/off 3rd order hadv
   real (kind=RKIND) :: coef_3rd_order !< Coefficient for 3rd order hadv

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_hadv3_tend
!
!> \brief   Computes tendency term for horizontal tracer advection 3rd order
!> \author  Phil Jones, Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the horizontal advection tendency for tracer
!>  based on current state using a 3rd order formulation.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hadv3_tend(grid, u, h_edge, tracers , tend, err)!{{{

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

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, nEdges, cell1, cell2, iTracer, num_tracers, i, k, &
                 boundaryMask, velMask

      integer, dimension(:), pointer :: maxLevelEdgeTop, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnCell, &
                                          cellMask, edgeMask

      real (kind=RKIND) :: flux, tracer_edge, d2fdx2_cell1, d2fdx2_cell2, &
                           invAreaCell1, invAreaCell2

      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, areaCell
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.hadv3On) return

      nEdges = grid % nEdges
      num_tracers = size(tracers, dim=1)
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      nEdgesOnCell => grid % nEdgesOnCell % array
      cellMask => grid % cellMask % array
      cellsOnEdge => grid % cellsOnEdge % array
      cellsOnCell => grid % cellsOnCell % array
      dvEdge => grid % dvEdge % array
      dcEdge => grid % dcEdge % array
      areaCell => grid % areaCell % array
      deriv_two => grid % deriv_two % array

      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         invAreaCell1 = 1.0 / areaCell(cell1)
         invAreaCell2 = 1.0 / areaCell(cell2)

         do k=1,maxLevelEdgeTop(iEdge)

            d2fdx2_cell1 = 0.0
            d2fdx2_cell2 = 0.0

            boundaryMask = abs(transfer(cellMask(k,cell1)  == 1 .and. cellMask(k,cell2) == 1,boundaryMask))

            do iTracer=1,num_tracers

               !-- if not a boundary cell
               d2fdx2_cell1 = deriv_two(1,1,iEdge) * tracers(iTracer,k,cell1) * boundaryMask
               d2fdx2_cell2 = deriv_two(1,2,iEdge) * tracers(iTracer,k,cell2) * boundaryMask

               !-- all edges of cell 1
               do i=1,nEdgesOnCell(cell1) * boundaryMask
                  d2fdx2_cell1 = d2fdx2_cell1 + &
                  deriv_two(i+1,1,iEdge) * tracers(iTracer,k,cellsOnCell(i,cell1))
               end do

               !-- all edges of cell 2
               do i=1,nEdgesOnCell(cell2) * boundaryMask
                  d2fdx2_cell2 = d2fdx2_cell2 + &
                  deriv_two(i+1,2,iEdge) * tracers(iTracer,k,cellsOnCell(i,cell2))
               end do


               velMask = 2*(abs(transfer(u(k,iEdge) <= 0, velMask))) - 1
               flux = dvEdge(iEdge) *  u(k,iEdge) * h_edge(k,iEdge) * (          &
                    0.5*(tracers(iTracer,k,cell1) + tracers(iTracer,k,cell2))      &
                    -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                    +velMask*(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12. )

               !-- update tendency
               tend(iTracer,k,cell1) = tend(iTracer,k,cell1) - flux*invAreaCell1
               tend(iTracer,k,cell2) = tend(iTracer,k,cell2) + flux*invAreaCell2
            enddo
         end do
      end do
   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hadv3_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_hadv3_init
!
!> \brief   Initializes ocean tracer horizontal advection quantities
!> \author  Phil Jones, Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  3rd order horizontal tracer advection in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hadv3_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0
      hadv3On = .false.

      if (config_horiz_tracer_adv_order == 3) then
          hadv3On = .true.

          coef_3rd_order = config_coef_3rd_order
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hadv3_init!}}}

!***********************************************************************

end module ocn_tracer_hadv3

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
