!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_thick_hadv
!
!> \brief MPAS ocean horizontal advection for thickness
!> \author Doug Jacobsen
!> \date   16 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the routine for computing 
!>  tendencies for thickness from horizontal advection
!
!-----------------------------------------------------------------------

module ocn_thick_hadv

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

   public :: ocn_thick_hadv_tend, &
             ocn_thick_hadv_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_thick_hadv_tend
!
!> \brief   Computes tendency term from horizontal advection of thickness
!> \author  Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the horizontal advection tendency for
!>  thicknes based on current state and user choices of forcings.
!
!-----------------------------------------------------------------------

   subroutine ocn_thick_hadv_tend(grid, u, h_edge, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         u    !< Input: velocity

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         h_edge     !< Input: thickness at edge

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
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

      integer :: iEdge, nEdges, cell1, cell2, nVertLevels, k
      integer :: iCell, nCells

      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND) :: flux, invAreaCell1, invAreaCell2
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      nEdges = grid % nEdges
      nCells = grid % nCells
      nVertLevels = grid % nVertLevels

      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      cellsOnEdge => grid % cellsOnEdge % array
      dvEdge => grid % dvEdge % array
      areaCell => grid % areaCell % array

      !dwj: 10/25/2011 - Need to explore isopycnal vs zlevel flags
      if (config_vert_grid_type.eq.'isopycnal') then
 
         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            invAreaCell1 = 1.0 / areaCell(cell1)
            invAreaCell2 = 1.0 / areaCell(cell2)

            do k=1,nVertLevels
               flux = u(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge)
               tend(k,cell1) = tend(k,cell1) - flux * invAreaCell1
               tend(k,cell2) = tend(k,cell2) + flux * invAreaCell2
            end do
         end do

      elseif (config_vert_grid_type.eq.'zlevel') then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            invAreaCell1 = 1.0 / areaCell(cell1)
            invAreaCell2 = 1.0 / areaCell(cell2)

            do k=1,min(1,maxLevelEdgeTop(iEdge))
               flux = u(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge)
               tend(k,cell1) = tend(k,cell1) - flux * invAreaCell1
               tend(k,cell2) = tend(k,cell2) + flux * invAreaCell2
            end do
         end do

      endif ! config_vert_grid_type


   !--------------------------------------------------------------------

   end subroutine ocn_thick_hadv_tend!}}}

!***********************************************************************
!
!  routine ocn_thick_hadv_init
!
!> \brief   Initializes ocean horizontal thickness advection
!> \author  Doug Jacobsen
!> \date    16 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes quantities related to horizontal thickness 
!>  advection in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_thick_hadv_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_thick_hadv_init!}}}

!***********************************************************************

end module ocn_thick_hadv

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
