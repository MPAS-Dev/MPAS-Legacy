!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_vadv
!
!> \brief MPAS ocean vertical advection 
!> \author Doug Jacobsen
!> \date   15 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the routine for computing 
!>  tendencies for vertical advection.
!>
!
!-----------------------------------------------------------------------

module ocn_vel_vadv

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

   public :: ocn_vel_vadv_tend, &
             ocn_vel_vadv_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: velVadvOn


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_vadv_tend
!
!> \brief   Computes tendency term for vertical advection
!> \author  Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the vertical advection tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vadv_tend(grid, u, wTop, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         u     !< Input: Horizontal velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         wTop  !< Input: Vertical velocity on top layer

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

      integer :: iEdge, nEdgesSolve, cell1, cell2, k
      integer :: nVertLevels
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND) :: wTopEdge
      real (kind=RKIND), dimension(:), allocatable :: w_dudzTopEdge
      real (kind=RKIND), dimension(:), pointer :: zMidZLevel

      if(.not.velVadvOn) return

      err = 0

      nVertLevels = grid % nVertLevels
      nEdgesSolve = grid % nEdgesSolve
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      cellsOnEdge => grid % cellsOnEdge % array
      zMidZLevel => grid % zMidZLevel % array

      allocate(w_dudzTopEdge(nVertLevels+1))
      w_dudzTopEdge = 0.0
      do iEdge=1,nEdgesSolve
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)

        do k=2,maxLevelEdgeTop(iEdge)
          ! Average w from cell center to edge
          wTopEdge = 0.5*(wTop(k,cell1)+wTop(k,cell2))

          ! compute dudz at vertical interface with first order derivative.
          w_dudzTopEdge(k) = wTopEdge * (u(k-1,iEdge)-u(k,iEdge)) &
                       / (zMidZLevel(k-1) - zMidZLevel(k))
        end do
        w_dudzTopEdge(maxLevelEdgeTop(iEdge)+1) = 0.0
        ! Average w*du/dz from vertical interface to vertical middle of cell
        do k=1,maxLevelEdgeTop(iEdge)

          tend(k,iEdge) = tend(k,iEdge) &
             - 0.5 * (w_dudzTopEdge(k) + w_dudzTopEdge(k+1))
        enddo
      enddo
      deallocate(w_dudzTopEdge)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vadv_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_vadv_init
!
!> \brief   Initializes ocean momentum vertical advection
!> \author  Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  vertical velocity advection in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vadv_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! Output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0
      velVadvOn = .false.

      if (config_vert_grid_type.eq.'zlevel') then
          velVadvOn = .true.
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vadv_init!}}}

!***********************************************************************

end module ocn_vel_vadv

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
