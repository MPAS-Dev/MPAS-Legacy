!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_pressure_grad
!
!> \brief MPAS ocean pressure gradient module
!> \author Doug Jacobsen
!> \date   15 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the routine for computing 
!>  tendencie from the horizontal pressure gradient.
!>
!
!-----------------------------------------------------------------------

module ocn_vel_pressure_grad

   use grid_types
   use configure

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

   public :: ocn_vel_pressure_grad_tend, &
             ocn_vel_pressure_grad_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real (kind=RKIND) :: rho0Inv


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_pressure_grad_tend
!
!> \brief   Computes tendency term for horizontal pressure gradient
!> \author  Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the pressure gradient tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_pressure_grad_tend(grid, pressure, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         pressure !< Input: Pressure field or Mongomery potential

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

      integer :: nEdgesSolve, iEdge, k, cell1, cell2
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(:), pointer :: dcEdge

      err = 0

      nEdgesSolve = grid % nEdgesSolve
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      cellsOnEdge => grid % cellsOnEdge % array
      dcEdge => grid % dcEdge % array

      if (config_vert_grid_type.eq.'isopycnal') then
        do iEdge=1,nEdgesSolve
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)
          do k=1,maxLevelEdgeTop(iEdge)
             tend(k,iEdge) = tend(k,iEdge)     &
               - (pressure(k,cell2) - pressure(k,cell1))/dcEdge(iEdge)
           end do
        enddo
      elseif (config_vert_grid_type.eq.'zlevel') then
        do iEdge=1,nEdgesSolve
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)
          do k=1,maxLevelEdgeTop(iEdge)

            tend(k,iEdge) = tend(k,iEdge)     &
              - rho0Inv*(  pressure(k,cell2) &
                         - pressure(k,cell1) )/dcEdge(iEdge)
          end do

        enddo
      endif


   !--------------------------------------------------------------------

   end subroutine ocn_vel_pressure_grad_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_pressure_grad_init
!
!> \brief   Initializes ocean momentum horizontal pressure gradient
!> \author  Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes parameters required for the computation of the
!>  horizontal pressure gradient.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_pressure_grad_init(err)!{{{

   !--------------------------------------------------------------------


      !-----------------------------------------------------------------
      !
      ! Output Variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag


      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      err = 0

      if (config_vert_grid_type.eq.'isopycnal') then
        rho0Inv = 1.0
      elseif (config_vert_grid_type.eq.'zlevel') then
        rho0Inv = 1.0/config_rho0
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_vel_pressure_grad_init!}}}

!***********************************************************************

end module ocn_vel_pressure_grad

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
