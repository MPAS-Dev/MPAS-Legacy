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

   use mpas_grid_types
   use mpas_configure
   use mpas_constants

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

   real (kind=RKIND) :: rho0Inv, grho0Inv


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

   subroutine ocn_vel_pressure_grad_tend(grid, pressure, zMid, rho, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         pressure, & !< Input: Pressure field or Mongomery potential
         zMid, &     !< Input: z-coordinate at mid-depth of layer
         rho         !< Input: density

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
      integer, dimension(:,:), pointer :: cellsOnEdge, edgeMask

      real (kind=RKIND), dimension(:), pointer :: dcEdge
      real (kind=RKIND) :: invdcEdge

      err = 0

      nEdgesSolve = grid % nEdgesSolve
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      cellsOnEdge => grid % cellsOnEdge % array
      dcEdge => grid % dcEdge % array
      edgeMask => grid % edgeMask % array

      ! pressure for generalized coordinates
      ! -1/rho_0 (grad p_k + rho g grad z_k^{mid})

      ! For pure isopycnal coordinates, we are still using 
      ! grad(M), the gradient of Montgomery Potential, because
      ! we have set rho0Inv=1 and grho0Inv=0 in the init routine,
      ! and pressure is passed in as MontPot.

      do iEdge=1,nEdgesSolve
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         invdcEdge = 1.0 / dcEdge(iEdge)

         do k=1,maxLevelEdgeTop(iEdge)
            tend(k,iEdge) = tend(k,iEdge)     &
              - edgeMask(k,iEdge) * rho0Inv*(  pressure(k,cell2) &
                         - pressure(k,cell1) )* invdcEdge &
              - edgeMask(k,iEdge) * grho0Inv*  0.5*(rho(k,cell1)+rho(k,cell2)) &
                        *(  zMid(k,cell2) &
                          - zMid(k,cell1) )* invdcEdge
                      
         end do

      end do



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
        grho0Inv = 0.0
      else 
        rho0Inv = 1.0/config_rho0
        grho0Inv = gravity/config_rho0
      end if



   !--------------------------------------------------------------------

   end subroutine ocn_vel_pressure_grad_init!}}}

!***********************************************************************

end module ocn_vel_pressure_grad

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
