!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_thick_vadv
!
!> \brief MPAS ocean vertical advection for thickness
!> \author Doug Jacobsen
!> \date   16 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the routine for computing 
!>  tendencies for thickness from vertical advection
!
!-----------------------------------------------------------------------

module ocn_thick_vadv

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

   public :: ocn_thick_vadv_tend, &
             ocn_thick_vadv_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_thick_vadv_tend
!
!> \brief   Computes tendency term from vertical advection of thickness
!> \author  Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the vertical advection tendency for
!>  thicknes based on current state and user choices of forcings.
!
!-----------------------------------------------------------------------

   subroutine ocn_thick_vadv_tend(grid, wTop, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         wTop     !< Input: vertical velocity on top layer

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

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, nCells

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      nCells = grid % nCells

      ! dwj 10/25/2011 - Need to explore isopycnal vs zlevel flags
      if (config_vert_grid_type.eq.'zlevel') then
        do iCell=1,nCells
           tend(1,iCell) =   tend(1,iCell) + wTop(2,iCell)
        end do
      endif ! coordinate type


   !--------------------------------------------------------------------

   end subroutine ocn_thick_vadv_tend!}}}

!***********************************************************************
!
!  routine ocn_thick_vadv_init
!
!> \brief   Initializes ocean thickness vertical advection
!> \author  Doug Jacobsen
!> \date    16 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes quantities related to vertical advection of 
!>  thickness in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_thick_vadv_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag
      
      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_thick_vadv_init!}}}

!***********************************************************************

end module ocn_thick_vadv

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
