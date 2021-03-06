!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_forcing
!
!> \brief MPAS ocean forcing driver
!> \author Doug Jacobsen
!> \date   16 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the main driver routine for computing 
!>  tendencies from forcings.  
!
!-----------------------------------------------------------------------

module ocn_vel_forcing

   use mpas_grid_types
   use mpas_configure

   use ocn_vel_forcing_windstress
   use ocn_vel_forcing_rayleigh

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

   public :: ocn_vel_forcing_tend, &
             ocn_vel_forcing_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_forcing_tend
!
!> \brief   Computes tendency term from forcings
!> \author  Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the forcing tendency for momentum
!>  based on current state and user choices of forcings.
!>  Multiple forcings may be chosen and added together.  These
!>  tendencies are generally computed by calling the specific routine
!>  for the chosen forcing, so this routine is primarily a
!>  driver for managing these choices.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_tend(grid, u, normalVelocityForcing, kineticEnergyEdge, layerThicknessEdge, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         u    !< Input: velocity

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocityForcing     !< Input: wind stress

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         kineticEnergyEdge     !< Input: kinetic energy at edge

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdge     !< Input: thickness at edge

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

      integer :: err1, err2

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      call ocn_vel_forcing_windstress_tend(grid, normalVelocityForcing, layerThicknessEdge, tend, err1)
      call ocn_vel_forcing_rayleigh_tend(grid, u, tend, err2)

      err = ior(err1, err2)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_forcing_init
!
!> \brief   Initializes ocean forcings
!> \author  Doug Jacobsen
!> \date    16 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes quantities related to forcings 
!>  in the ocean. Since a multiple forcings are available, 
!>  this routine primarily calls the
!>  individual init routines for each forcing. 
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      integer :: err1, err2

      call ocn_vel_forcing_windstress_init(err1)
      call ocn_vel_forcing_rayleigh_init(err2)

      err = ior(err1, err2)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_init!}}}

!***********************************************************************

end module ocn_vel_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
