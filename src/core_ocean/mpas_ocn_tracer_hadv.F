!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_hadv
!
!> \brief MPAS ocean horizontal tracer advection driver
!> \author Doug Jacobsen
!> \date   16 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the main driver routine for computing 
!>  horizontal advection tendencies.  
!
!-----------------------------------------------------------------------

module ocn_tracer_hadv

   use mpas_grid_types
   use mpas_configure

   use ocn_tracer_hadv2
   use ocn_tracer_hadv3
   use ocn_tracer_hadv4

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

   public :: ocn_tracer_hadv_tend, &
             ocn_tracer_hadv_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_hadv_tend
!
!> \brief   Computes tendency term for horizontal tracer advection
!> \author  Phil Jones, Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the horizontal advection tendency for tracer
!>  based on current state and user choices of advection parameterization.
!>  Multiple parameterizations may be chosen and added together.  These
!>  tendencies are generally computed by calling the specific routine
!>  for the chosen parameterization, so this routine is primarily a
!>  driver for managing these choices.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hadv_tend(grid, u, h_edge, tracers, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         u    !< Input: velocity

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

      integer :: err1, err2, err3

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      call ocn_tracer_hadv2_tend(grid, u, h_edge, tracers, tend, err1)
      call ocn_tracer_hadv3_tend(grid, u, h_edge, tracers, tend, err2)
      call ocn_tracer_hadv4_tend(grid, u, h_edge, tracers, tend, err3)

      err = ior(err1, ior(err2, err3))

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hadv_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_hadv_init
!
!> \brief   Initializes ocean tracer horizontal advection quantities
!> \author  Phil Jones, Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  horizontal velocity advection in the ocean. Since a variety of 
!>  parameterizations are available, this routine primarily calls the
!>  individual init routines for each parameterization. 
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_hadv_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      integer :: err1, err2, err3

      call ocn_tracer_hadv2_init(err1)
      call ocn_tracer_hadv3_init(err2)
      call ocn_tracer_hadv4_init(err3)

      err = ior(err1, ior(err2, err3))

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hadv_init!}}}

!***********************************************************************

end module ocn_tracer_hadv

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
