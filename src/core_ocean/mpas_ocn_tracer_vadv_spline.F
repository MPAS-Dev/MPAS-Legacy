!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_vadv_spline
!
!> \brief MPAS ocean vertical tracer advection driver
!> \author Doug Jacobsen
!> \date   16 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the main driver routine for computing 
!>  vertical advection tendencies.  
!
!-----------------------------------------------------------------------

module ocn_tracer_vadv_spline

   use mpas_grid_types
   use mpas_configure
   use mpas_timer

   use ocn_tracer_vadv_spline2
   use ocn_tracer_vadv_spline3

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

   public :: ocn_tracer_vadv_spline_tend, &
             ocn_tracer_vadv_spline_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   type (timer_node), pointer :: spline2_timer, spline3_timer
   logical :: splineOn


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_vadv_spline_tend
!
!> \brief   Computes tendency term for vertical tracer advection
!> \author  Doug Jacobsen
!> \date    16 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the vertical advection tendency for tracer
!>  based on current state and user choices of advection parameterization.
!>  Multiple parameterizations may be chosen and added together.  These
!>  tendencies are generally computed by calling the specific routine
!>  for the chosen parameterization, so this routine is primarily a
!>  driver for managing these choices.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vadv_spline_tend(grid, h, wTop, tracers, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         h, &    !< Input: layer thickness
         wTop    !< Input: vertical tracer in top layer

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

      integer :: err1, err2

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.splineOn) return

      call mpas_timer_start("spline 2", .false., spline2_timer)
      call ocn_tracer_vadv_spline2_tend(grid, h, wTop, tracers, tend, err1)
      call mpas_timer_stop("spline 2", spline2_timer)

      call mpas_timer_start("spline 3", .false., spline3_timer)
      call ocn_tracer_vadv_spline3_tend(grid, h, wTop, tracers, tend, err2)
      call mpas_timer_stop("spline 3", spline3_timer)

      err = ior(err1, err2)

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vadv_spline_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_vadv_spline_init
!
!> \brief   Initializes ocean tracer vertical advection quantities
!> \author  Doug Jacobsen
!> \date    16 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  vertical tracer advection in the ocean. Since a variety of 
!>  parameterizations are available, this routine primarily calls the
!>  individual init routines for each parameterization. 
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vadv_spline_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      integer :: err1, err2

      err = 0

      splineOn = .false.

      if(config_vert_tracer_adv.eq.'spline') then
         splineOn = .true.

         call ocn_tracer_vadv_spline2_init(err1)
         call ocn_tracer_vadv_spline3_init(err2)

         err = ior(err1, err2)
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vadv_spline_init!}}}

!***********************************************************************

end module ocn_tracer_vadv_spline

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
