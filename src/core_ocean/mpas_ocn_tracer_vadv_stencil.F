!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_vadv_stencil
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

module ocn_tracer_vadv_stencil

   use mpas_grid_types
   use mpas_configure
   use mpas_timer

   use ocn_tracer_vadv_stencil2
   use ocn_tracer_vadv_stencil3
   use ocn_tracer_vadv_stencil4

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

   public :: ocn_tracer_vadv_stencil_tend, &
             ocn_tracer_vadv_stencil_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   type (timer_node), pointer :: stencil2_timer, stencil3_timer, stencil4_timer

   logical :: stencilOn


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_vadv_stencil_tend
!
!> \brief   Computes tendency term for vertical tracer advection
!> \author  Phil Jones, Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the vertical advection tendency for tracers
!>  based on current state and user choices of stencil based advection parameterization.
!>  Multiple parameterizations may be chosen and added together.  These
!>  tendencies are generally computed by calling the specific routine
!>  for the chosen parameterization, so this routine is primarily a
!>  driver for managing these choices.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vadv_stencil_tend(grid, h, wTop, tracers, tend, err)!{{{

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

      integer :: err1, err2, err3

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not. stencilOn) return

      call mpas_timer_start("stencil 2", .false., stencil2_timer)
      call ocn_tracer_vadv_stencil2_tend(grid, wTop, tracers, tend, err1)
      call mpas_timer_stop("stencil 2", stencil2_timer)
      call mpas_timer_start("stencil 3", .false., stencil3_timer)
      call ocn_tracer_vadv_stencil3_tend(grid, h, wTop, tracers, tend, err2)
      call mpas_timer_stop("stencil 3", stencil3_timer)
      call mpas_timer_start("stencil 4", .false., stencil4_timer)
      call ocn_tracer_vadv_stencil4_tend(grid, h, wTop, tracers, tend, err3)
      call mpas_timer_stop("stencil 4", stencil4_timer)

      err = ior(err1, ior(err2, err3))

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vadv_stencil_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_vadv_stencil_init
!
!> \brief   Initializes ocean tracer vertical advection quantities
!> \author  Phil Jones, Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  stencil based vertical tracer advection in the ocean. Since a variety of 
!>  parameterizations are available, this routine primarily calls the
!>  individual init routines for each parameterization. 
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vadv_stencil_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      integer :: err1, err2, err3

      err = 0

      stencilOn = .false.

      if (config_vert_tracer_adv.eq.'stencil') then
         stencilOn = .true.

         call ocn_tracer_vadv_stencil2_init(err1)
         call ocn_tracer_vadv_stencil3_init(err2)
         call ocn_tracer_vadv_stencil4_init(err3)

         err = ior(err1, ior(err2, err3))
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vadv_stencil_init!}}}

!***********************************************************************

end module ocn_tracer_vadv_stencil

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
