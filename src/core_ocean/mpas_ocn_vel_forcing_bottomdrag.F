!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_forcing_bottomdrag
!
!> \brief MPAS ocean bottom drag
!> \author Doug Jacobsen
!> \date   16 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the routine for computing 
!>  tendencies from bottom drag.  
!
!-----------------------------------------------------------------------

module ocn_vel_forcing_bottomdrag

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

   public :: ocn_vel_forcing_bottomdrag_tend, &
             ocn_vel_forcing_bottomdrag_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: bottomDragOn
   real (kind=RKIND) :: bottomDragCoef


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_forcing_bottomdrag_tend
!
!> \brief   Computes tendency term from bottom drag
!> \author  Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the bottom drag tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_bottomdrag_tend(grid, u, ke_edge, h_edge, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         u    !< Input: velocity 

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         ke_edge     !< Input: kinetic energy at edge

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

      integer :: iEdge, nEdgesSolve, k
      integer, dimension(:), pointer :: maxLevelEdgeTop

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.bottomDragOn) return

      nEdgesSolve = grid % nEdgesSolve
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array

      do iEdge=1,grid % nEdgesSolve

        k = maxLevelEdgeTop(iEdge)

        ! efficiency note: it would be nice to avoid this
        ! if within a do.  This could be done with
        ! k =  max(maxLevelEdgeTop(iEdge),1)
        ! and then tend_u(1,iEdge) is just not used for land cells.

        if (k>0) then
           ! bottom drag is the same as POP:
           ! -c |u| u  where c is unitless and 1.0e-3.
           ! see POP Reference guide, section 3.4.4.

           tend(k,iEdge) = tend(k,iEdge)  &
               -bottomDragCoef*u(k,iEdge) &
               *sqrt(2.0*ke_edge(k,iEdge))/h_edge(k,iEdge)

        endif

      enddo



   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_bottomdrag_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_forcing_bottomdrag_init
!
!> \brief   Initializes ocean bottom drag
!> \author  Doug Jacobsen
!> \date    16 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes quantities related to bottom drag 
!>  in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_bottomdrag_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag


      err = 0

      bottomDragOn = .false.

      if (.not.config_implicit_vertical_mix) then
          bottomDragOn = .true.
          bottomDragCoef = config_bottom_drag_coeff
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_bottomdrag_init!}}}

!***********************************************************************

end module ocn_vel_forcing_bottomdrag

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
