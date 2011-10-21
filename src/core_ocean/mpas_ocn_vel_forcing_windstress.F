!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_forcing_windstress
!
!> \brief MPAS ocean wind stress
!> \author Doug Jacobsen
!> \date   16 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the routine for computing 
!>  tendencies from wind stress.  
!
!-----------------------------------------------------------------------

module ocn_vel_forcing_windstress

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

   public :: ocn_vel_forcing_windstress_tend, &
             ocn_vel_forcing_windstress_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: windStressOn
   real (kind=RKIND) :: rho_ref


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_forcing_windstress_tend
!
!> \brief   Computes tendency term from wind stress
!> \author  Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the wind stress tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_windstress_tend(grid, u_src, h_edge, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         u_src    !< Input: wind stress

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

      if(.not.windStressOn) return

      nEdgesSolve = grid % nEdgesSolve
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array

      do iEdge=1,nEdgesSolve

        k = maxLevelEdgeTop(iEdge)

        ! efficiency note: it would be nice to avoid this
        ! if within a do.  This could be done with
        ! k =  max(maxLevelEdgeTop(iEdge),1)
        ! and then tend_u(1,iEdge) is just not used for land cells.

        if (k>0) then
           ! forcing in top layer only
           tend(1,iEdge) =  tend(1,iEdge) &
              + u_src(1,iEdge)/rho_ref/h_edge(1,iEdge)
        endif

      enddo


   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_windstress_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_forcing_windstress_init
!
!> \brief   Initializes ocean wind stress forcing
!> \author  Doug Jacobsen
!> \date    16 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes quantities related to wind stress 
!>  in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_windstress_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag


      windStressOn = .true.
      rho_ref = 1000.0

      err = 0


   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_windstress_init!}}}

!***********************************************************************

end module ocn_vel_forcing_windstress

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
