!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_forcing_rayleigh
!
!> \brief MPAS ocean Rayleigh Friction (to be used to smooth "shocks" from cold starts)
!> \author Todd Ringler
!> \date   5 January 2012
!> \version SVN:$Id:$
!> \details
!>  This module contains the routine for computing 
!>  tendencies based on linear Rayleigh friction.
!
!-----------------------------------------------------------------------

module ocn_vel_forcing_rayleigh

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

   public :: ocn_vel_forcing_rayleigh_tend, &
             ocn_vel_forcing_rayleigh_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: rayleighFrictionOn
   real (kind=RKIND) :: rayleighDampingCoef


!***********************************************************************

contains

!***********************************************************************
!
!  ocn_vel_forcing_rayleigh_tend
!
!> \brief   Computes tendency term from Rayleigh friction
!> \author  Todd Ringler
!> \date    5 January 2012
!> \version SVN:$Id$
!> \details 
!>  This routine computes the Rayleigh friction tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_rayleigh_tend(grid, normalVelocity, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity    !< Input: velocity 

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

      if(.not.rayleighFrictionOn) return

      nEdgesSolve = grid % nEdgesSolve
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array

      do iEdge=1,nEdgesSolve
        do k=1,maxLevelEdgeTop(iEdge)

           tend(k,iEdge) = tend(k,iEdge) - rayleighDampingCoef * normalVelocity(k,iEdge)

        enddo
      enddo


   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_rayleigh_tend!}}}

!***********************************************************************
!
!  ocn_vel_forcing_rayleigh_init
!
!> \brief   Initializes ocean Rayleigh friction
!> \author  Todd Ringler
!> \date    5 January 2012
!> \version SVN:$Id$
!> \details 
!>  This routine initializes quantities related to 
!>  in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_rayleigh_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag


      err = 0

      rayleighDampingCoef = 0.0

      if (config_rayleigh_friction) then
          rayleighFrictionOn = .true.
          rayleighDampingCoef = config_rayleigh_damping_coeff
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_rayleigh_init!}}}

!***********************************************************************

end module ocn_vel_forcing_rayleigh

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
