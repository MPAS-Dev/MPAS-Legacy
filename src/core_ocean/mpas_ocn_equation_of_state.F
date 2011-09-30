!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_equation_of_state
!
!> \brief MPAS ocean equation of state driver
!> \author Doug Jacobsen
!> \date   29 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the main driver routine for calling
!>  the equation of state.
!
!-----------------------------------------------------------------------

module ocn_equation_of_state

   use grid_types
   use configure
   use timer
   use ocn_equation_of_state_linear
   use ocn_equation_of_state_jm

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

   public :: ocn_equation_of_state_rho, &
             ocn_equation_of_state_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: eosON
   logical :: linearEos, jmEos


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_equation_of_state
!
!> \brief   Calls equation of state
!> \author  Doug Jacobsen
!> \date    29 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine calls the equation of state to update the density
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_rho(s, grid, k_displaced, displacement_type, err)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This module contains routines necessary for computing the density
   !  from model temperature and salinity using an equation of state.
   !
   ! Input: grid - grid metadata
   !        s - state: tracers
   !        k_displaced 
   !  If k_displaced<=0, state % rho is returned with no displaced
   !  If k_displaced>0,the state % rhoDisplaced is returned, and is for
   !  a parcel adiabatically displaced from its original level to level 
   !  k_displaced.  This does not effect the linear EOS.
   !
   ! Output: s - state: computed density
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none

      type (state_type), intent(inout) :: s
      type (mesh_type), intent(in) :: grid
      integer, intent(out) :: err
      integer :: k_displaced
      character(len=8), intent(in) :: displacement_type

      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:,:), pointer :: rho
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      integer :: nCells, iCell, k, indexT, indexS
      type (dm_info) :: dminfo

      err = 0

      if(.not.eosOn) return

      call timer_start("ocn_equation_of_state_rho")

      tracers => s % tracers % array
      indexT = s % index_temperature
      indexS = s % index_salinity

      if (linearEos) then
         rho => s % rho % array

         call ocn_equation_of_state_linear_rho(grid, indexT, indexS, tracers, rho, err)

      elseif (jmEos) then

         if(k_displaced == 0) then
             rho => s % rho % array
         else
             rho => s % rhoDisplaced % array
         endif

         call ocn_equation_of_state_jm_rho(grid, k_displaced, displacement_type, indexT, indexS, tracers, rho, err)

      endif

      call timer_stop("ocn_equation_of_state_rho")

   end subroutine ocn_equation_of_state_rho!}}}

!***********************************************************************
!
!  routine ocn_equation_of_stateInit
!
!> \brief   Initializes ocean momentum horizontal mixing quantities
!> \author  Doug Jacobsen
!> \date    29 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  horizontal velocity mixing in the ocean. Since a variety of 
!>  parameterizations are available, this routine primarily calls the
!>  individual init routines for each parameterization. 
!
!----------------------------------------------------------------------

   subroutine ocn_equation_of_state_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err

      err = 0
      eosON = .false.
      linearEos = .false.
      jmEos = .false.

      if(config_vert_grid_type.eq.'zlevel') then
          eosON = .true.

          if (config_eos_type.eq.'linear') then
              linearEos = .true.
          elseif (config_eos_type.eq.'jm') then
              jmEos = .true.
          else
              print *,'Invalid choice for config_eos_type.'
              print *,'  Choices are: linear, jm'
              err = 1
          endif
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_equation_of_state_init!}}}

!***********************************************************************

end module ocn_equation_of_state

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker