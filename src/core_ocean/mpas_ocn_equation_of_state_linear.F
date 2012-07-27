!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_equation_of_state_linear
!
!> \brief MPAS ocean equation of state driver
!> \author Doug Jacobsen
!> \date   28 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the main driver routine for calling
!>  the equation of state.
!
!-----------------------------------------------------------------------

module ocn_equation_of_state_linear

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

   public :: ocn_equation_of_state_linear_rho, &
             ocn_equation_of_state_linear_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_equation_of_state_linear_rho
!
!> \brief   Calls equation of state
!> \author  Doug Jacobsen
!> \date    28 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine uses a linear equation of state to update the density
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_linear_rho(grid, indexT, indexS, tracers, rho, err)!{{{
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

      type (mesh_type), intent(in) :: grid
      real (kind=RKIND), dimension(:,:), intent(inout) :: rho
      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers
      integer, intent(in) :: indexT, indexS
      integer, intent(out) :: err

      real (kind=RKIND), parameter :: rho_ref = 1025.022 ! kg / m^3
      real (kind=RKIND), parameter :: alpha =  2.55e-1 ! kg / m^3 / K (dT/dRho)
      real (kind=RKIND), parameter :: beta = 7.64e-1 ! kg / m^3 / psu (dS/dRho)
      real (kind=RKIND), parameter :: T_ref = 19.0 ! K
      real (kind=RKIND), parameter :: S_ref = 35.0 ! psu
      real (kind=RKIND), parameter :: rho_prime_ref = rho_ref + alpha * T_ref - beta * S_ref

      integer, dimension(:), pointer :: maxLevelCell
      integer :: nCells, iCell, k
      type (dm_info) :: dminfo

      maxLevelCell      => grid % maxLevelCell % array
      nCells      = grid % nCells

      err = 0

      do iCell=1,nCells
         do k=1,maxLevelCell(iCell)
            ! Linear equation of state
            ! rho = rho_ref - alpha * (T - T_ref) + beta * (S - S_ref)
            rho(k,iCell) = rho_prime_ref - alpha*tracers(indexT,k,iCell) + beta*tracers(indexS,k,iCell)
         end do
      end do

   end subroutine ocn_equation_of_state_linear_rho!}}}

!***********************************************************************
!
!  routine ocn_equation_of_state_linear_init
!
!> \brief   Initializes ocean momentum horizontal mixing quantities
!> \author  Doug Jacobsen
!> \date    28 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  horizontal velocity mixing in the ocean. Since a variety of 
!>  parameterizations are available, this routine primarily calls the
!>  individual init routines for each parameterization. 
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_linear_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err

      integer :: err1, err2

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_equation_of_state_linear_init!}}}

!***********************************************************************

end module ocn_equation_of_state_linear

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
