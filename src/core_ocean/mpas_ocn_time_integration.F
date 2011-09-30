!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration
!
!> \brief MPAS ocean time integration driver
!> \author Doug Jacobsen
!> \date   26 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the main driver routine for calling
!>  the time integration scheme
!
!-----------------------------------------------------------------------

module ocn_time_integration

   use grid_types
   use configure
   use constants
   use dmpar
   use vector_reconstruction
   use spline_interpolation
   use timer

   use ocn_time_integration_rk4
   use ocn_time_integration_split

   implicit none
   private
   save

   public :: ocn_timestep, &
             ocn_timestep_init

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

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

    logical :: rk4On, splitOn

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_timestep
!
!> \brief MPAS ocean time integration driver
!> \author Doug Jacobsen
!> \date   26 September 2011
!> \version SVN:$Id:$
!> \details
!>  This routine handles a single timestep for the ocean. It determines
!>  the time integrator that will be used for the run, and calls the
!>  appropriate one.
!
!-----------------------------------------------------------------------

   subroutine ocn_timestep(domain, dt, timeStamp)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Advance model state forward in time by the specified time step
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
   !                 plus grid meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(in) :: timeStamp

      type (dm_info) :: dminfo
      type (block_type), pointer :: block

      if (rk4On) then
         call ocn_time_integrator_rk4(domain, dt)
      elseif (splitOn) then
         call ocn_time_integrator_split(domain, dt)
     endif

     block => domain % blocklist
     do while (associated(block))
        block % state % time_levs(2) % state % xtime % scalar = timeStamp

        if (isNaN(sum(block % state % time_levs(2) % state % u % array))) then
           write(0,*) 'Abort: NaN detected'
           call dmpar_abort(dminfo)
        endif

        block => block % next
     end do

   end subroutine ocn_timestep!}}}

   subroutine ocn_timestep_init(err)!{{{

      integer, intent(out) :: err

      rk4On = .false.
      splitOn = .false.

      if (trim(config_time_integration) == 'RK4') then
          rk4On = .true.
      elseif (trim(config_time_integration) == 'split_explicit' &
          .or.trim(config_time_integration) == 'unsplit_explicit') then
          splitOn = .true.
      else
          err = 1
          write(*,*) 'Incorrect choice for config_time_integration:', trim(config_time_integration)
          write(*,*) '   choices are: RK4, split_explicit, unsplit_explicit'
      endif


   end subroutine ocn_timestep_init!}}}

end module ocn_time_integration

! vim: foldmethod=marker
