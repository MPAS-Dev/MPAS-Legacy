!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_rk4
!
!> \brief MPAS ocean RK4 Time integration scheme
!> \author Doug Jacobsen
!> \date   26 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the RK4 time integration routine.
!
!-----------------------------------------------------------------------

module ocn_time_integration_rk4

   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_dmpar
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer

   use ocn_tendency
   use ocn_diagnostics

   use ocn_equation_of_state
   use ocn_vmix
   use ocn_time_average

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

   public :: ocn_time_integrator_rk4

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integrator_rk4
!
!> \brief MPAS ocean RK4 Time integration scheme
!> \author Doug Jacobsen
!> \date   26 September 2011
!> \version SVN:$Id:$
!> \details
!>  This routine integrates one timestep (dt) using an RK4 time integrator.
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integrator_rk4(domain, dt)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Advance model state forward in time by the specified time step using 
   !   4th order Runge-Kutta
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
   !                 plus grid meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information
      real (kind=RKIND), intent(in) :: dt !< Input: timestep

      integer :: iCell, k, i, err
      type (block_type), pointer :: block

      integer :: rk_step

      real (kind=RKIND), dimension(4) :: rk_weights, rk_substep_weights

      integer :: nCells, nEdges, nVertLevels, num_tracers
      real (kind=RKIND) :: coef
      real (kind=RKIND), dimension(:,:), pointer :: &
        u, layerThickness, layerThicknessEdge, vertViscTopOfEdge, vertDiffTopOfCell, kineticEnergyEdge
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      integer, dimension(:), pointer :: & 
        maxLevelCell, maxLevelEdgeTop
      real (kind=RKIND), dimension(:), allocatable:: A,C,uTemp
      real (kind=RKIND), dimension(:,:), allocatable:: tracersTemp

      call mpas_setup_provis_states(domain % blocklist)

      !
      ! Initialize time_levs(2) with state at current time
      ! Initialize first RK state
      ! Couple tracers time_levs(2) with layerThickness in time-levels
      ! Initialize RK weights
      !
      block => domain % blocklist
      do while (associated(block))
        block % state % time_levs(2) % state % normalVelocity % array(:,:) = block % state % time_levs(1) % state % normalVelocity % array(:,:)
        block % state % time_levs(2) % state % layerThickness % array(:,:) = block % state % time_levs(1) % state % layerThickness % array(:,:)
        do iCell=1,block % mesh % nCells  ! couple tracers to h
          do k=1,block % mesh % maxLevelCell % array(iCell)
            block % state % time_levs(2) % state % tracers % array(:,k,iCell) = block % state % time_levs(1) % state % tracers % array(:,k,iCell) &
                                                                      * block % state % time_levs(1) % state % layerThickness % array(k,iCell)
           end do
        end do

        call mpas_copy_state(block % provis, block % state % time_levs(1) % state)

        block => block % next
      end do

      rk_weights(1) = dt/6.
      rk_weights(2) = dt/3.
      rk_weights(3) = dt/3.
      rk_weights(4) = dt/6.

      rk_substep_weights(1) = dt/2.
      rk_substep_weights(2) = dt/2.
      rk_substep_weights(3) = dt
      rk_substep_weights(4) = 0.


      call mpas_timer_start("RK4-main loop")
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN RK loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do rk_step = 1, 4
! ---  update halos for diagnostic variables

        call mpas_timer_start("RK4-diagnostic halo update")
        call mpas_dmpar_exch_halo_field(domain % blocklist % provis % potentialVorticityEdge)
        if (config_mom_del4 > 0.0) then
           call mpas_dmpar_exch_halo_field(domain % blocklist % provis % divergence)
           call mpas_dmpar_exch_halo_field(domain % blocklist % provis % relativeVorticity)
        end if
        call mpas_timer_stop("RK4-diagnostic halo update")

! ---  compute tendencies

        call mpas_timer_start("RK4-tendency computations")
        block => domain % blocklist
        do while (associated(block))
           ! advection of u uses u, while advection of layerThickness and tracers use uTransport.
           call ocn_vert_transport_velocity_top(block % mesh, block % provis % layerThickness % array, block % provis % layerThicknessEdge % array, &
              block % provis % normalVelocity % array, block % provis % vertTransportVelocityTop % array, err)
           call ocn_tend_vel(block % tend, block % provis, block % diagnostics, block % mesh)

           call ocn_vert_transport_velocity_top(block % mesh, block % provis % layerThickness % array, block % provis % layerThicknessEdge % array, &
              block % provis % uTransport % array, block % provis % vertTransportVelocityTop % array, err)
           call ocn_tend_thick(block % tend, block % provis, block % mesh)

           if (config_filter_btr_mode) then
               call ocn_filter_btr_mode_tend_vel(block % tend, block % provis, block % mesh)
           endif

           call ocn_tend_tracer(block % tend, block % provis, block % diagnostics, block % mesh, dt)
           block => block % next
        end do
        call mpas_timer_stop("RK4-tendency computations")

! ---  update halos for prognostic variables

        call mpas_timer_start("RK4-pronostic halo update")
        call mpas_dmpar_exch_halo_field(domain % blocklist % tend % normalVelocity)
        call mpas_dmpar_exch_halo_field(domain % blocklist % tend % layerThickness)
        call mpas_dmpar_exch_halo_field(domain % blocklist % tend % tracers)
        call mpas_timer_stop("RK4-pronostic halo update")

! ---  compute next substep state

        call mpas_timer_start("RK4-update diagnostic variables")
        if (rk_step < 4) then
           block => domain % blocklist
           do while (associated(block))

              block % provis % normalVelocity % array(:,:) = block % state % time_levs(1) % state % normalVelocity % array(:,:)  &
                                                    + rk_substep_weights(rk_step) * block % tend % normalVelocity % array(:,:)

              block % provis % layerThickness % array(:,:) = block % state % time_levs(1) % state % layerThickness % array(:,:)  &
                                              + rk_substep_weights(rk_step) * block % tend % layerThickness % array(:,:)
              do iCell=1,block % mesh % nCells
                 do k=1,block % mesh % maxLevelCell % array(iCell)
                 block % provis % tracers % array(:,k,iCell) = ( block % state % time_levs(1) % state % layerThickness % array(k,iCell) * &
                                                                 block % state % time_levs(1) % state % tracers % array(:,k,iCell)  &
                                                             + rk_substep_weights(rk_step) * block % tend % tracers % array(:,k,iCell) &
                                                               ) / block % provis % layerThickness % array(k,iCell)
                 end do

              end do
              if (config_prescribe_velocity) then
                 block % provis % normalVelocity % array(:,:) = block % state % time_levs(1) % state % normalVelocity % array(:,:)
              end if

              if (config_prescribe_thickness) then
                 block % provis % layerThickness % array(:,:) = block % state % time_levs(1) % state % layerThickness % array(:,:)
              end if

              call ocn_diagnostic_solve(dt, block % provis, block % mesh)

              ! Compute velocity transport, used in advection terms of layerThickness and tracer tendency
              block % provis % uTransport % array(:,:) &
                    = block % provis % normalVelocity % array(:,:) &
                    + block % provis % uBolusGM   % array(:,:)

              block => block % next
           end do
        end if
        call mpas_timer_stop("RK4-update diagnostic variables")

!--- accumulate update (for RK4)

        call mpas_timer_start("RK4-RK4 accumulate update")
        block => domain % blocklist
        do while (associated(block))
           block % state % time_levs(2) % state % normalVelocity % array(:,:) = block % state % time_levs(2) % state % normalVelocity % array(:,:) &
                                   + rk_weights(rk_step) * block % tend % normalVelocity % array(:,:) 

           block % state % time_levs(2) % state % layerThickness % array(:,:) = block % state % time_levs(2) % state % layerThickness % array(:,:) &
                                   + rk_weights(rk_step) * block % tend % layerThickness % array(:,:) 

           do iCell=1,block % mesh % nCells
              do k=1,block % mesh % maxLevelCell % array(iCell)
                 block % state % time_levs(2) % state % tracers % array(:,k,iCell) =  &
                                                                        block % state % time_levs(2) % state % tracers % array(:,k,iCell) &
                                                                        + rk_weights(rk_step) * block % tend % tracers % array(:,k,iCell)
              end do
           end do

           block => block % next
        end do
        call mpas_timer_stop("RK4-RK4 accumulate update")

      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END RK loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call mpas_timer_stop("RK4-main loop")

      !
      !  A little clean up at the end: rescale tracer fields and compute diagnostics for new state
      !
      call mpas_timer_start("RK4-cleaup phase")

      ! Rescale tracers
      block => domain % blocklist
      do while(associated(block))
        do iCell = 1, block % mesh % nCells
          do k = 1, block % mesh % maxLevelCell % array(iCell)
            block % state % time_levs(2) % state % tracers % array(:, k, iCell) = block % state % time_levs(2) % state % tracers % array(:, k, iCell) &
                                                                                / block % state % time_levs(2) % state % layerThickness % array(k, iCell)
          end do
        end do
        block => block % next
      end do

      call mpas_timer_start("RK4-implicit vert mix")
      block => domain % blocklist
      do while(associated(block))

        ! Call ocean diagnostic solve in preparation for vertical mixing.  Note 
        ! it is called again after vertical mixing, because u and tracers change.
        ! For Richardson vertical mixing, only density, layerThicknessEdge, and kineticEnergyEdge need to 
        ! be computed.  For kpp, more variables may be needed.  Either way, this
        ! could be made more efficient by only computing what is needed for the
        ! implicit vmix routine that follows. 
        call ocn_diagnostic_solve(dt, block % state % time_levs(2) % state, block % mesh)

        call ocn_vmix_implicit(dt, block % mesh, block % diagnostics, block % state % time_levs(2) % state, err)
        block => block % next
      end do

      ! Update halo on u and tracers, which were just updated for implicit vertical mixing.  If not done, 
      ! this leads to lack of volume conservation.  It is required because halo updates in RK4 are only
      ! conducted on tendencies, not on the velocity and tracer fields.  So this update is required to 
      ! communicate the change due to implicit vertical mixing across the boundary.
      call mpas_timer_start("RK4-implicit vert mix halos")
      call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % normalVelocity)
      call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % tracers)
      call mpas_timer_stop("RK4-implicit vert mix halos")

      call mpas_timer_stop("RK4-implicit vert mix")

      block => domain % blocklist
      do while (associated(block))
         if (config_prescribe_velocity) then
            block % state % time_levs(2) % state % normalVelocity % array(:,:) = block % state % time_levs(1) % state % normalVelocity % array(:,:)
         end if

         if (config_prescribe_thickness) then
            block % state % time_levs(2) % state % layerThickness % array(:,:) = block % state % time_levs(1) % state % layerThickness % array(:,:)
         end if

         call ocn_diagnostic_solve(dt, block % state % time_levs(2) % state, block % mesh)

         ! Compute velocity transport, used in advection terms of layerThickness and tracer tendency
            block % state % time_levs(2) % state % uTransport % array(:,:) &
          = block % state % time_levs(2) % state % normalVelocity % array(:,:) &
          + block % state % time_levs(2) % state % uBolusGM % array(:,:)

         call mpas_reconstruct(block % mesh, block % state % time_levs(2) % state % normalVelocity % array,          &
                          block % state % time_levs(2) % state % normalVelocityX % array,            &
                          block % state % time_levs(2) % state % normalVelocityY % array,            &
                          block % state % time_levs(2) % state % normalVelocityZ % array,            &
                          block % state % time_levs(2) % state % normalVelocityZonal % array,        &
                          block % state % time_levs(2) % state % normalVelocityMeridional % array    &
                         )

         call mpas_reconstruct(block % mesh, block % mesh % normalVelocityForcing % array,          &
                          block % state % time_levs(2) % state % normalVelocityForcingReconstructX % array,            &
                          block % state % time_levs(2) % state % normalVelocityForcingReconstructY % array,            &
                          block % state % time_levs(2) % state % normalVelocityForcingReconstructZ % array,            &
                          block % state % time_levs(2) % state % normalVelocityForcingReconstructZonal % array,        &
                          block % state % time_levs(2) % state % normalVelocityForcingReconstructMeridional % array    &
                         )

         call ocn_time_average_accumulate(block % state % time_levs(2) % state, block % state % time_levs(1) % state)

         block => block % next
      end do
      call mpas_timer_stop("RK4-cleaup phase")

      call mpas_deallocate_provis_states(domain % blocklist)

   end subroutine ocn_time_integrator_rk4!}}}

end module ocn_time_integration_rk4

! vim: foldmethod=marker
