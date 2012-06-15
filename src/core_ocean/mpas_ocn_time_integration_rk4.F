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
      type (state_type), target :: provis
      type (state_type), pointer :: provis_ptr

      integer :: rk_step

      real (kind=RKIND), dimension(4) :: rk_weights, rk_substep_weights

      integer :: nCells
      real (kind=RKIND), dimension(:,:), pointer :: u, h, h_edge, vertViscTopOfEdge, vertDiffTopOfCell, ke_edge
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      integer, dimension(:), pointer :: maxLevelCell

      block => domain % blocklist
      call mpas_allocate_state(block, provis, &
                          block % mesh % nCells, block % mesh % nEdges, block % mesh % maxEdges, block % mesh % maxEdges2, &
                          block % mesh % nVertices, block % mesh % vertexDegree, block % mesh % nVertLevels, block % mesh % nMonths )

      provis_ptr => provis
      call mpas_create_state_links(provis_ptr)

      !
      ! Initialize time_levs(2) with state at current time
      ! Initialize first RK state
      ! Couple tracers time_levs(2) with h in time-levels
      ! Initialize RK weights
      !
      block => domain % blocklist
      do while (associated(block))

         block % state % time_levs(2) % state % u % array(:,:) = block % state % time_levs(1) % state % u % array(:,:)
         block % state % time_levs(2) % state % h % array(:,:) = block % state % time_levs(1) % state % h % array(:,:)
         do iCell=1,block % mesh % nCells  ! couple tracers to h
           do k=1,block % mesh % maxLevelCell % array(iCell)
             block % state % time_levs(2) % state % tracers % array(:,k,iCell) = block % state % time_levs(1) % state % tracers % array(:,k,iCell) &
                                                                       * block % state % time_levs(1) % state % h % array(k,iCell)
            end do
         end do

         call mpas_copy_state(provis, block % state % time_levs(1) % state)

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
        call mpas_dmpar_exch_halo_field(provis % Vor_edge)
        if (config_h_mom_eddy_visc4 > 0.0) then
           call mpas_dmpar_exch_halo_field(provis % divergence)
           call mpas_dmpar_exch_halo_field(provis % vorticity)
        end if
        call mpas_timer_stop("RK4-diagnostic halo update")

! ---  compute tendencies

        call mpas_timer_start("RK4-tendency computations")
        block => domain % blocklist
        do while (associated(block))

           ! mrp 111206 put ocn_wtop call at top for ALE
           call ocn_wtop(provis, provis, block % mesh)

           if (.not.config_implicit_vertical_mix) then
              call ocn_vmix_coefs(block % mesh, provis, block % diagnostics, err)
           end if
           call ocn_tend_h(block % tend, provis, block % mesh)
           call ocn_tend_u(block % tend, provis, block % diagnostics, block % mesh)

           ! mrp 110718 filter btr mode out of u_tend
           ! still got h perturbations with just this alone.  Try to set uBtr=0 after full u computation
           if (config_rk_filter_btr_mode) then
               call ocn_filter_btr_mode_tend_u(block % tend, provis, block % mesh)
           endif

           call ocn_tend_scalar(block % tend, provis, block % diagnostics, block % mesh, dt)
           block => block % next
        end do
        call mpas_timer_stop("RK4-tendency computations")

! ---  update halos for prognostic variables

        call mpas_timer_start("RK4-pronostic halo update")
        call mpas_dmpar_exch_halo_field(domain % blocklist % tend % u)
        call mpas_dmpar_exch_halo_field(domain % blocklist % tend % h)
        call mpas_dmpar_exch_halo_field(domain % blocklist % tend % tracers)
        call mpas_timer_stop("RK4-pronostic halo update")

! ---  compute next substep state

        call mpas_timer_start("RK4-update diagnostic variables")
        if (rk_step < 4) then
           block => domain % blocklist
           do while (associated(block))

              provis % u % array(:,:)       = block % state % time_levs(1) % state % u % array(:,:)  &
                                         + rk_substep_weights(rk_step) * block % tend % u % array(:,:)

              provis % h % array(:,:)       = block % state % time_levs(1) % state % h % array(:,:)  &
                                         + rk_substep_weights(rk_step) * block % tend % h % array(:,:)
              do iCell=1,block % mesh % nCells
                 do k=1,block % mesh % maxLevelCell % array(iCell)
                    provis % tracers % array(:,k,iCell) = ( &
                                                                      block % state % time_levs(1) % state % h % array(k,iCell) * &
                                                                      block % state % time_levs(1) % state % tracers % array(:,k,iCell)  &
                                      + rk_substep_weights(rk_step) * block % tend % tracers % array(:,k,iCell) &
                                                                     ) / provis % h % array(k,iCell)
                 end do

              end do
              if (config_test_case == 1) then    ! For case 1, wind field should be fixed
                 provis % u % array(:,:) = block % state % time_levs(1) % state % u % array(:,:)
              end if

              if (config_prescribe_velocity) then
                 provis % u % array(:,:) = block % state % time_levs(1) % state % u % array(:,:)
              end if

              if (config_prescribe_thickness) then
                 provis % h % array(:,:) = block % state % time_levs(1) % state % h % array(:,:)
              end if

              call ocn_diagnostic_solve(dt, provis, block % mesh)

              ! Compute velocity transport, used in advection terms of h and tracer tendancy
                 provis % uTransport % array(:,:) &
               = provis % u          % array(:,:) &
               + provis % uBolusGM   % array(:,:)

              block => block % next
           end do
        end if
        call mpas_timer_stop("RK4-update diagnostic variables")



!--- accumulate update (for RK4)

        call mpas_timer_start("RK4-RK4 accumulate update")
        block => domain % blocklist
        do while (associated(block))
           block % state % time_levs(2) % state % u % array(:,:) = block % state % time_levs(2) % state % u % array(:,:) &
                                   + rk_weights(rk_step) * block % tend % u % array(:,:) 

           block % state % time_levs(2) % state % h % array(:,:) = block % state % time_levs(2) % state % h % array(:,:) &
                                   + rk_weights(rk_step) * block % tend % h % array(:,:) 

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
      !  A little clean up at the end: decouple new scalar fields and compute diagnostics for new state
      !
      call mpas_timer_start("RK4-cleaup phase")
      block => domain % blocklist
      do while (associated(block))

         u           => block % state % time_levs(2) % state % u % array
         tracers     => block % state % time_levs(2) % state % tracers % array
         h           => block % state % time_levs(2) % state % h % array
         h_edge      => block % state % time_levs(2) % state % h_edge % array
         ke_edge     => block % state % time_levs(2) % state % ke_edge % array
         vertViscTopOfEdge => block % diagnostics % vertViscTopOfEdge % array
         vertDiffTopOfCell => block % diagnostics % vertDiffTopOfCell % array
         maxLevelCell    => block % mesh % maxLevelCell % array
                  
         nCells      = block % mesh % nCells

         do iCell=1,nCells
            do k=1,maxLevelCell(iCell)
               tracers(:,k,iCell) = tracers(:,k,iCell) / h(k,iCell)
            end do
         end do

         if (config_implicit_vertical_mix) then
            call mpas_timer_start("RK4-implicit vert mix")

            call ocn_vmix_coefs(block % mesh, block % state % time_levs(2) % state, block % diagnostics, err)

            !
            !  Implicit vertical solve for momentum
            !
            call ocn_vel_vmix_tend_implicit(block % mesh, dt, ke_edge, vertViscTopOfEdge, h, h_edge, u, err)

          !  mrp 110718 filter btr mode out of u
           if (config_rk_filter_btr_mode) then
               call ocn_filter_btr_mode_u(block % state % time_levs(2) % state, block % mesh)
               !block % tend % h % array(:,:) = 0.0 ! I should not need this
           endif

            !
            !  Implicit vertical solve for tracers
            !

            call ocn_tracer_vmix_tend_implicit(block % mesh, dt, vertDiffTopOfCell, h, tracers, err)
            call mpas_timer_stop("RK4-implicit vert mix")
         end if

         block => block % next
      end do

      ! Update halo on u and tracers, which weres just updated for implicit vertical mixing.  If not done, 
      ! this leads to lack of volume conservation.  It is required because halo updates in RK4 are only
      ! conducted on tendencies, not on the velocity and tracer fields.  So this update is required to 
      ! communicate the change due to implicit vertical mixing across the boundary.

      if (config_implicit_vertical_mix) then
         call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % u)
         call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % tracers)
      end if

      block => domain % blocklist
      do while (associated(block))

         if (config_test_case == 1) then    ! For case 1, wind field should be fixed
            block % state % time_levs(2) % state % u % array(:,:) = block % state % time_levs(1) % state % u % array(:,:)
         end if

         if (config_prescribe_velocity) then
            block % state % time_levs(2) % state % u % array(:,:) = block % state % time_levs(1) % state % u % array(:,:)
         end if

         if (config_prescribe_thickness) then
            block % state % time_levs(2) % state % h % array(:,:) = block % state % time_levs(1) % state % h % array(:,:)
         end if

         call ocn_diagnostic_solve(dt, block % state % time_levs(2) % state, block % mesh)

         ! Compute velocity transport, used in advection terms of h and tracer tendancy
            block % state % time_levs(2) % state % uTransport % array(:,:) &
          = block % state % time_levs(2) % state % u % array(:,:) &
          + block % state % time_levs(2) % state % uBolusGM % array(:,:)

         call mpas_reconstruct(block % mesh, block % state % time_levs(2) % state % u % array,          &
                          block % state % time_levs(2) % state % uReconstructX % array,            &
                          block % state % time_levs(2) % state % uReconstructY % array,            &
                          block % state % time_levs(2) % state % uReconstructZ % array,            &
                          block % state % time_levs(2) % state % uReconstructZonal % array,        &
                          block % state % time_levs(2) % state % uReconstructMeridional % array    &
                         )

!TDR
         call mpas_reconstruct(block % mesh, block % mesh % u_src % array,          &
                          block % state % time_levs(2) % state % uSrcReconstructX % array,            &
                          block % state % time_levs(2) % state % uSrcReconstructY % array,            &
                          block % state % time_levs(2) % state % uSrcReconstructZ % array,            &
                          block % state % time_levs(2) % state % uSrcReconstructZonal % array,        &
                          block % state % time_levs(2) % state % uSrcReconstructMeridional % array    &
                         )
!TDR

         call ocn_time_average_accumulate(block % state % time_levs(2) % state, block % state % time_levs(1) % state)

         block => block % next
      end do
      call mpas_timer_stop("RK4-cleaup phase")

      call mpas_deallocate_state(provis)

   end subroutine ocn_time_integrator_rk4!}}}

end module ocn_time_integration_rk4

! vim: foldmethod=marker
