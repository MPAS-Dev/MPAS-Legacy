module mpas_core

   use mpas_configure
   use mpas_framework
   use mpas_timekeeping
   use mpas_dmpar
   use mpas_timer

   use ocn_global_diagnostics
   use ocn_test_cases
   use ocn_time_integration
   use ocn_tendency

   use ocn_vel_pressure_grad
   use ocn_vel_vadv
   use ocn_vel_hmix
   use ocn_vel_forcing

   use ocn_tracer_hadv
   use ocn_tracer_vadv
   use ocn_tracer_hmix
   use ocn_gm
   use ocn_restoring

   use ocn_equation_of_state

   use ocn_vmix

   use ocn_time_average

   type (io_output_object) :: restart_obj

   integer :: current_outfile_frames

   type (MPAS_Clock_type) :: clock

   integer, parameter :: outputAlarmID = 1
   integer, parameter :: restartAlarmID = 2
   integer, parameter :: statsAlarmID = 3

   type (timer_node), pointer :: globalDiagTimer, timeIntTimer
   type (timer_node), pointer :: initDiagSolveTimer

   contains

   subroutine mpas_core_init(domain, startTimeStamp)!{{{

      use mpas_grid_types
      use mpas_ocn_tracer_advection

      implicit none

      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp

      real (kind=RKIND) :: dt
      type (block_type), pointer :: block
      type (dm_info) :: dminfo

      integer :: err, err_tmp

      ! Initialize submodules before initializing blocks.
      call ocn_timestep_init(err)

      call ocn_vel_pressure_grad_init(err_tmp)
      err = ior(err, err_tmp)

      call ocn_vel_vadv_init(err_tmp)
      err = ior(err, err_tmp)
      call ocn_vel_hmix_init(err_tmp)
      err = ior(err, err_tmp)
      call ocn_vel_forcing_init(err_tmp)
      err = ior(err, err_tmp)

      call ocn_tracer_hadv_init(err_tmp)
      err = ior(err, err_tmp)
      call ocn_tracer_vadv_init(err_tmp)
      err = ior(err, err_tmp)
      call ocn_tracer_hmix_init(err_tmp)
      err = ior(err, err_tmp)
      call ocn_restoring_init(err_tmp)
      err = ior(err, err_tmp)

      call ocn_vmix_init(err_tmp)
      err = ior(err, err_tmp)

      call ocn_equation_of_state_init(err_tmp)
      err = ior(err, err_tmp)

      call ocn_tendency_init(err_tmp)
      err = ior(err,err_tmp)

      call mpas_ocn_tracer_advection_init(err_tmp)
      err = ior(err,err_tmp)

      call mpas_timer_init(domain)

      if(err.eq.1) then
          call mpas_dmpar_abort(dminfo)
      endif

      if (.not. config_do_restart) call setup_sw_test_case(domain)

      call ocn_compute_max_level(domain)

      call ocn_init_z_level(domain)

      print *, ' Vertical grid type is: ',config_vert_grid_type

      if (config_vert_grid_type.ne.'isopycnal'.and. &
          config_vert_grid_type.ne.'zlevel'.and. &
          config_vert_grid_type.ne.'zstar1'.and. &
          config_vert_grid_type.ne.'zstar'.and. &
          config_vert_grid_type.ne.'zstarWeights') then
         print *, ' Incorrect choice of config_vert_grid_type.'
         call mpas_dmpar_abort(dminfo)
      endif

      print *, ' Pressure type is: ',config_pressure_type
      if (config_pressure_type.ne.'pressure'.and. &
          config_pressure_type.ne.'MontgomeryPotential') then
         print *, ' Incorrect choice of config_pressure_type.'
         call mpas_dmpar_abort(dminfo)
      endif

      if (config_filter_btr_mode.and. &
          config_vert_grid_type.ne.'zlevel')then
         print *, 'filter_btr_mode has only been tested with'// &
            ' config_vert_grid_type=zlevel.'
         call mpas_dmpar_abort(dminfo)
      endif

      !
      ! Initialize core
      !
      dt = config_dt

      call ocn_simulation_clock_init(domain, dt, startTimeStamp)

      block => domain % blocklist
      do while (associated(block))
         call mpas_init_block(block, block % mesh, dt, err)
         if(err.eq.1) then
             call mpas_dmpar_abort(dminfo)
         endif
         block % state % time_levs(1) % state % xtime % scalar = startTimeStamp 
         block => block % next
      end do

   ! mrp 100316 In order for this to work, we need to pass domain % dminfo as an 
   ! input arguement into mpas_init.  Ask about that later.  For now, there will be
   ! no initial statistics write.
   
      if (config_initial_stats) then
          call mpas_timer_start("global diagnostics", .false., globalDiagTimer)
          call ocn_compute_global_diagnostics(domain, 1 , 0, dt)
          call mpas_timer_stop("global diagnostics", globalDiagTimer)
!         call mpas_output_state_init(output_obj, domain, "OUTPUT")
!         call ocn_write_output_frame(output_obj, output_frame, domain)
      endif

      current_outfile_frames = 0

   end subroutine mpas_core_init!}}}

   subroutine ocn_simulation_clock_init(domain, dt, startTimeStamp)!{{{

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(out) :: startTimeStamp

      type (MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      integer :: ierr

      call mpas_set_time(curr_time=startTime, dateTimeString=config_start_time, ierr=ierr)
      call mpas_set_timeInterval(timeStep, dt=dt, ierr=ierr)

      if (trim(config_run_duration) /= "none") then
         call mpas_set_timeInterval(runDuration, timeString=config_run_duration, ierr=ierr)
         call mpas_create_clock(clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=ierr)

         if (trim(config_stop_time) /= "none") then
            call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=ierr)
            if(startTime + runduration /= stopTime) then
               write(0,*) 'Warning: config_run_duration and config_stop_time are inconsitent: using config_run_duration.'
            end if
         end if
      else if (trim(config_stop_time) /= "none") then
         call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=ierr)
         call mpas_create_clock(clock, startTime=startTime, timeStep=timeStep, stopTime=stopTime, ierr=ierr)
      else
          write(0,*) 'Error: Neither config_run_duration nor config_stop_time were specified.'
          call mpas_dmpar_finalize(domain % dminfo)
      end if

      ! set output alarm
      call mpas_set_timeInterval(alarmTimeStep, timeString=config_output_interval, ierr=ierr)
      alarmStartTime = startTime + alarmTimeStep
      call mpas_add_clock_alarm(clock, outputAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)

      ! set restart alarm, if necessary
      if (trim(config_restart_interval) /= "none") then
         call mpas_set_timeInterval(alarmTimeStep, timeString=config_restart_interval, ierr=ierr)
         alarmStartTime = startTime + alarmTimeStep
         call mpas_add_clock_alarm(clock, restartAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
      end if

      !TODO: use this code if we desire to convert config_stats_interval to alarms 
      !(must also change config_stats_interval type to character) 
      ! set stats alarm, if necessary
      !if (trim(config_stats_interval) /= "none") then      
      !   call mpas_set_timeInterval(alarmTimeStep, timeString=config_stats_interval, ierr=ierr)
      !   alarmStartTime = startTime + alarmTimeStep
      !   call mpas_add_clock_alarm(clock, statsAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
      !end if

      call mpas_get_time(curr_time=startTime, dateTimeString=startTimeStamp, ierr=ierr)

   end subroutine ocn_simulation_clock_init!}}}

   subroutine mpas_init_block(block, mesh, dt, err)!{{{
   
      use mpas_grid_types
      use mpas_rbf_interpolation
      use mpas_vector_reconstruction
      use mpas_ocn_tracer_advection
      use ocn_advection
   
      implicit none
   
      type (block_type), intent(inout) :: block
      type (mesh_type), intent(inout) :: mesh
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err
      integer :: i, iEdge, iCell, k
      integer :: err1
   
      call ocn_initialize_advection_rk(mesh, err)
      call mpas_ocn_tracer_advection_coefficients(mesh, err1)
      err = ior(err, err1)

      call ocn_time_average_init(block % state % time_levs(1) % state)
   
      call mpas_timer_start("diagnostic solve", .false., initDiagSolveTimer)
      call ocn_diagnostic_solve(dt, block % state % time_levs(1) % state, mesh)
      call mpas_timer_stop("diagnostic solve", initDiagSolveTimer)

      ! Compute velocity transport, used in advection terms of h and tracer tendancy
        block % state % time_levs(1) % state % uTransport % array(:,:) &
      = block % state % time_levs(1) % state % u % array(:,:) &
      + block % state % time_levs(1) % state % uBolusGM % array(:,:)

      call ocn_wtop(block % state % time_levs(1) % state,block % state % time_levs(1) % state, mesh)

      call ocn_compute_mesh_scaling(mesh)
 
      call mpas_rbf_interp_initialize(mesh)
      call mpas_init_reconstruct(mesh)
      call mpas_reconstruct(mesh, block % state % time_levs(1) % state % u % array,                  &
                       block % state % time_levs(1) % state % uReconstructX % array,            &
                       block % state % time_levs(1) % state % uReconstructY % array,            &
                       block % state % time_levs(1) % state % uReconstructZ % array,            &
                       block % state % time_levs(1) % state % uReconstructZonal % array,        &
                       block % state % time_levs(1) % state % uReconstructMeridional % array    &
                      )

      ! initialize velocities and tracers on land to be -1e34
      ! The reconstructed velocity on land will have values not exactly
      ! -1e34 due to the interpolation of reconstruction.

      block % mesh % areaCell % array(block % mesh % nCells+1) = -1.0e34

      do iEdge=1,block % mesh % nEdges
         ! mrp 101115 note: in order to include flux boundary conditions, the following
         ! line will need to change.  Right now, set boundary edges between land and 
         ! water to have zero velocity.
         block % state % time_levs(1) % state % u % array( &
             block % mesh % maxLevelEdgeTop % array(iEdge)+1 &
            :block % mesh % maxLevelEdgeBot % array(iEdge), iEdge) = 0.0

         block % state % time_levs(1) % state % u % array( &
             block % mesh % maxLevelEdgeBot % array(iEdge)+1: &
             block % mesh % nVertLevels,iEdge) = 0.0
! mrp changed to 0
!             block % mesh % nVertLevels,iEdge) = -1e34
      end do
      do iCell=1,block % mesh % nCells
         block % state % time_levs(1) % state % tracers % array( &
            :, block % mesh % maxLevelCell % array(iCell)+1 &
              :block % mesh % nVertLevels,iCell) =  0.0
! mrp changed to 0
!              :block % mesh % nVertLevels,iCell) =  -1e34
      end do

      do i=2,nTimeLevs
        call mpas_copy_state(block % state % time_levs(i) % state, &
                             block % state % time_levs(1) % state)
      end do

   end subroutine mpas_init_block!}}}
   
   subroutine mpas_core_run(domain, output_obj, output_frame)!{{{
   
      use mpas_grid_types
      use mpas_io_output
      use mpas_timer
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (io_output_object), intent(inout) :: output_obj
      integer, intent(inout) :: output_frame
   
      integer :: itimestep
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime
      character(len=32) :: timeStamp
      integer :: ierr
   
      ! Eventually, dt should be domain specific
      dt = config_dt

      currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
      write(0,*) 'Initial time ', timeStamp

      call ocn_write_output_frame(output_obj, output_frame, domain)
      block_ptr => domain % blocklist

      do while(associated(block_ptr))
        call ocn_time_average_init(block_ptr % state % time_levs(1) % state)
        block_ptr => block_ptr % next
      end do

      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      itimestep = 0
      do while (.not. mpas_is_clock_stop_time(clock))

         itimestep = itimestep + 1
         call mpas_advance_clock(clock)

         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
         write(0,*) 'Doing timestep ', timeStamp

         call mpas_timer_start("time integration", .false., timeIntTimer)
         call mpas_timestep(domain, itimestep, dt, timeStamp)
         call mpas_timer_stop("time integration", timeIntTimer)
   
         ! Move time level 2 fields back into time level 1 for next time step
         call mpas_shift_time_levels_state(domain % blocklist % state)
      
         if (mpas_is_alarm_ringing(clock, outputAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, outputAlarmID, ierr=ierr)
            ! output_frame will always be > 1 here unless it was reset after the maximum number of frames per outfile was reached
            if(output_frame == 1) then
               call mpas_output_state_finalize(output_obj, domain % dminfo)
               call mpas_output_state_init(output_obj, domain, "OUTPUT", trim(timeStamp))
            end if

            block_ptr => domain % blocklist
            do while (associated(block_ptr))
                call ocn_time_average_normalize(block_ptr % state % time_levs(1) % state)
                block_ptr => block_ptr % next
            end do

            call ocn_write_output_frame(output_obj, output_frame, domain)

            block_ptr => domain % blocklist
            do while (associated(block_ptr))
                call ocn_time_average_init(block_ptr % state % time_levs(1) % state)
                block_ptr => block_ptr % next
            end do
         end if

         if (mpas_is_alarm_ringing(clock, restartAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, restartAlarmID, ierr=ierr)

            ! Write one restart time per file
            call mpas_output_state_init(restart_obj, domain, "RESTART", trim(timeStamp))
            call mpas_output_state_for_domain(restart_obj, domain, 1)
            call mpas_output_state_finalize(restart_obj, domain % dminfo)
         end if

      end do

   end subroutine mpas_core_run!}}}
   
   subroutine ocn_write_output_frame(output_obj, output_frame, domain)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain and write model state to output file
   !
   ! Input/Output: domain - contains model state; diagnostic field are computed
   !                        before returning
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use mpas_grid_types
      use mpas_io_output
   
      implicit none
   
      integer, intent(inout) :: output_frame
      type (domain_type), intent(inout) :: domain
      type (io_output_object), intent(inout) :: output_obj
   
      integer :: i, j, k
      integer :: eoe
      type (block_type), pointer :: block_ptr
   
      block_ptr => domain % blocklist
      do while (associated(block_ptr))
         call ocn_compute_output_diagnostics(block_ptr % state % time_levs(1) % state, block_ptr % mesh)
         block_ptr => block_ptr % next
      end do
   
      call mpas_output_state_for_domain(output_obj, domain, output_frame)
      output_frame = output_frame + 1

      ! reset frame if the maximum number of frames per outfile has been reached
      if (config_frames_per_outfile > 0) then
         current_outfile_frames = current_outfile_frames + 1            
         if(current_outfile_frames >= config_frames_per_outfile) then
            current_outfile_frames = 0
            output_frame = 1
         end if
      end if
   
   end subroutine ocn_write_output_frame!}}}
   
   subroutine ocn_compute_output_diagnostics(state, grid)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain
   !
   ! Input: state - contains model prognostic fields
   !        grid  - contains grid metadata
   !
   ! Output: state - upon returning, diagnostic fields will have be computed
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use mpas_grid_types
   
      implicit none
   
      type (state_type), intent(inout) :: state
      type (mesh_type), intent(in) :: grid
   
      integer :: i, eoe
      integer :: iEdge, k
   
   end subroutine ocn_compute_output_diagnostics!}}}
   
   subroutine mpas_timestep(domain, itimestep, dt, timeStamp)!{{{
   
      use mpas_grid_types
   
      implicit none
   
      type (domain_type), intent(inout) :: domain 
      integer, intent(in) :: itimestep
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(in) :: timeStamp

      type (block_type), pointer :: block_ptr
      integer :: ierr
   
      call ocn_timestep(domain, dt, timeStamp)

      if (config_stats_interval > 0) then
          if (mod(itimestep, config_stats_interval) == 0) then
             call mpas_timer_start("global diagnostics", .false., globalDiagTimer)
             call ocn_compute_global_diagnostics(domain, 2, itimestep, dt);
             call mpas_timer_stop("global diagnostics", globalDiagTimer)
          end if
      end if

      !TODO: replace the above code block with this if we desire to convert config_stats_interval to use alarms
      !if (mpas_is_alarm_ringing(clock, statsAlarmID, ierr=ierr)) then
      !   call mpas_reset_clock_alarm(clock, statsAlarmID, ierr=ierr)

      !   block_ptr => domain % blocklist
      !   if (associated(block_ptr % next)) then
      !      write(0,*) 'Error: computeGlobalDiagnostics assumes ',&
      !                 'that there is only one block per processor.'
      !   end if
   
      !   call mpas_timer_start("global diagnostics")
      !   call ocn_compute_global_diagnostics(domain % dminfo, &
      !            block_ptr % state % time_levs(2) % state, block_ptr % mesh, &
      !            timeStamp, dt)
      !   call mpas_timer_stop("global diagnostics")
      !end if

   end subroutine mpas_timestep!}}}

   subroutine ocn_init_z_level(domain)!{{{
   ! Initialize maxLevel and bouncary grid variables.

      use mpas_grid_types
      use mpas_configure

      implicit none

      type (domain_type), intent(inout) :: domain

      integer :: i, iCell, iEdge, iVertex, k
      type (block_type), pointer :: block

      integer :: iTracer, cell, cell1, cell2
      real (kind=RKIND) :: uhSum, hSum, hEdge1
      real (kind=RKIND), dimension(:), pointer :: referenceBottomDepth, referenceBottomDepthTopOfCell
         
      real (kind=RKIND), dimension(:,:), pointer :: h
      integer :: nVertLevels

      ! Initialize z-level grid variables from h, read in from input file.
      block => domain % blocklist
      do while (associated(block))

         h          => block % state % time_levs(1) % state % h % array
         referenceBottomDepth => block % mesh % referenceBottomDepth % array
         referenceBottomDepthTopOfCell => block % mesh % referenceBottomDepthTopOfCell % array
         nVertLevels = block % mesh % nVertLevels

         ! mrp 120208 right now hZLevel is in the grid.nc file.
         ! We would like to transition to using referenceBottomDepth
         ! as the defining variable instead, and will transition soon.
         ! When the transition is done, hZLevel can be removed from
         ! registry and the following four lines deleted.
         referenceBottomDepth(1) = block % mesh % hZLevel % array(1)
         do k = 2,nVertLevels
            referenceBottomDepth(k) = referenceBottomDepth(k-1) + block % mesh % hZLevel % array(k)
         end do

         ! TopOfCell needed where zero depth for the very top may be referenced.
         referenceBottomDepthTopOfCell(1) = 0.0
         do k = 1,nVertLevels
            referenceBottomDepthTopOfCell(k+1) = referenceBottomDepth(k)
         end do

         ! Compute barotropic velocity at first timestep
         ! This is only done upon start-up.
         if (trim(config_time_integration) == 'unsplit_explicit') then
            block % state % time_levs(1) % state % uBtr % array(:) = 0.0

              block % state % time_levs(1) % state % uBcl % array(:,:) &
            = block % state % time_levs(1) % state % u % array(:,:) 

         elseif (trim(config_time_integration) == 'split_explicit') then

            if (config_filter_btr_mode) then
               do iCell=1,block % mesh % nCells
                  block % state % time_levs(1) % state % h % array(1,iCell) & 
                = block % mesh % referenceBottomDepth % array(1)
               enddo
            endif 

            do iEdge=1,block % mesh % nEdges
               cell1 = block % mesh % cellsOnEdge % array(1,iEdge)
               cell2 = block % mesh % cellsOnEdge % array(2,iEdge)

               ! uBtr = sum(u)/sum(h) on each column
               ! ocn_diagnostic_solve has not yet been called, so compute hEdge 
               ! just for this edge.

               ! hSum is initialized outside the loop because on land boundaries 
               ! maxLevelEdgeTop=0, but I want to initialize hSum with a 
               ! nonzero value to avoid a NaN.
               hEdge1 = 0.5*( &
                   block % state % time_levs(1) % state % h % array(1,cell1) & 
                 + block % state % time_levs(1) % state % h % array(1,cell2) ) 
               uhSum = hEdge1*block % state % time_levs(1) % state % u % array(1,iEdge)
               hSum = hEdge1

               do k=2,block % mesh % maxLevelEdgeTop % array(iEdge)
                  ! ocn_diagnostic_solve has not yet been called, so compute hEdge 
                  ! just for this edge.
                  hEdge1 = 0.5*( &
                      block % state % time_levs(1) % state % h % array(k,cell1) & 
                    + block % state % time_levs(1) % state % h % array(k,cell2) ) 

                  uhSum = uhSum &
                     + hEdge1*block % state % time_levs(1) % state % u % array(k,iEdge)
                  hSum = hSum + hEdge1

               enddo
               block % state % time_levs(1) % state % uBtr % array(iEdge) = uhSum/hsum

               ! uBcl(k,iEdge) = u(k,iEdge) - uBtr(iEdge)
               do k=1,block % mesh % maxLevelEdgeTop % array(iEdge)
                 block % state % time_levs(1) % state % uBcl % array(k,iEdge) &
                 = block % state % time_levs(1) % state % u % array(k,iEdge) &
                 - block % state % time_levs(1) % state % uBtr % array(iEdge)
               enddo

               ! uBcl=0, u=0 on land cells
               do k=block % mesh % maxLevelEdgeTop % array(iEdge)+1, block % mesh % nVertLevels
                 block % state % time_levs(1) % state % uBcl % array(k,iEdge) = 0.0
                 block % state % time_levs(1) % state % u % array(k,iEdge) = 0.0
               enddo
            enddo

            if (config_filter_btr_mode) then
               ! filter uBtr out of initial condition
                block % state % time_levs(1) % state % u % array(:,:) &
              = block % state % time_levs(1) % state % uBcl % array(:,:)

               block % state % time_levs(1) % state % uBtr % array(:) = 0.0
            endif 

         endif

      block => block % next
      end do

   end subroutine ocn_init_z_level!}}}

subroutine ocn_compute_max_level(domain)!{{{
! Initialize maxLevel and bouncary grid variables.

   use mpas_grid_types
   use mpas_configure
   use mpas_constants

   implicit none

   type (domain_type), intent(inout) :: domain

   integer :: i, iCell, iEdge, iVertex, k
   type (block_type), pointer :: block

   real (kind=RKIND), dimension(:,:), pointer :: h, u, u_src, rho
   real (kind=RKIND), dimension(:,:,:), pointer :: tracers
   real (kind=RKIND) :: delta_rho, pi, latCenter, lonCenter, dist
   real (kind=RKIND) :: centerx, centery
   integer :: nCells, nEdges, nVertices, nVertLevels, vertexDegree

   integer, dimension(:), pointer :: &
      maxLevelCell, maxLevelEdgeTop, maxLevelEdgeBot, &
      maxLevelVertexTop, maxLevelVertexBot
   integer, dimension(:,:), pointer :: &
      cellsOnEdge, cellsOnVertex, boundaryEdge, boundaryCell, &
      boundaryVertex, verticesOnEdge, edgeMask, cellMask, vertexMask

   ! Initialize z-level grid variables from h, read in from input file.
   block => domain % blocklist
   do while (associated(block))

      maxLevelCell => block % mesh % maxLevelCell % array
      maxLevelEdgeTop => block % mesh % maxLevelEdgeTop % array
      maxLevelEdgeBot => block % mesh % maxLevelEdgeBot % array
      maxLevelVertexTop => block % mesh % maxLevelVertexTop % array
      maxLevelVertexBot => block % mesh % maxLevelVertexBot % array
      cellsOnEdge    => block % mesh % cellsOnEdge % array
      cellsOnVertex  => block % mesh % cellsOnVertex % array
      verticesOnEdge => block % mesh % verticesOnEdge % array
      boundaryEdge   => block % mesh % boundaryEdge % array
      boundaryCell   => block % mesh % boundaryCell % array
      boundaryVertex => block % mesh % boundaryVertex % array
      edgeMask => block % mesh % edgeMask % array
      cellMask => block % mesh % cellMask % array
      vertexMask => block % mesh % vertexMask % array

      nCells      = block % mesh % nCells
      nEdges      = block % mesh % nEdges
      nVertices   = block % mesh % nVertices
      nVertLevels = block % mesh % nVertLevels
      vertexDegree = block % mesh % vertexDegree

      ! for z-grids, maxLevelCell should be in input state
      ! Isopycnal grid uses all vertical cells
      if (config_vert_grid_type.eq.'isopycnal') then
         maxLevelCell(1:nCells) = nVertLevels
      endif
      maxLevelCell(nCells+1) = 0

      ! maxLevelEdgeTop is the minimum (shallowest) of the surrounding cells
      do iEdge=1,nEdges
         maxLevelEdgeTop(iEdge) = &
            min( maxLevelCell(cellsOnEdge(1,iEdge)), &
                 maxLevelCell(cellsOnEdge(2,iEdge)) )
      end do 
      maxLevelEdgeTop(nEdges+1) = 0

      ! maxLevelEdgeBot is the maximum (deepest) of the surrounding cells
      do iEdge=1,nEdges
         maxLevelEdgeBot(iEdge) = &
            max( maxLevelCell(cellsOnEdge(1,iEdge)), &
                 maxLevelCell(cellsOnEdge(2,iEdge)) )
      end do 
      maxLevelEdgeBot(nEdges+1) = 0

      ! maxLevelVertexBot is the maximum (deepest) of the surrounding cells
      do iVertex = 1,nVertices
         maxLevelVertexBot(iVertex) = maxLevelCell(cellsOnVertex(1,iVertex))
         do i=2,vertexDegree
            maxLevelVertexBot(iVertex) = &
               max( maxLevelVertexBot(iVertex), &
                    maxLevelCell(cellsOnVertex(i,iVertex)))
         end do
      end do 
      maxLevelVertexBot(nVertices+1) = 0

      ! maxLevelVertexTop is the minimum (shallowest) of the surrounding cells
      do iVertex = 1,nVertices
         maxLevelVertexTop(iVertex) = maxLevelCell(cellsOnVertex(1,iVertex))
         do i=2,vertexDegree
            maxLevelVertexTop(iVertex) = &
               min( maxLevelVertexTop(iVertex), &
                    maxLevelCell(cellsOnVertex(i,iVertex)))
         end do
      end do 
      maxLevelVertexTop(nVertices+1) = 0

      ! set boundary edge
      boundaryEdge(:,1:nEdges+1)=1
      edgeMask(:,1:nEdges+1)=0
      do iEdge=1,nEdges
         boundaryEdge(1:maxLevelEdgeTop(iEdge),iEdge)=0
         edgeMask(1:maxLevelEdgeTop(iEdge),iEdge)=1
      end do 

      !
      ! Find cells and vertices that have an edge on the boundary
      !
      boundaryCell(:,1:nCells+1) = 0
      cellMask(:,1:nCells+1) = 1
      boundaryVertex(:,1:nVertices+1) = 0
      vertexMask(:,1:nVertices+1) = 1
      do iEdge=1,nEdges
         do k=1,nVertLevels
            if (boundaryEdge(k,iEdge).eq.1) then
               boundaryCell(k,cellsOnEdge(1,iEdge)) = 1
               boundaryCell(k,cellsOnEdge(2,iEdge)) = 1
               if(maxLevelCell(cellsOnEdge(1,iEdge)) > k) then
                   cellMask(k, cellsOnEdge(1,iEdge)) = 0
               end if
               if(maxLevelCell(cellsOnEdge(2,iEdge)) > k) then
                   cellMask(k, cellsOnEdge(2,iEdge)) = 0
               end if
               boundaryVertex(k,verticesOnEdge(1,iEdge)) = 1
               boundaryVertex(k,verticesOnEdge(2,iEdge)) = 1
               vertexMask(k,verticesOnEdge(1,iEdge)) = 0
               vertexMask(k,verticesOnEdge(2,iEdge)) = 0
            endif
         end do
      end do

      block => block % next
   end do

   ! Note: We do not update halos on maxLevel* variables.  I want the
   ! outside edge of a halo to be zero on each processor.

end subroutine ocn_compute_max_level!}}}
   
   subroutine mpas_core_finalize(domain)!{{{
   
      use mpas_grid_types
   
      implicit none

      type (domain_type), intent(inout) :: domain 
      integer :: ierr

      call mpas_destroy_clock(clock, ierr)

   end subroutine mpas_core_finalize!}}}

   subroutine ocn_compute_mesh_scaling(mesh)!{{{

      use mpas_grid_types
      use mpas_configure

      implicit none

      type (mesh_type), intent(inout) :: mesh

      integer :: iEdge, cell1, cell2
      real (kind=RKIND), dimension(:), pointer :: meshDensity, meshScalingDel2, meshScalingDel4

      meshDensity => mesh % meshDensity % array
      meshScalingDel2 => mesh % meshScalingDel2 % array
      meshScalingDel4 => mesh % meshScalingDel4 % array

      !
      ! Compute the scaling factors to be used in the del2 and del4 dissipation
      !
      meshScalingDel2(:) = 1.0
      meshScalingDel4(:) = 1.0
      if (config_h_ScaleWithMesh) then
         do iEdge=1,mesh%nEdges
            cell1 = mesh % cellsOnEdge % array(1,iEdge)
            cell2 = mesh % cellsOnEdge % array(2,iEdge)
            meshScalingDel2(iEdge) = 1.0 / ( (meshDensity(cell1) + meshDensity(cell2) )/2.0)**(5.0/12.0)
            meshScalingDel4(iEdge) = 1.0 / ( (meshDensity(cell1) + meshDensity(cell2) )/2.0)**(5.0/6.0)
         end do
      end if

   end subroutine ocn_compute_mesh_scaling!}}}

end module mpas_core

! vim: foldmethod=marker
