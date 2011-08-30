module mpas_core

   use mpas_framework
   use mpas_timekeeping

   type (io_output_object) :: restart_obj
   integer :: restart_frame
   integer :: current_outfile_frames

   type (MPAS_Clock_type) :: clock

   integer, parameter :: outputAlarmID = 1
   integer, parameter :: restartAlarmID = 2
   !integer, parameter :: statsAlarmID = 3

   contains

   subroutine mpas_core_init(domain, startTimeStamp)
   
      use configure
      use grid_types
      use test_cases
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp
   
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block


      if (.not. config_do_restart) call setup_sw_test_case(domain)

      !
      ! Initialize core
      !
      dt = config_dt

      call simulation_clock_init(domain, dt, startTimeStamp)

      block => domain % blocklist
      do while (associated(block))
         call mpas_init_block(block, block % mesh, dt)
         block % state % time_levs(1) % state % xtime % scalar = startTimeStamp
         block => block % next
      end do

      restart_frame = 1
      current_outfile_frames = 0

   end subroutine mpas_core_init


   subroutine simulation_clock_init(domain, dt, startTimeStamp)

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(out) :: startTimeStamp

      type (MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      integer :: ierr

      call MPAS_setTime(curr_time=startTime, dateTimeString=config_start_time, ierr=ierr)
      call MPAS_setTimeInterval(timeStep, dt=dt, ierr=ierr)

      if (trim(config_run_duration) /= "none") then
         call MPAS_setTimeInterval(runDuration, timeString=config_run_duration, ierr=ierr)
         call MPAS_createClock(clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=ierr)

         if (trim(config_stop_time) /= "none") then
            call MPAS_setTime(curr_time=stopTime, dateTimeString=config_stop_time, ierr=ierr)
            if(startTime + runduration /= stopTime) then
               write(0,*) 'Warning: config_run_duration and config_stop_time are inconsitent: using config_run_duration.'
            end if
         end if
      else if (trim(config_stop_time) /= "none") then
         call MPAS_setTime(curr_time=stopTime, dateTimeString=config_stop_time, ierr=ierr)
         call MPAS_createClock(clock, startTime=startTime, timeStep=timeStep, stopTime=stopTime, ierr=ierr)
      else
          write(0,*) 'Error: Neither config_run_duration nor config_stop_time were specified.'
          call dmpar_abort(domain % dminfo)
      end if

      ! set output alarm
      call MPAS_setTimeInterval(alarmTimeStep, timeString=config_output_interval, ierr=ierr)
      alarmStartTime = startTime + alarmTimeStep
      call MPAS_addClockAlarm(clock, outputAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)

      ! set restart alarm, if necessary
      if (trim(config_restart_interval) /= "none") then
         call MPAS_setTimeInterval(alarmTimeStep, timeString=config_restart_interval, ierr=ierr)
         alarmStartTime = startTime + alarmTimeStep
         call MPAS_addClockAlarm(clock, restartAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
      end if

      !TODO: use this code if we desire to convert config_stats_interval to alarms 
      !(must also change config_stats_interval type to character) 
      ! set stats alarm, if necessary
      !if (trim(config_stats_interval) /= "none") then      
      !   call MPAS_setTimeInterval(alarmTimeStep, timeString=config_stats_interval, ierr=ierr)
      !   alarmStartTime = startTime + alarmTimeStep
      !   call MPAS_addClockAlarm(clock, statsAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
      !end if

      call MPAS_getTime(curr_time=startTime, dateTimeString=startTimeStamp, ierr=ierr)

   end subroutine simulation_clock_init


   subroutine mpas_init_block(block, mesh, dt)
   
      use grid_types
      use time_integration
      use RBF_interpolation
      use vector_reconstruction
   
      implicit none
   
      type (block_type), intent(inout) :: block
      type (mesh_type), intent(inout) :: mesh
      real (kind=RKIND), intent(in) :: dt
   

      call compute_solve_diagnostics(dt, block % state % time_levs(1) % state, mesh)
      call compute_mesh_scaling(mesh) 

      call rbfInterp_initialize(mesh)
      call init_reconstruct(mesh)
      call reconstruct(block % state % time_levs(1) % state, mesh)
   
   end subroutine mpas_init_block
   
   
   subroutine mpas_core_run(domain, output_obj, output_frame)
   
      use grid_types
      use io_output
      use timer
   
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

      currTime = MPAS_getClockTime(clock, MPAS_NOW, ierr)
      call MPAS_getTime(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         
      write(0,*) 'Initial timestep ', timeStamp

      call write_output_frame(output_obj, output_frame, domain)

      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      itimestep = 0
      do while (.not. MPAS_isClockStopTime(clock))

         itimestep = itimestep + 1
         call MPAS_advanceClock(clock)

         currTime = MPAS_getClockTime(clock, MPAS_NOW, ierr)
         call MPAS_getTime(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         
         write(0,*) 'Doing timestep ', timeStamp

         call timer_start("time integration")
         call mpas_timestep(domain, itimestep, dt, timeStamp)
         call timer_stop("time integration")

         ! Move time level 2 fields back into time level 1 for next time step
         call shift_time_levels_state(domain % blocklist % state)

         !TODO: MPAS_getClockRingingAlarms is probably faster than multiple MPAS_isAlarmRinging...

         if (MPAS_isAlarmRinging(clock, outputAlarmID, ierr=ierr)) then
            call MPAS_resetClockAlarm(clock, outputAlarmID, ierr=ierr)
            if(output_frame == 1) call output_state_init(output_obj, domain, "OUTPUT", trim(timeStamp)) ! output_frame will always be > 1 here unless it is reset after the output file is finalized
            call write_output_frame(output_obj, output_frame, domain)
         end if

         if (MPAS_isAlarmRinging(clock, restartAlarmID, ierr=ierr)) then
            call MPAS_resetClockAlarm(clock, restartAlarmID, ierr=ierr)
            if (restart_frame == 1) call output_state_init(restart_obj, domain, "RESTART")
            call output_state_for_domain(restart_obj, domain, restart_frame)
            restart_frame = restart_frame + 1
         end if

      end do

   end subroutine mpas_core_run
   
   
   subroutine write_output_frame(output_obj, output_frame, domain)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain and write model state to output file
   !
   ! Input/Output: domain - contains model state; diagnostic field are computed
   !                        before returning
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use grid_types
      use io_output
   
      implicit none

      type (io_output_object), intent(inout) :: output_obj
      integer, intent(inout) :: output_frame
      type (domain_type), intent(inout) :: domain
   
      integer :: i, j, k
      integer :: eoe
      type (block_type), pointer :: block_ptr
   
      block_ptr => domain % blocklist
      do while (associated(block_ptr))
         call compute_output_diagnostics(block_ptr % state % time_levs(1) % state, block_ptr % mesh)
         block_ptr => block_ptr % next
      end do
   
      call output_state_for_domain(output_obj, domain, output_frame)
      output_frame = output_frame + 1

      ! if the maximum number of frames per outfile has been reached, finalize outfile and reset frame   
      if (config_frames_per_outfile > 0) then
         current_outfile_frames = current_outfile_frames + 1            
         if(current_outfile_frames >= config_frames_per_outfile) then
            current_outfile_frames = 0
            call output_state_finalize(output_obj, domain % dminfo)
            output_frame = 1
         end if
      end if

   end subroutine write_output_frame
   
   
   subroutine compute_output_diagnostics(state, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain
   !
   ! Input: state - contains model prognostic fields
   !        grid  - contains grid metadata
   !
   ! Output: state - upon returning, diagnostic fields will have be computed
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use grid_types
   
      implicit none
   
      type (state_type), intent(inout) :: state
      type (mesh_type), intent(in) :: grid
   
      integer :: i, eoe
      integer :: iEdge, k
   
   end subroutine compute_output_diagnostics
   
   
   subroutine mpas_timestep(domain, itimestep, dt, timeStamp)
   
      use grid_types
      use time_integration
      use timer
      use global_diagnostics
   
      implicit none
   
      type (domain_type), intent(inout) :: domain 
      integer, intent(in) :: itimestep
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(in) :: timeStamp
      
      type (block_type), pointer :: block_ptr
      integer :: ierr
   
      call timestep(domain, dt, timeStamp)
   
      if(config_stats_interval .gt. 0) then
          if(mod(itimestep, config_stats_interval) == 0) then
              block_ptr => domain % blocklist
              if(associated(block_ptr % next)) then
                  write(0,*) 'Error: computeGlobalDiagnostics assumes ',&
                             'that there is only one block per processor.'
              end if
   
              call timer_start("global_diagnostics")
              call computeGlobalDiagnostics(domain % dminfo, &
                       block_ptr % state % time_levs(2) % state, block_ptr % mesh, &
                       itimestep, dt)
              call timer_stop("global_diagnostics")
          end if
      end if

      !TODO: replace the above code block with this if we desire to convert config_stats_interval to use alarms
      !if (MPAS_isAlarmRinging(clock, statsAlarmID, ierr=ierr)) then
      !   call MPAS_resetClockAlarm(clock, statsAlarmID, ierr=ierr)

      !   block_ptr => domain % blocklist
      !   if(associated(block_ptr % next)) then
      !      write(0,*) 'Error: computeGlobalDiagnostics assumes ',&
      !                 'that there is only one block per processor.'
      !   end if

      !   call timer_start("global_diagnostics")
      !   call computeGlobalDiagnostics(domain % dminfo, &
      !            block_ptr % state % time_levs(2) % state, block_ptr % mesh, &
      !            timeStamp, dt)
      !   call timer_stop("global_diagnostics")
      !end if
   
   end subroutine mpas_timestep
   
   
   subroutine mpas_core_finalize(domain)
   
      use grid_types
   
      implicit none

      integer :: ierr

      type (domain_type), intent(inout) :: domain 

      if (restart_frame > 1) call output_state_finalize(restart_obj, domain % dminfo)

      call MPAS_destroyClock(clock, ierr)

   end subroutine mpas_core_finalize


   subroutine compute_mesh_scaling(mesh)

      use grid_types

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

   end subroutine compute_mesh_scaling

end module mpas_core
