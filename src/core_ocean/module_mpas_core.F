module mpas_core

   use mpas_framework
   use mpas_timekeeping
   use dmpar
   use test_cases

   type (io_output_object) :: restart_obj
   integer :: restart_frame

   integer :: current_outfile_frames

   type (MPAS_Clock_type) :: clock

   integer, parameter :: outputAlarmID = 1
   integer, parameter :: restartAlarmID = 2
   integer, parameter :: statsAlarmID = 3

   contains

   subroutine mpas_core_init(domain, startTimeStamp)

      use configure
      use grid_types

      implicit none

      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp

      real (kind=RKIND) :: dt
      type (block_type), pointer :: block
      type (dm_info) :: dminfo

      if (.not. config_do_restart) call setup_sw_test_case(domain)

      if (config_vert_grid_type.eq.'isopycnal') then
         print *, ' Using isopycnal coordinates'
      elseif (config_vert_grid_type.eq.'zlevel') then
         print *, ' Using z-level coordinates'
         call init_ZLevel(domain)
      else 
         print *, ' Incorrect choice of config_vert_grid_type:',&
           config_vert_grid_type
         call dmpar_abort(dminfo)
      endif

      call compute_maxLevel(domain)

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

   ! mrp 100316 In order for this to work, we need to pass domain % dminfo as an 
   ! input arguement into mpas_init.  Ask about that later.  For now, there will be
   ! no initial statistics write.
   
   !   call timer_start("global diagnostics")
   !   call computeGlobalDiagnostics(domain % dminfo, block % state % time_levs(1) % state, mesh, 0, dt)
   !   call timer_stop("global diagnostics")
   !   call output_state_init(output_obj, domain, "OUTPUT")
   !   call write_output_frame(output_obj, domain)

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
          call dmpar_finalize(domain % dminfo)
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
      integer :: i, iEdge, iCell, k
   
   
      call compute_solve_diagnostics(dt, block % state % time_levs(1) % state, mesh)
 
      call rbfInterp_initialize(mesh)
      call init_reconstruct(mesh)
      call reconstruct(block % state % time_levs(1) % state, mesh)

      ! initialize velocities and tracers on land to be -1e34
      ! The reconstructed velocity on land will have values not exactly
      ! -1e34 due to the interpolation of reconstruction.

      do iEdge=1,block % mesh % nEdges
         ! mrp 101115 note: in order to include flux boundary conditions, the following
         ! line will need to change.  Right now, set boundary edges between land and 
         ! water to have zero velocity.
         block % state % time_levs(1) % state % u % array( &
             block % mesh % maxLevelEdgeTop % array(iEdge)+1 &
            :block % mesh % maxLevelEdgeBot % array(iEdge), iEdge) = 0.0

         block % state % time_levs(1) % state % u % array( &
             block % mesh % maxLevelEdgeBot % array(iEdge)+1: &
             block % mesh % nVertLevels,iEdge) = -1e34
      end do
      do iCell=1,block % mesh % nCells
         block % state % time_levs(1) % state % tracers % array( &
            :, block % mesh % maxLevelCell % array(iCell)+1 &
              :block % mesh % nVertLevels,iCell) =  -1e34
      end do
   
      if (config_do_restart) then 
          do i=2,nTimeLevs
             call copy_state(block % state % time_levs(i) % state, &
                             block % state % time_levs(1) % state)
          end do
      endif
   
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
      write(0,*) 'Initial time ', timeStamp

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
   
      integer, intent(inout) :: output_frame
      type (domain_type), intent(inout) :: domain
      type (io_output_object), intent(inout) :: output_obj
   
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

      if (config_stats_interval > 0) then
          if (mod(itimestep, config_stats_interval) == 0) then
              block_ptr => domain % blocklist
              if (associated(block_ptr % next)) then
                  write(0,*) 'Error: computeGlobalDiagnostics assumes ',&
                     'that there is only one block per processor.'
              end if

          call timer_start("global diagnostics")
          call computeGlobalDiagnostics(domain % dminfo, &
             block_ptr % state % time_levs(2) % state, block_ptr % mesh, &
             itimestep, dt)
          call timer_stop("global diagnostics")
          end if
      end if

      !TODO: replace the above code block with this if we desire to convert config_stats_interval to use alarms
      !if (MPAS_isAlarmRinging(clock, statsAlarmID, ierr=ierr)) then
      !   call MPAS_resetClockAlarm(clock, statsAlarmID, ierr=ierr)

      !   block_ptr => domain % blocklist
      !   if (associated(block_ptr % next)) then
      !      write(0,*) 'Error: computeGlobalDiagnostics assumes ',&
      !                 'that there is only one block per processor.'
      !   end if
   
      !   call timer_start("global diagnostics")
      !   call computeGlobalDiagnostics(domain % dminfo, &
      !            block_ptr % state % time_levs(2) % state, block_ptr % mesh, &
      !            timeStamp, dt)
      !   call timer_stop("global diagnostics")
      !end if

   end subroutine mpas_timestep


subroutine init_ZLevel(domain)
! Initialize maxLevel and bouncary grid variables.

   use grid_types
   use configure

   implicit none

   type (domain_type), intent(inout) :: domain

   integer :: i, iCell, iEdge, iVertex, k
   type (block_type), pointer :: block

   integer :: iTracer, cell
   real (kind=RKIND), dimension(:), pointer :: &
      hZLevel, zMidZLevel, zTopZLevel, &
      hMeanTopZLevel, hRatioZLevelK, hRatioZLevelKm1
   real (kind=RKIND), dimension(:,:), pointer :: h
   integer :: nVertLevels

   ! Initialize z-level grid variables from h, read in from input file.
   block => domain % blocklist
   do while (associated(block))

      h          => block % state % time_levs(1) % state % h % array
      hZLevel    => block % mesh % hZLevel % array
      zMidZLevel => block % mesh % zMidZLevel % array
      zTopZLevel => block % mesh % zTopZLevel % array
      nVertLevels = block % mesh % nVertLevels
      hMeanTopZLevel    => block % mesh % hMeanTopZLevel % array
      hRatioZLevelK    => block % mesh % hRatioZLevelK % array
      hRatioZLevelKm1    => block % mesh % hRatioZLevelKm1 % array

      ! These should eventually be in an input file.  For now
      ! I just read them in from h(:,1).
      ! Upon restart, the correct hZLevel should be in restart.nc
      if (.not. config_do_restart) hZLevel = h(:,1)

      ! hZLevel should be in the grid.nc and restart.nc file, 
      ! and h for k=1 must be specified there as well.
 
      zTopZLevel(1) = 0.0
      do k = 1,nVertLevels
         zMidZLevel(k) = zTopZLevel(k)-0.5*hZLevel(k)
         zTopZLevel(k+1) = zTopZLevel(k)-  hZLevel(k)
      end do

      hMeanTopZLevel(1) = 0.0
      hRatioZLevelK(1) = 0.0
      hRatioZLevelKm1(1) = 0.0
      do k = 2,nVertLevels
         hMeanTopZLevel(k) = 0.5*(hZLevel(k-1) + hZLevel(k))
         hRatioZLevelK(k) = 0.5*hZLevel(k)/hMeanTopZLevel(k)
         hRatioZLevelKm1(k) = 0.5*hZLevel(k-1)/hMeanTopZLevel(k)
      end do

      block => block % next

   end do

end subroutine init_ZLevel


subroutine compute_maxLevel(domain)
! Initialize maxLevel and bouncary grid variables.

   use grid_types
   use configure
   use constants

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
      boundaryVertex, verticesOnEdge

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
      boundaryEdge=1
      do iEdge=1,nEdges
         boundaryEdge(1:maxLevelEdgeTop(iEdge),iEdge)=0
      end do 

      !
      ! Find cells and vertices that have an edge on the boundary
      !
      boundaryCell(:,:) = 0
      do iEdge=1,nEdges
         do k=1,nVertLevels
            if (boundaryEdge(k,iEdge).eq.1) then
               boundaryCell(k,cellsOnEdge(1,iEdge)) = 1
               boundaryCell(k,cellsOnEdge(2,iEdge)) = 1
               boundaryVertex(k,verticesOnEdge(1,iEdge)) = 1
               boundaryVertex(k,verticesOnEdge(2,iEdge)) = 1
            endif
         end do
      end do

      block => block % next
   end do

   ! Note: We do not update halos on maxLevel* variables.  I want the
   ! outside edge of a halo to be zero on each processor.

end subroutine compute_maxLevel
   
   
   subroutine mpas_core_finalize(domain)
   
      use grid_types
   
      implicit none

      integer :: ierr

      type (domain_type), intent(inout) :: domain 

      if (restart_frame > 1) call output_state_finalize(restart_obj, domain % dminfo)

      call MPAS_destroyClock(clock, ierr)

   end subroutine mpas_core_finalize

end module mpas_core
