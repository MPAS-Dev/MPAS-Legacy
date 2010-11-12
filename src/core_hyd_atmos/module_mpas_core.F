module mpas_core

   use mpas_framework

   type (io_output_object) :: restart_obj
   integer :: restart_frame

  
   contains


     subroutine mpas_core_init(domain)

      use configure
      use grid_types
      use test_cases

      implicit none

      type (domain_type), intent(inout) :: domain

      real (kind=RKIND) :: dt
      type (block_type), pointer :: block


      if (.not. config_do_restart) call setup_hyd_test_case(domain)

      !
      ! Initialize core
      !
      dt = config_dt
      block => domain % blocklist
      do while (associated(block))
         call mpas_init_block(block, block % mesh, dt)
         block => block % next
      end do

      restart_frame = 1

   end subroutine mpas_core_init
 

   subroutine mpas_init_block(block, mesh, dt)
   
      use grid_types
      use advection
      use time_integration
      use RBF_interpolation
      use vector_reconstruction
   
      implicit none
   
      type (block_type), intent(inout) :: block
      type (mesh_type), intent(inout) :: mesh
      real (kind=RKIND), intent(in) :: dt
   
   
      call compute_solver_constants(block % state % time_levs(1) % state, mesh)
      call compute_state_diagnostics(block % state % time_levs(1) % state, mesh)
      call compute_solve_diagnostics(dt, block % state % time_levs(1) % state, mesh)
      call initialize_advection_rk(mesh)
      call rbfInterp_initialize(mesh)
      call init_reconstruct(mesh)
      call reconstruct(block % state % time_levs(1) % state, mesh)
   
      if (.not. config_do_restart) block % state % time_levs(1) % state % xtime % scalar = 0.0
   
   end subroutine mpas_init_block
   
   
   subroutine mpas_core_run(domain, output_obj, output_frame)
   
      use grid_types
      use io_output
      use timer
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (io_output_object), intent(inout) :: output_obj
      integer, intent(inout) :: output_frame
   
      integer :: ntimesteps, itimestep
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block_ptr
   
      ! Eventually, dt should be domain specific
      dt = config_dt
      ntimesteps = config_ntimesteps
   
      call write_output_frame(output_obj, output_frame, domain)
   
      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      do itimestep = 1,ntimesteps
         write(0,*) 'Doing timestep ', itimestep
         call timer_start("time integration")
         call mpas_timestep(domain, itimestep, dt)
         call timer_stop("time integration")
   
         ! Move time level 2 fields back into time level 1 for next time step
         call shift_time_levels_state(domain % blocklist % state)
   
         if (mod(itimestep, config_output_interval) == 0) then
            call write_output_frame(output_obj, output_frame, domain)
         end if
         if (mod(itimestep, config_restart_interval) == 0 .and. config_restart_interval > 0) then
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
   
   
   subroutine mpas_timestep(domain, itimestep, dt)
   
      use grid_types
      use time_integration
   
      implicit none
   
      type (domain_type), intent(inout) :: domain 
      integer, intent(in) :: itimestep
      real (kind=RKIND), intent(in) :: dt
   
      call timestep(domain, dt)
   
   end subroutine mpas_timestep
   
   
   subroutine mpas_core_finalize(domain)
   
      use grid_types
   
      implicit none
   
      type (domain_type), intent(inout) :: domain 

      if (restart_frame > 1) call output_state_finalize(restart_obj, domain % dminfo)
   
   end subroutine mpas_core_finalize

end module mpas_core
