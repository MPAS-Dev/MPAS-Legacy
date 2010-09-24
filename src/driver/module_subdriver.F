module subdriver

   use grid_types
   use io_output
   use configure
   use dmpar
   use timer

   integer :: output_frame
   integer :: restart_frame

   interface
      subroutine mpas_init(block, mesh, dt)
         use grid_types
         type (block_type), intent(inout) :: block
         type (mesh_type), intent(inout) :: mesh
         real (kind=RKIND), intent(in) :: dt
      end subroutine mpas_init

      subroutine mpas_timestep(domain, itimestep, dt)
         use grid_types
         type (domain_type), intent(inout) :: domain
         integer, intent(in) :: itimestep
         real (kind=RKIND), intent(in) :: dt
      end subroutine mpas_timestep
   end interface


   contains


   subroutine solve(domain)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Solve 2d shallow water equations over the time range specified in the 
   !   namelist, writing model state to an output file periodically
   !
   ! Input/Output: domain - grid metadata and model state, to be advanced forward 
   !                       in time
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain
     
      integer :: ntimesteps, itimestep
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block_ptr
      type (io_output_object) :: output_obj
      type (io_output_object) :: restart_obj


      ! Eventually, dt should be domain specific
      dt = config_dt
      ntimesteps = config_ntimesteps


      ! Compute diagnostic fields needed in solve loop, and initialize 
      !    simulation time to 0 for all blocks
      block_ptr => domain % blocklist
      do while (associated(block_ptr))
         call mpas_init(block_ptr, block_ptr % mesh, dt)
         if (.not. config_do_restart) block_ptr % state % time_levs(1) % state % xtime % scalar = 0.0
         block_ptr => block_ptr % next
      end do

      ! Before integrating, write out the initial state
      output_frame = 1
      restart_frame = 1
      call output_state_init(output_obj, domain, "OUTPUT")
      call write_output_frame(output_obj, domain)


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
            call write_output_frame(output_obj, domain)
         end if
         if (mod(itimestep, config_restart_interval) == 0 .and. config_restart_interval > 0) then 
            if (restart_frame == 1) call output_state_init(restart_obj, domain, "RESTART")
            call output_state_for_domain(restart_obj, domain, restart_frame)
            restart_frame = restart_frame + 1
         end if
      end do

      call output_state_finalize(output_obj, domain % dminfo)
      if (restart_frame > 1) call output_state_finalize(restart_obj, domain % dminfo)

   end subroutine solve


   subroutine write_output_frame(output_obj, domain)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain and write model state to output file
   !
   ! Input/Output: domain - contains model state; diagnostic field are computed 
   !                        before returning
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

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

      implicit none

      type (state_type), intent(inout) :: state
      type (mesh_type), intent(in) :: grid

      integer :: i, eoe
      integer :: iEdge, k

   end subroutine compute_output_diagnostics


end module subdriver
