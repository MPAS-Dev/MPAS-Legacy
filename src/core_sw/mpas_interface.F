subroutine mpas_setup_test_case(domain)

   use grid_types
   use test_cases

   implicit none

   type (domain_type), intent(inout) :: domain

   call setup_sw_test_case(domain)

end subroutine mpas_setup_test_case


subroutine mpas_init(block, mesh, dt)

   use grid_types
   use time_integration
   use RBF_interpolation
   use vector_reconstruction

   implicit none

   type (block_type), intent(inout) :: block
   type (grid_meta), intent(inout) :: mesh
   real (kind=RKIND), intent(in) :: dt

   call compute_solve_diagnostics(dt, block % time_levs(1) % state, mesh)

   call rbfInterp_initialize(mesh)
   call init_reconstruct(mesh)
   call reconstruct(block % time_levs(1) % state, mesh)

end subroutine mpas_init


subroutine mpas_query(key, ivalue)

   implicit none

   character (len=256), intent(in) :: key
   integer, intent(out) :: ivalue

   if (index(key,'STORAGE_FACTOR') /= 0) then
      ivalue = 2
   end if

end subroutine mpas_query


subroutine mpas_timestep(domain, itimestep, dt)

   use grid_types
   use time_integration
   use timer
   use global_diagnostics

   implicit none

   type (domain_type), intent(inout) :: domain 
   integer, intent(in) :: itimestep
   real (kind=RKIND), intent(in) :: dt
   type (block_type), pointer :: block_ptr

   call timestep(domain, dt)

   if (mod(itimestep, config_stats_interval) == 0) then
       block_ptr => domain % blocklist
       if(associated(block_ptr % next)) then
           write(0,*) 'Error: computeGlobalDiagnostics assumes ',&
                      'that there is only one block per processor.'
       end if

       call timer_start("global_diagnostics")
       call computeGlobalDiagnostics(domain % dminfo, &
                block_ptr % time_levs(2) % state, block_ptr % mesh, &
                itimestep, dt)
       call timer_stop("global_diagnostics")
   end if

end subroutine mpas_timestep


subroutine mpas_finalize()

   implicit none

end subroutine mpas_finalize
