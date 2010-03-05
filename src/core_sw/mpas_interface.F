subroutine mpas_setup_test_case(domain)

   use grid_types
   use test_cases
   use vector_reconstruction

   implicit none

   type (domain_type), intent(inout) :: domain

   call setup_sw_test_case(domain)

end subroutine mpas_setup_test_case


subroutine mpas_init(block, mesh, dt)

   use grid_types
   use time_integration

   implicit none

   type (block_type), intent(inout) :: block
   type (grid_meta), intent(inout) :: mesh
   real (kind=RKIND), intent(in) :: dt

   call compute_solve_diagnostics(dt, block % time_levs(1) % state, mesh)
   call init_reconstruct(mesh)

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

   implicit none

   type (domain_type), intent(inout) :: domain 
   integer, intent(in) :: itimestep
   real (kind=RKIND), intent(in) :: dt

   call timestep(domain, dt)

end subroutine mpas_timestep


subroutine mpas_finalize()

   implicit none

end subroutine mpas_finalize
