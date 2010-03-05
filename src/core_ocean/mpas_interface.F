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

   implicit none

   type (block_type), intent(inout) :: block
   type (grid_meta), intent(inout) :: mesh
   real (kind=RKIND), intent(in) :: dt

   call compute_solve_diagnostics(dt, block % time_levs(1) % state, mesh)

   ! xsad 10-02-09:
   ! commenting this out until we incorporate the necessary lapack routines into mpas
   ! call init_reconstruct(block_ptr % mesh)
   ! xsad 10-02-09 end

!   call timer_start("global diagnostics")
!   call computeGlobalDiagnostics(domain % dminfo, block_ptr % time_levs(1) % state, block_ptr % mesh, 0, dt)
!   call timer_stop("global diagnostics")
!   ! xsad 10-02-08 end
!   call output_state_init(output_obj, domain, "OUTPUT")
!   call write_output_frame(output_obj, domain)

end subroutine mpas_init


subroutine mpas_query(key, ivalue)

   implicit none

   character (len=256), intent(in) :: key
   integer, intent(out) :: ivalue

   if (index(key,'STORAGE_FACTOR') /= 0) then
      ivalue = 1
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

!   ! mrp 100120:
!   if (mod(itimestep, config_stats_interval) == 0) then
!      ! xsad 10-02-08:
!      !call write_stats(domain, itimestep, dt)
!      block_ptr => domain % blocklist
!      if(associated(block_ptr % next)) then
!         write(0,*) 'Error: computeGlobalDiagnostics assumes that there is only one block per processor.'
!         call dmpar_abort(domain % dminfo)
!      end if
!
!      call timer_start("global diagnostics")
!      call computeGlobalDiagnostics(domain % dminfo, block_ptr % time_levs(2) % state, block_ptr % mesh, &
!         itimestep, dt)
!      call timer_stop("global diagnostics")
!      ! xsad 10-02-08 end
!   end if
!   ! mrp 100120 end

end subroutine mpas_timestep


subroutine mpas_finalize()

   implicit none

end subroutine mpas_finalize
