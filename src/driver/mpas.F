program mpas

   use grid_types
   use configure
   use io_input
   use dmpar
   use timer
   use subdriver

   implicit none

   interface mpas_setup_test_case
      subroutine mpas_setup_test_case(domain)
         use grid_types
         type (domain_type), intent(inout) :: domain
      end subroutine mpas_setup_test_case
   end interface mpas_setup_test_case


   type (dm_info), pointer :: dminfo
   type (domain_type), pointer :: domain

   allocate(dminfo)
   call dmpar_init(dminfo)

   call read_namelist(dminfo)

   call timer_start("total time")

   call timer_start("initialize")
   call allocate_domain(domain, dminfo)

   call input_state_for_domain(domain)

   if (.not. config_do_restart) call mpas_setup_test_case(domain)
   call timer_stop("initialize")

   call solve(domain) 

   call deallocate_domain(domain)

   call timer_stop("total time")
   call timer_write()

   call dmpar_finalize(dminfo)

   stop

end program mpas
