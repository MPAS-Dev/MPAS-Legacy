module mpas_framework

   use dmpar
   use grid_types
   use io_input
   use io_output
   use configure
   use timer
   use mpas_timekeeping


   contains

   
   subroutine mpas_framework_init(dminfo, domain)

      implicit none

      type (dm_info), pointer :: dminfo
      type (domain_type), pointer :: domain

      allocate(dminfo)
      call dmpar_init(dminfo)

      call read_namelist(dminfo)

      call allocate_domain(domain, dminfo)
      
      call mpas_timekeeping_init(config_calendar_type)

   end subroutine mpas_framework_init

   
   subroutine mpas_framework_finalize(dminfo, domain)
  
      implicit none

      type (dm_info), pointer :: dminfo
      type (domain_type), pointer :: domain

      call deallocate_domain(domain)

      call dmpar_finalize(dminfo)

      call mpas_timekeeping_finalize()

   end subroutine mpas_framework_finalize

end module mpas_framework
