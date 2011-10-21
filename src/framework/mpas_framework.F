module mpas_framework

   use mpas_dmpar
   use mpas_grid_types
   use mpas_io_input
   use mpas_io_output
   use mpas_configure
   use mpas_timer
   use mpas_timekeeping


   contains

   
   subroutine mpas_framework_init(dminfo, domain)

      implicit none

      type (dm_info), pointer :: dminfo
      type (domain_type), pointer :: domain

      allocate(dminfo)
      call mpas_dmpar_init(dminfo)

      call mpas_read_namelist(dminfo)

      call mpas_allocate_domain(domain, dminfo)
      
      call mpas_timekeeping_init(config_calendar_type)

   end subroutine mpas_framework_init

   
   subroutine mpas_framework_finalize(dminfo, domain)
  
      implicit none

      type (dm_info), pointer :: dminfo
      type (domain_type), pointer :: domain

      call mpas_deallocate_domain(domain)

      call mpas_dmpar_finalize(dminfo)

      call mpas_timekeeping_finalize()

   end subroutine mpas_framework_finalize

end module mpas_framework
