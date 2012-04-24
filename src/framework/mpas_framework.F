module mpas_framework

   use mpas_dmpar
   use mpas_grid_types
   use mpas_io_input
   use mpas_io_output
   use mpas_configure
   use mpas_timer
   use mpas_timekeeping
   use mpas_io


   contains

   
   subroutine mpas_framework_init(dminfo, domain)

      implicit none

      type (dm_info), pointer :: dminfo
      type (domain_type), pointer :: domain

      integer :: pio_num_iotasks
      integer :: pio_stride

      allocate(dminfo)
      call mpas_dmpar_init(dminfo)

      call mpas_read_namelist(dminfo)

      call mpas_allocate_domain(domain, dminfo)
      
      call mpas_timekeeping_init(config_calendar_type)

      pio_num_iotasks = config_pio_num_iotasks
      pio_stride = config_pio_stride
      if (pio_num_iotasks == 0) then
         pio_num_iotasks = domain % dminfo % nprocs
      end if
      call MPAS_io_init(dminfo, pio_num_iotasks, pio_stride)

   end subroutine mpas_framework_init

   
   subroutine mpas_framework_finalize(dminfo, domain)
  
      implicit none

      type (dm_info), pointer :: dminfo
      type (domain_type), pointer :: domain

      call MPAS_io_finalize()

      call mpas_deallocate_domain(domain)

      call mpas_dmpar_finalize(dminfo)

      call mpas_timekeeping_finalize()

   end subroutine mpas_framework_finalize

end module mpas_framework
