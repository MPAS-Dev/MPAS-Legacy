module mpas_subdriver

   use mpas_framework
   use mpas_core

   type (dm_info), pointer :: dminfo
   type (domain_type), pointer :: domain
   type (io_output_object) :: output_obj
   integer :: output_frame


   contains


   subroutine mpas_init()
 
      implicit none

      real (kind=RKIND) :: dt

      call timer_start("total time")
      call timer_start("initialize")


      !
      ! Initialize infrastructure
      !
      call mpas_framework_init(dminfo, domain)


      call input_state_for_domain(domain)


      !
      ! Initialize core
      !
      call mpas_core_init(domain)

      call timer_stop("initialize")


      !
      ! Set up output streams to be written to by the MPAS core
      !
      output_frame = 1
      call output_state_init(output_obj, domain, "OUTPUT")

   end subroutine mpas_init


   subroutine mpas_run()

      implicit none

      call mpas_core_run(domain, output_obj, output_frame)

   end subroutine mpas_run


   subroutine mpas_finalize()
   
      implicit none

      !
      ! Finalize output streams
      !
      call output_state_finalize(output_obj, domain % dminfo)


      !
      ! Finalize core
      !
      call mpas_core_finalize(domain)

      call timer_stop("total time")
      call timer_write()


      !
      ! Finalize infrastructure
      !
      call mpas_framework_finalize(dminfo, domain)

   end subroutine mpas_finalize

end module mpas_subdriver
