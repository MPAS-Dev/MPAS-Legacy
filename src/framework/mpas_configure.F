module mpas_configure

   use mpas_dmpar

#include "config_defs.inc"

   contains


   subroutine mpas_read_namelist(dminfo)

      implicit none

      type (dm_info), intent(in) :: dminfo

      integer :: funit, ierr

#include "config_namelist_defs.inc"

      funit = 21

      ! Set default values for namelist options
#include "config_set_defaults.inc"

      if (dminfo % my_proc_id == IO_NODE) then
         write(0,*) 'Reading namelist.input'
         open(funit,file='namelist.input',status='old',form='formatted')

#include "config_namelist_reads.inc"
         close(funit)
         write(0,*) ' '
      end if

#include "config_bcast_namelist.inc"

   end subroutine mpas_read_namelist

end module mpas_configure
