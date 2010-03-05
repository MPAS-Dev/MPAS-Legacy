module configure

   use dmpar

#include "config_defs.inc"

   contains


   subroutine read_namelist(dminfo)

      implicit none

      type (dm_info), intent(in) :: dminfo

      integer :: funit

#include "config_namelist_defs.inc"

      funit = 21

      ! Set default values for namelist options
#include "config_set_defaults.inc"

      if (dminfo % my_proc_id == IO_NODE) then
         open(funit,file='namelist.input',status='old',form='formatted')

#include "config_namelist_reads.inc"
         close(funit)
      end if

#include "config_bcast_namelist.inc"

   end subroutine read_namelist

end module configure
