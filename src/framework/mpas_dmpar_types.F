module mpas_dmpar_types

   use mpas_kind_types

   type dm_info
      integer :: nprocs, my_proc_id, comm, info
      logical :: using_external_comm
   end type dm_info


   type exchange_list
      integer :: procID
      integer :: blockID
      integer :: nlist
      integer, dimension(:), pointer :: list
      type (exchange_list), pointer :: next
      real (kind=RKIND), dimension(:), pointer :: rbuffer
      integer, dimension(:), pointer           :: ibuffer
      integer :: reqID

   end type exchange_list

end module mpas_dmpar_types
