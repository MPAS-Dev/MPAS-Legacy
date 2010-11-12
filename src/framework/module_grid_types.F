module grid_types

   use dmpar

   integer, parameter :: nTimeLevs = 2


   ! Derived type describing info for doing I/O specific to a field
   type io_info
      character (len=1024) :: fieldName
      integer, dimension(4) :: start
      integer, dimension(4) :: count
      logical :: input
      logical :: restart
      logical :: output
   end type io_info


   ! Derived type for storing fields
   type field3DReal
      type (block_type), pointer :: block
      real (kind=RKIND), dimension(:,:,:), pointer :: array
      type (io_info), pointer :: ioinfo
   end type field3DReal


   ! Derived type for storing fields
   type field2DReal
      type (block_type), pointer :: block
      real (kind=RKIND), dimension(:,:), pointer :: array
      type (io_info), pointer :: ioinfo
   end type field2DReal


   ! Derived type for storing fields
   type field1DReal
      type (block_type), pointer :: block
      real (kind=RKIND), dimension(:), pointer :: array
      type (io_info), pointer :: ioinfo
   end type field1DReal


   ! Derived type for storing fields
   type field0DReal
      type (block_type), pointer :: block
      real (kind=RKIND) :: scalar
      type (io_info), pointer :: ioinfo
   end type field0DReal


   ! Derived type for storing fields
   type field2DInteger
      type (block_type), pointer :: block
      integer, dimension(:,:), pointer :: array
      type (io_info), pointer :: ioinfo
   end type field2DInteger


   ! Derived type for storing fields
   type field1DInteger
      type (block_type), pointer :: block
      integer, dimension(:), pointer :: array
      type (io_info), pointer :: ioinfo
   end type field1DInteger


   ! Derived type for storing grid meta-data
   type mesh_type

#include "field_dimensions.inc"

      logical :: on_a_sphere
      real (kind=RKIND) :: sphere_radius

#include "time_invariant_fields.inc"

   end type mesh_type


#include "variable_groups.inc"


   ! Type for storing (possibly architecture specific) information concerning to parallelism
   type parallel_info
      type (exchange_list), pointer :: cellsToSend            ! List of types describing which cells to send to other blocks
      type (exchange_list), pointer :: cellsToRecv            ! List of types describing which cells to receive from other blocks
      type (exchange_list), pointer :: edgesToSend            ! List of types describing which edges to send to other blocks
      type (exchange_list), pointer :: edgesToRecv            ! List of types describing which edges to receive from other blocks
      type (exchange_list), pointer :: verticesToSend         ! List of types describing which vertices to send to other blocks
      type (exchange_list), pointer :: verticesToRecv         ! List of types describing which vertices to receive from other blocks
   end type parallel_info


   ! Derived type for storing part of a domain; used as a basic unit of work for a process
   type block_type

#include "block_group_members.inc"

      type (domain_type), pointer :: domain

      type (parallel_info), pointer :: parinfo

      type (block_type), pointer :: prev, next
   end type block_type


   ! Derived type for storing list of blocks from a domain to be handled by a process
   type domain_type
      type (block_type), pointer :: blocklist
   
      ! Also store parallelization info here
      type (dm_info), pointer :: dminfo
   end type domain_type


   contains


   subroutine allocate_domain(dom, dminfo)

      implicit none

      type (domain_type), pointer :: dom
      type (dm_info), pointer :: dminfo

      allocate(dom)
      nullify(dom % blocklist)
      dom % dminfo => dminfo

   end subroutine allocate_domain


   subroutine allocate_block(b, dom, &
#include "dim_dummy_args.inc"
                            )

      implicit none

      type (block_type), pointer :: b
      type (domain_type), pointer :: dom
#include "dim_dummy_decls.inc"

      integer :: i

      nullify(b % prev)
      nullify(b % next)

      allocate(b % parinfo)

      b % domain => dom

#include "block_allocs.inc"

   end subroutine allocate_block


#include "group_alloc_routines.inc"


   subroutine deallocate_domain(dom)

      implicit none

      type (domain_type), pointer :: dom

      type (block_type), pointer :: block_ptr

      block_ptr => dom % blocklist
      do while (associated(block_ptr))
         call deallocate_block(block_ptr)
         block_ptr => block_ptr % next
      end do

      deallocate(dom) 

   end subroutine deallocate_domain


   subroutine deallocate_block(b)
 
      implicit none

      type (block_type), intent(inout) :: b

      integer :: i

      deallocate(b % parinfo)

#include "block_deallocs.inc"

   end subroutine deallocate_block


#include "group_dealloc_routines.inc"


#include "group_copy_routines.inc"


#include "group_shift_level_routines.inc"

end module grid_types
