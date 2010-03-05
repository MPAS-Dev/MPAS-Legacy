module grid_types

   use dmpar

   integer, parameter :: nTimeLevs = 2

#include "super_array_indices.inc"
  
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
   type grid_meta

#include "field_dimensions.inc"

#include "time_invariant_fields.inc"

   end type grid_meta


   ! Derived type for storing model state
   type grid_state

#include "time_varying_fields.inc"

   end type grid_state

   type grid_state_ptr
      type (grid_state), pointer :: state
   end type grid_state_ptr


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
      integer :: storageFactor    ! Additional storage used by time integration scheme

      type (grid_meta), pointer :: mesh
      type (grid_state_ptr), pointer, dimension(:) :: time_levs

      type (grid_state), allocatable, dimension(:) :: intermediate_step

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

   !
   ! Solver interface routines provided by specific dycore
   !
   interface
      subroutine mpas_query(key, ivalue)
         character (len=256), intent(in) :: key
         integer, intent(out) :: ivalue
      end subroutine mpas_query
   end interface


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
      character (len=256) :: key

      nullify(b % prev)
      nullify(b % next)

      allocate(b % time_levs(nTimeLevs))

      allocate(b % mesh)
      call allocate_grid_meta(b % mesh, &
#include "dim_dummy_args.inc"
                             )
  
      do i=1,nTimeLevs
         allocate(b % time_levs(i) % state)
         call allocate_grid_state(b % time_levs(i) % state, b)
      end do

      key = 'STORAGE_FACTOR'
      call mpas_query(key, b % storageFactor)

      ! Allocate storage for intermediate steps used by time integration scheme
      allocate(b % intermediate_step(b % storageFactor))
      do i=1,b % storageFactor
         call allocate_grid_state(b % intermediate_step(i), b)
      end do

      allocate(b % parinfo)

      b % domain => dom

   end subroutine allocate_block


   subroutine allocate_grid_meta(g, &
#include "dim_dummy_args.inc"
                                )

      implicit none

      type (grid_meta), intent(inout) :: g
#include "dim_dummy_decls.inc"

#include "grid_meta_allocs.inc"

   end subroutine allocate_grid_meta


   subroutine allocate_grid_state(s, b)

      implicit none

      type (grid_state), intent(inout) :: s
      type (block_type), pointer :: b

#include "grid_state_allocs.inc"

   end subroutine allocate_grid_state


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

      call deallocate_grid_meta(b % mesh)
      deallocate(b % mesh)
      do i=1,nTimeLevs
         call deallocate_grid_state(b % time_levs(i) % state)
         deallocate(b % time_levs(i) % state)
      end do
      deallocate(b % time_levs)
      do i=1,b % storageFactor
         call deallocate_grid_state(b % intermediate_step(i))
      end do
      deallocate(b % intermediate_step)
      deallocate(b % parinfo)

   end subroutine deallocate_block


   subroutine deallocate_grid_meta(g)

      implicit none

      type (grid_meta), intent(inout) :: g

#include "grid_meta_deallocs.inc"

   end subroutine deallocate_grid_meta


   subroutine deallocate_grid_state(s)

      implicit none

      type (grid_state), intent(inout) :: s

#include "grid_state_deallocs.inc"

   end subroutine deallocate_grid_state


   subroutine copy_state(src, dest)

      implicit none

      type (grid_state), intent(in) :: src
      type (grid_state), intent(inout) :: dest

#include "copy_state.inc"
  
   end subroutine copy_state


   subroutine shift_time_levels(domain)

      implicit none

      type (domain_type), intent(inout) :: domain

      integer :: i
      type (block_type), pointer :: block_ptr
      type (grid_state), pointer :: sptr

      block_ptr => domain % blocklist
      do while (associated(block_ptr))

         sptr => block_ptr % time_levs(1) % state
         do i=1,nTimeLevs-1
            block_ptr % time_levs(i) % state => block_ptr % time_levs(i+1) % state
         end do
         block_ptr % time_levs(nTimeLevs) % state => sptr

         block_ptr => block_ptr % next
      end do

   end subroutine shift_time_levels

end module grid_types
