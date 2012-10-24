module mpas_grid_types

   use mpas_kind_types
   use mpas_dmpar_types
   use mpas_attlist

   integer, parameter :: nTimeLevs = 2


   ! Derived type describing info for doing I/O specific to a field
   type io_info
      character (len=StrKIND) :: fieldName
      character (len=StrKIND) :: units
      character (len=StrKIND) :: description
      integer, dimension(4) :: start
      integer, dimension(4) :: count
      logical :: input
      logical :: sfc
      logical :: restart
      logical :: output
   end type io_info


   ! Derived type for storing fields
   type field3DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:,:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(3) :: dimNames
      integer, dimension(3) :: dimSizes
      logical :: hasTimeDimension
      logical :: isSuperArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field3DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field3DReal


   ! Derived type for storing fields
   type field2DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(2) :: dimNames
      integer, dimension(2) :: dimSizes
      logical :: hasTimeDimension
      logical :: isSuperArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field2DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field2DReal


   ! Derived type for storing fields
   type field1DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(1) :: dimNames
      integer, dimension(1) :: dimSizes
      logical :: hasTimeDimension
      logical :: isSuperArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field1DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field1DReal


   ! Derived type for storing fields
   type field0DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND) :: scalar

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      logical :: hasTimeDimension
      logical :: isSuperArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DReal


   ! Derived type for storing fields
   type field3DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer, dimension(:,:,:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(3) :: dimNames
      integer, dimension(3) :: dimSizes
      logical :: hasTimeDimension
      logical :: isSuperArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field3DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field3DInteger


   ! Derived type for storing fields
   type field2DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer, dimension(:,:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(2) :: dimNames
      integer, dimension(2) :: dimSizes
      logical :: hasTimeDimension
      logical :: isSuperArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field2DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field2DInteger


   ! Derived type for storing fields
   type field1DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer, dimension(:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(1) :: dimNames
      integer, dimension(1) :: dimSizes
      logical :: hasTimeDimension
      logical :: isSuperArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field1DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field1DInteger


   ! Derived type for storing fields
   type field0DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer :: scalar

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      logical :: hasTimeDimension
      logical :: isSuperArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DInteger


   ! Derived type for storing fields
   type field1DChar
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      character (len=StrKIND), dimension(:), pointer :: array

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(1) :: dimNames
      integer, dimension(1) :: dimSizes
      logical :: hasTimeDimension
      logical :: isSuperArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field1DChar), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field1DChar


   ! Derived type for storing fields
   type field0DChar
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      character (len=StrKIND) :: scalar

      ! Information used by the I/O layer
      type (io_info), pointer :: ioinfo       ! to be removed later
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      logical :: hasTimeDimension
      logical :: isSuperArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DChar), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DChar


   ! Derived type for storing grid meta-data
   type mesh_type

      type (block_type), pointer :: block

#include "field_dimensions.inc"

      logical :: on_a_sphere
      real (kind=RKIND) :: sphere_radius

#include "time_invariant_fields.inc"

   end type mesh_type


#include "variable_groups.inc"


   ! Type for storing (possibly architecture specific) information concerning to parallelism
   type parallel_info
      type (mpas_multihalo_exchange_list), pointer :: cellsToSend            ! List of types describing which cells to send to other blocks
      type (mpas_multihalo_exchange_list), pointer :: cellsToRecv            ! List of types describing which cells to receive from other blocks
      type (mpas_multihalo_exchange_list), pointer :: cellsToCopy            ! List of types describing which cells to copy from other blocks

      type (mpas_multihalo_exchange_list), pointer :: edgesToSend            ! List of types describing which edges to send to other blocks
      type (mpas_multihalo_exchange_list), pointer :: edgesToRecv            ! List of types describing which edges to receive from other blocks
      type (mpas_multihalo_exchange_list), pointer :: edgesToCopy            ! List of types describing which edges to copy from other blocks

      type (mpas_multihalo_exchange_list), pointer :: verticesToSend         ! List of types describing which vertices to send to other blocks
      type (mpas_multihalo_exchange_list), pointer :: verticesToRecv         ! List of types describing which vertices to receive from other blocks
      type (mpas_multihalo_exchange_list), pointer :: verticesToCopy         ! List of types describing which vertices to copy from other blocks
   end type parallel_info


   ! Derived type for storing part of a domain; used as a basic unit of work for a process
   type block_type

#include "block_group_members.inc"

      integer :: blockID   ! Unique global ID number for this block
      integer :: localBlockID  ! Unique local ID number for this block

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

   interface mpas_allocate_scratch_field
     module procedure mpas_allocate_scratch_field1d_integer
     module procedure mpas_allocate_scratch_field2d_integer
     module procedure mpas_allocate_scratch_field3d_integer
     module procedure mpas_allocate_scratch_field1d_real
     module procedure mpas_allocate_scratch_field2d_real
     module procedure mpas_allocate_scratch_field3d_real
     module procedure mpas_allocate_scratch_field1d_char
   end interface

   interface mpas_deallocate_scratch_field
     module procedure mpas_deallocate_scratch_field1d_integer
     module procedure mpas_deallocate_scratch_field2d_integer
     module procedure mpas_deallocate_scratch_field3d_integer
     module procedure mpas_deallocate_scratch_field1d_real
     module procedure mpas_deallocate_scratch_field2d_real
     module procedure mpas_deallocate_scratch_field3d_real
     module procedure mpas_deallocate_scratch_field1d_char
   end interface

   interface mpas_deallocate_field
     module procedure mpas_deallocate_field0d_integer
     module procedure mpas_deallocate_field1d_integer
     module procedure mpas_deallocate_field2d_integer
     module procedure mpas_deallocate_field3d_integer
     module procedure mpas_deallocate_field0d_real
     module procedure mpas_deallocate_field1d_real
     module procedure mpas_deallocate_field2d_real
     module procedure mpas_deallocate_field3d_real
     module procedure mpas_deallocate_field0d_char
     module procedure mpas_deallocate_field1d_char
   end interface


   contains


   subroutine mpas_allocate_domain(dom, dminfo)

      implicit none

      type (domain_type), pointer :: dom
      type (dm_info), pointer :: dminfo

      allocate(dom)
      nullify(dom % blocklist)
      dom % dminfo => dminfo

   end subroutine mpas_allocate_domain


   subroutine mpas_allocate_block(nHaloLayers, b, dom, blockID, &
#include "dim_dummy_args.inc"
                            )

      implicit none

      integer, intent(in) :: nHaloLayers
      type (block_type), pointer :: b
      type (domain_type), pointer :: dom
      integer, intent(in) :: blockID
#include "dim_dummy_decls.inc"



      integer :: i

      b % blockID = blockID

      allocate(b % parinfo)

      b % domain => dom

#include "block_allocs.inc"

   end subroutine mpas_allocate_block


#include "group_alloc_routines.inc"

#include "provis_alloc_routines.inc"


   subroutine mpas_deallocate_domain(dom)!{{{

      implicit none

      type (domain_type), pointer :: dom

      type (block_type), pointer :: block_ptr

      block_ptr => dom % blocklist
      do while (associated(block_ptr))
         call mpas_deallocate_block(block_ptr)
         block_ptr => block_ptr % next
      end do

      deallocate(dom) 

   end subroutine mpas_deallocate_domain!}}}

   subroutine mpas_allocate_scratch_field1d_integer(f, single_block_in)!{{{
       type (field1dInteger), pointer :: f
       logical, intent(in), optional :: single_block_in
       logical :: single_block
       type (field1dInteger), pointer :: f_cursor

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field1d_integer!}}}

   subroutine mpas_allocate_scratch_field2d_integer(f, single_block_in)!{{{
       type (field2dInteger), pointer :: f
       logical, intent(in), optional :: single_block_in
       logical :: single_block
       type (field2dInteger), pointer :: f_cursor

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1), f % dimSizes(2)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field2d_integer!}}}

   subroutine mpas_allocate_scratch_field3d_integer(f, single_block_in)!{{{
       type (field3dInteger), pointer :: f
       logical, intent(in), optional :: single_block_in
       logical :: single_block
       type (field3dInteger), pointer :: f_cursor

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2), f_cursor % dimSizes(3)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1), f % dimSizes(2), f % dimSizes(3)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field3d_integer!}}}

   subroutine mpas_allocate_scratch_field1d_real(f, single_block_in)!{{{
       type (field1dReal), pointer :: f
       logical, intent(in), optional :: single_block_in
       logical :: single_block
       type (field1dReal), pointer :: f_cursor

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field1d_real!}}}

   subroutine mpas_allocate_scratch_field2d_real(f, single_block_in)!{{{
       type (field2dReal), pointer :: f
       logical, intent(in), optional :: single_block_in
       logical :: single_block
       type (field2dReal), pointer :: f_cursor

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1), f % dimSizes(2)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field2d_real!}}}

   subroutine mpas_allocate_scratch_field3d_real(f, single_block_in)!{{{
       type (field3dReal), pointer :: f
       logical, intent(in), optional :: single_block_in
       logical :: single_block
       type (field3dReal), pointer :: f_cursor

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1), f_cursor % dimSizes(2), f_cursor % dimSizes(3)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1), f % dimSizes(2), f % dimSizes(3)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field3d_real!}}}

   subroutine mpas_allocate_scratch_field1d_char(f, single_block_in)!{{{
       type (field1dChar), pointer :: f
       logical, intent(in), optional :: single_block_in
       logical :: single_block
       type (field1dChar), pointer :: f_cursor

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not. single_block) then
          f_cursor => f
          do while(associated(f_cursor))
             if(.not.associated(f_cursor % array)) then
                allocate(f_cursor % array(f_cursor % dimSizes(1)))
             end if
             f_cursor => f_cursor % next
          end do
       else
          if(.not.associated(f % array)) then
            allocate(f % array(f % dimSizes(1)))
          end if
       end if

   end subroutine mpas_allocate_scratch_field1d_char!}}}

   subroutine mpas_deallocate_scratch_field1d_integer(f, single_block_in)!{{{
       type (field1dInteger), pointer :: f
       logical, intent(in), optional :: single_block_in
       logical :: single_block
       type (field1dInteger), pointer :: f_cursor

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field1d_integer!}}}

   subroutine mpas_deallocate_scratch_field2d_integer(f, single_block_in)!{{{
       type (field2dInteger), pointer :: f
       logical, intent(in), optional :: single_block_in
       logical :: single_block
       type (field2dInteger), pointer :: f_cursor

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field2d_integer!}}}

   subroutine mpas_deallocate_scratch_field3d_integer(f, single_block_in)!{{{
       type (field3dInteger), pointer :: f
       logical, intent(in), optional :: single_block_in
       logical :: single_block
       type (field3dInteger), pointer :: f_cursor

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field3d_integer!}}}

   subroutine mpas_deallocate_scratch_field1d_real(f, single_block_in)!{{{
       type (field1dReal), pointer :: f
       logical, intent(in), optional :: single_block_in
       logical :: single_block
       type (field1dReal), pointer :: f_cursor

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field1d_real!}}}

   subroutine mpas_deallocate_scratch_field2d_real(f, single_block_in)!{{{
       type (field2dReal), pointer :: f
       logical, intent(in), optional :: single_block_in
       logical :: single_block
       type (field2dReal), pointer :: f_cursor

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field2d_real!}}}

   subroutine mpas_deallocate_scratch_field3d_real(f, single_block_in)!{{{
       type (field3dReal), pointer :: f
       logical, intent(in), optional :: single_block_in
       logical :: single_block
       type (field3dReal), pointer :: f_cursor

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field3d_real!}}}

   subroutine mpas_deallocate_scratch_field1d_char(f, single_block_in)!{{{
       type (field1dChar), pointer :: f
       logical, intent(in), optional :: single_block_in
       logical :: single_block
       type (field1dChar), pointer :: f_cursor

       if(present(single_block_in)) then
          single_block = single_block_in
       else
          single_block = .false.
       end if

       if(.not.single_block) then
          f_cursor => f
          do while(associated(f_cursor))
            if(associated(f_cursor % array)) then
              deallocate(f_cursor % array)
            end if
   
            f_cursor => f_cursor % next
          end do
       else
          if(associated(f % array)) then
             deallocate(f % array)
          end if
       end if

   end subroutine mpas_deallocate_scratch_field1d_char!}}}


   subroutine mpas_deallocate_field0d_integer(f)!{{{
       type (field0dInteger), pointer :: f
       type (field0dInteger), pointer :: f_cursor

       f_cursor => f

       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         deallocate(f_cursor)
         f_cursor => f
       end do

   end subroutine mpas_deallocate_field0d_integer!}}}

   subroutine mpas_deallocate_field1d_integer(f)!{{{
       type (field1dInteger), pointer :: f
       type (field1dInteger), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field1d_integer!}}}

   subroutine mpas_deallocate_field2d_integer(f)!{{{
       type (field2dInteger), pointer :: f
       type (field2dInteger), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field2d_integer!}}}

   subroutine mpas_deallocate_field3d_integer(f)!{{{
       type (field3dInteger), pointer :: f
       type (field3dInteger), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field3d_integer!}}}

   subroutine mpas_deallocate_field0d_real(f)!{{{
       type (field0dReal), pointer :: f
       type (field0dReal), pointer :: f_cursor

       f_cursor => f

       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field0d_real!}}}

   subroutine mpas_deallocate_field1d_real(f)!{{{
       type (field1dReal), pointer :: f
       type (field1dReal), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field1d_real!}}}

   subroutine mpas_deallocate_field2d_real(f)!{{{
       type (field2dReal), pointer :: f
       type (field2dReal), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field2d_real!}}}

   subroutine mpas_deallocate_field3d_real(f)!{{{
       type (field3dReal), pointer :: f
       type (field3dReal), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field3d_real!}}}

   subroutine mpas_deallocate_field0d_char(f)!{{{
       type (field0dChar), pointer :: f
       type (field0dChar), pointer :: f_cursor

       f_cursor => f

       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         deallocate(f_cursor)
         f_cursor => f
       end do

   end subroutine mpas_deallocate_field0d_char!}}}

   subroutine mpas_deallocate_field1d_char(f)!{{{
       type (field1dChar), pointer :: f
       type (field1dChar), pointer :: f_cursor

       f_cursor => f
       do while(associated(f_cursor))
         if(associated(f % next)) then
           f => f % next
         else
           nullify(f)
         end if

         if(associated(f_cursor % ioinfo)) then
           deallocate(f_cursor % ioinfo)
         end if

         if(associated(f_cursor % array)) then
           deallocate(f_cursor % array)
         end if

         deallocate(f_cursor)

         f_cursor => f
       end do

   end subroutine mpas_deallocate_field1d_char!}}}

   subroutine mpas_deallocate_block(b)!{{{
 
      implicit none

      type (block_type), intent(inout) :: b

      integer :: i

      ! BUG: It seems like we should be deallocating the exchange lists before we 
      !      deallocate the array of head pointers and the parinfo type...

      deallocate(b % parinfo % cellsToSend)
      deallocate(b % parinfo % cellsToRecv)
      deallocate(b % parinfo % cellsToCopy)

      deallocate(b % parinfo % edgesToSend)
      deallocate(b % parinfo % edgesToRecv)
      deallocate(b % parinfo % edgesToCopy)

      deallocate(b % parinfo % verticesToSend)
      deallocate(b % parinfo % verticesToRecv)
      deallocate(b % parinfo % verticesToCopy)

      deallocate(b % parinfo)

#include "block_deallocs.inc"

   end subroutine mpas_deallocate_block!}}}


#include "group_dealloc_routines.inc"


#include "group_copy_routines.inc"


#include "group_shift_level_routines.inc"


#include "field_links.inc"


end module mpas_grid_types
