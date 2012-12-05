module mpas_dmpar

   use mpas_dmpar_types
   use mpas_grid_types
   use mpas_sort
   use mpas_hash

#ifdef _MPI
include 'mpif.h'
   integer, parameter :: MPI_INTEGERKIND = MPI_INTEGER

#ifdef SINGLE_PRECISION
   integer, parameter :: MPI_REALKIND = MPI_REAL
#else
   integer, parameter :: MPI_REALKIND = MPI_DOUBLE_PRECISION
#endif
#endif

   integer, parameter :: IO_NODE = 0
   integer, parameter :: BUFSIZE = 6000

   interface mpas_dmpar_alltoall_field
      module procedure mpas_dmpar_alltoall_field1d_integer
      module procedure mpas_dmpar_alltoall_field2d_integer
      module procedure mpas_dmpar_alltoall_field1d_real
      module procedure mpas_dmpar_alltoall_field2d_real
      module procedure mpas_dmpar_alltoall_field3d_real
   end interface

   private :: mpas_dmpar_alltoall_field1d_integer
   private :: mpas_dmpar_alltoall_field2d_integer
   private :: mpas_dmpar_alltoall_field1d_real
   private :: mpas_dmpar_alltoall_field2d_real
   private :: mpas_dmpar_alltoall_field3d_real


   interface mpas_dmpar_exch_halo_field
      module procedure mpas_dmpar_exch_halo_field1d_integer
      module procedure mpas_dmpar_exch_halo_field2d_integer
      module procedure mpas_dmpar_exch_halo_field3d_integer
      module procedure mpas_dmpar_exch_halo_field1d_real
      module procedure mpas_dmpar_exch_halo_field2d_real
      module procedure mpas_dmpar_exch_halo_field3d_real
   end interface

   private :: mpas_dmpar_exch_halo_field1d_integer
   private :: mpas_dmpar_exch_halo_field2d_integer
   private :: mpas_dmpar_exch_halo_field3d_integer
   private :: mpas_dmpar_exch_halo_field1d_real
   private :: mpas_dmpar_exch_halo_field2d_real
   private :: mpas_dmpar_exch_halo_field3d_real

   interface mpas_dmpar_copy_field
      module procedure mpas_dmpar_copy_field1d_integer
      module procedure mpas_dmpar_copy_field2d_integer
      module procedure mpas_dmpar_copy_field3d_integer
      module procedure mpas_dmpar_copy_field1d_real
      module procedure mpas_dmpar_copy_field2d_real
      module procedure mpas_dmpar_copy_field3d_real
   end interface

   private :: mpas_dmpar_copy_field1d_integer
   private :: mpas_dmpar_copy_field2d_integer
   private :: mpas_dmpar_copy_field3d_integer
   private :: mpas_dmpar_copy_field1d_real
   private :: mpas_dmpar_copy_field2d_real
   private :: mpas_dmpar_copy_field3d_real

   contains

   subroutine mpas_dmpar_init(dminfo, mpi_comm)!{{{

      implicit none

      type (dm_info), intent(inout) :: dminfo
      integer, intent(in), optional :: mpi_comm     ! Optional: externally-supplied MPI communicator

#ifdef _MPI
      integer :: mpi_rank, mpi_size
      integer :: mpi_ierr

      if (present(mpi_comm)) then
         dminfo % comm = mpi_comm
         dminfo % using_external_comm = .true.
      else
         call MPI_Init(mpi_ierr)
         dminfo % comm = MPI_COMM_WORLD
         dminfo % using_external_comm = .false.
      end if

      ! Find out our rank and the total number of processors
      call MPI_Comm_rank(dminfo % comm, mpi_rank, mpi_ierr)
      call MPI_Comm_size(dminfo % comm, mpi_size, mpi_ierr)

      dminfo % nprocs = mpi_size
      dminfo % my_proc_id = mpi_rank

      write(0,'(a,i5,a,i5,a)') 'task ', mpi_rank, ' of ', mpi_size, &
        ' is running'

#ifndef MPAS_CESM
      call open_streams(dminfo % my_proc_id)
#endif

      dminfo % info = MPI_INFO_NULL
#else
      dminfo % comm = 0
      dminfo % my_proc_id = IO_NODE
      dminfo % nprocs = 1
      dminfo % using_external_comm = .false.
#endif

   end subroutine mpas_dmpar_init!}}}

   subroutine mpas_dmpar_finalize(dminfo)!{{{

      implicit none

      type (dm_info), intent(inout) :: dminfo

#ifdef _MPI
      integer :: mpi_ierr

      if (.not. dminfo % using_external_comm) then
         call MPI_Finalize(mpi_ierr)
      end if
#endif

   end subroutine mpas_dmpar_finalize!}}}

   subroutine mpas_dmpar_abort(dminfo)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo

#ifdef _MPI
      integer :: mpi_ierr, mpi_errcode

      call MPI_Abort(dminfo % comm, mpi_errcode, mpi_ierr)
#endif

      stop

   end subroutine mpas_dmpar_abort!}}}

   subroutine mpas_dmpar_global_abort(mesg)!{{{

      implicit none

      character (len=*), intent(in) :: mesg

#ifdef _MPI
      integer :: mpi_ierr, mpi_errcode

      write(0,*) trim(mesg)
      call MPI_Abort(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
#endif

      write(0,*) trim(mesg)
      stop

   end subroutine mpas_dmpar_global_abort!}}}

   subroutine mpas_dmpar_bcast_int(dminfo, i)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(inout) :: i

#ifdef _MPI
      integer :: mpi_ierr

      call MPI_Bcast(i, 1, MPI_INTEGERKIND, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_bcast_int!}}}

   subroutine mpas_dmpar_bcast_ints(dminfo, n, iarray)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: n
      integer, dimension(n), intent(inout) :: iarray

#ifdef _MPI
      integer :: mpi_ierr

      call MPI_Bcast(iarray, n, MPI_INTEGERKIND, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_bcast_ints!}}}

   subroutine mpas_dmpar_bcast_real(dminfo, r)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      real (kind=RKIND), intent(inout) :: r

#ifdef _MPI
      integer :: mpi_ierr

      call MPI_Bcast(r, 1, MPI_REALKIND, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_bcast_real!}}}

   subroutine mpas_dmpar_bcast_reals(dminfo, n, rarray)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: n
      real (kind=RKIND), dimension(n), intent(inout) :: rarray

#ifdef _MPI
      integer :: mpi_ierr

      call MPI_Bcast(rarray, n, MPI_REALKIND, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_bcast_reals!}}}

   subroutine mpas_dmpar_bcast_logical(dminfo, l)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      logical, intent(inout) :: l

#ifdef _MPI
      integer :: mpi_ierr
      integer :: itemp

      if (dminfo % my_proc_id == IO_NODE) then
         if (l) then
            itemp = 1
         else
            itemp = 0
         end if
      end if

      call MPI_Bcast(itemp, 1, MPI_INTEGERKIND, IO_NODE, dminfo % comm, mpi_ierr)

      if (itemp == 1) then
         l = .true.
      else
         l = .false.
      end if
#endif

   end subroutine mpas_dmpar_bcast_logical!}}}

   subroutine mpas_dmpar_bcast_char(dminfo, c)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      character (len=*), intent(inout) :: c

#ifdef _MPI
      integer :: mpi_ierr

      call MPI_Bcast(c, len(c), MPI_CHARACTER, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_bcast_char!}}}

   subroutine mpas_dmpar_sum_int(dminfo, i, isum)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: i
      integer, intent(out) :: isum

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(i, isum, 1, MPI_INTEGERKIND, MPI_SUM, dminfo % comm, mpi_ierr)
#else
      isum = i
#endif

   end subroutine mpas_dmpar_sum_int!}}}

   subroutine mpas_dmpar_sum_real(dminfo, r, rsum)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      real(kind=RKIND), intent(in) :: r
      real(kind=RKIND), intent(out) :: rsum

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(r, rsum, 1, MPI_REALKIND, MPI_SUM, dminfo % comm, mpi_ierr)
#else
      rsum = r
#endif

   end subroutine mpas_dmpar_sum_real!}}}

   subroutine mpas_dmpar_min_int(dminfo, i, imin)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: i
      integer, intent(out) :: imin

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(i, imin, 1, MPI_INTEGERKIND, MPI_MIN, dminfo % comm, mpi_ierr)
#else
      imin = i
#endif

   end subroutine mpas_dmpar_min_int!}}}

   subroutine mpas_dmpar_min_real(dminfo, r, rmin)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      real(kind=RKIND), intent(in) :: r
      real(kind=RKIND), intent(out) :: rmin

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(r, rmin, 1, MPI_REALKIND, MPI_MIN, dminfo % comm, mpi_ierr)
#else
      rmin = r
#endif

   end subroutine mpas_dmpar_min_real!}}}

   subroutine mpas_dmpar_max_int(dminfo, i, imax)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: i
      integer, intent(out) :: imax
      
      integer :: mpi_ierr 
      
#ifdef _MPI
      call MPI_Allreduce(i, imax, 1, MPI_INTEGERKIND, MPI_MAX, dminfo % comm, mpi_ierr)
#else
      imax = i
#endif

   end subroutine mpas_dmpar_max_int!}}}

   subroutine mpas_dmpar_max_real(dminfo, r, rmax)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      real(kind=RKIND), intent(in) :: r
      real(kind=RKIND), intent(out) :: rmax

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(r, rmax, 1, MPI_REALKIND, MPI_MAX, dminfo % comm, mpi_ierr)
#else
      rmax = r
#endif

   end subroutine mpas_dmpar_max_real!}}}

   subroutine mpas_dmpar_sum_int_array(dminfo, nElements, inArray, outArray)!{{{

      implicit none
   
      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nElements
      integer, dimension(nElements), intent(in) :: inArray
      integer, dimension(nElements), intent(out) :: outArray
      
      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(inArray, outArray, nElements, MPI_INTEGERKIND, MPI_SUM, dminfo % comm, mpi_ierr)
#else
      outArray = inArray
#endif

   end subroutine mpas_dmpar_sum_int_array!}}}

   subroutine mpas_dmpar_min_int_array(dminfo, nElements, inArray, outArray)!{{{
   
      implicit none
      
      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nElements
      integer, dimension(nElements), intent(in) :: inArray
      integer, dimension(nElements), intent(out) :: outArray

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(inArray, outArray, nElements, MPI_INTEGERKIND, MPI_MIN, dminfo % comm, mpi_ierr)
#else
      outArray = inArray
#endif

   end subroutine mpas_dmpar_min_int_array!}}}

   subroutine mpas_dmpar_max_int_array(dminfo, nElements, inArray, outArray)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nElements
      integer, dimension(nElements), intent(in) :: inArray
      integer, dimension(nElements), intent(out) :: outArray

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(inArray, outArray, nElements, MPI_INTEGERKIND, MPI_MAX, dminfo % comm, mpi_ierr)
#else
      outArray = inArray
#endif

   end subroutine mpas_dmpar_max_int_array!}}}

   subroutine mpas_dmpar_sum_real_array(dminfo, nElements, inArray, outArray)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nElements
      real(kind=RKIND), dimension(nElements), intent(in) :: inArray
      real(kind=RKIND), dimension(nElements), intent(out) :: outArray

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(inArray, outArray, nElements, MPI_REALKIND, MPI_SUM, dminfo % comm, mpi_ierr)
#else
      outArray = inArray
#endif

   end subroutine mpas_dmpar_sum_real_array!}}}

   subroutine mpas_dmpar_min_real_array(dminfo, nElements, inArray, outArray)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nElements
      real(kind=RKIND), dimension(nElements), intent(in) :: inArray
      real(kind=RKIND), dimension(nElements), intent(out) :: outArray

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(inArray, outArray, nElements, MPI_REALKIND, MPI_MIN, dminfo % comm, mpi_ierr)
#else
      outArray = inArray
#endif

   end subroutine mpas_dmpar_min_real_array!}}}

   subroutine mpas_dmpar_max_real_array(dminfo, nElements, inArray, outArray)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nElements
      real(kind=RKIND), dimension(nElements), intent(in) :: inArray
      real(kind=RKIND), dimension(nElements), intent(out) :: outArray

      integer :: mpi_ierr

#ifdef _MPI
      call MPI_Allreduce(inArray, outArray, nElements, MPI_REALKIND, MPI_MAX, dminfo % comm, mpi_ierr)
#else
      outArray = inArray
#endif

   end subroutine mpas_dmpar_max_real_array!}}}

   subroutine mpas_dmpar_scatter_ints(dminfo, nprocs, noutlist, displs, counts, inlist, outlist)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nprocs, noutlist
      integer, dimension(nprocs), intent(in) :: displs, counts
      integer, dimension(:), pointer :: inlist
      integer, dimension(noutlist), intent(inout) :: outlist

#ifdef _MPI
      integer :: mpi_ierr
      
      call MPI_Scatterv(inlist, counts, displs, MPI_INTEGERKIND, outlist, noutlist, MPI_INTEGERKIND, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_scatter_ints!}}}

   subroutine mpas_dmpar_get_index_range(dminfo, &!{{{
                                    global_start, global_end, &
                                    local_start, local_end)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: global_start, global_end
      integer, intent(out) :: local_start, local_end

      local_start = nint(real(dminfo % my_proc_id) * real(global_end - global_start + 1) / real(dminfo % nprocs)) + 1
      local_end   = nint(real(dminfo % my_proc_id + 1) * real(global_end - global_start + 1) / real(dminfo % nprocs)) 

   end subroutine mpas_dmpar_get_index_range!}}}

   subroutine mpas_dmpar_compute_index_range(dminfo, &!{{{
                                        local_start, local_end, &
                                        global_start, global_end)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: local_start, local_end
      integer, intent(inout) :: global_start, global_end

      integer :: n
      integer :: mpi_ierr

      n = local_end - local_start + 1

      if (dminfo % my_proc_id == 0) then
         global_start = 1
         global_end = global_start + n - 1
         
#ifdef _MPI
      else if (dminfo % my_proc_id == dminfo % nprocs - 1) then
         call MPI_Recv(global_start, 1, MPI_INTEGERKIND, dminfo % my_proc_id - 1, 0, dminfo % comm, MPI_STATUS_IGNORE, mpi_ierr)
         global_end = global_start + n - 1

      else
         call MPI_Recv(global_start, 1, MPI_INTEGERKIND, dminfo % my_proc_id - 1, 0, dminfo % comm, MPI_STATUS_IGNORE, mpi_ierr)
         global_end = global_start + n
         call MPI_Send(global_end, 1, MPI_INTEGERKIND, dminfo % my_proc_id + 1, 0, dminfo % comm, mpi_ierr)
         global_end = global_end - 1
#endif

      end if
      
   
   end subroutine mpas_dmpar_compute_index_range!}}}

   ! ----- NEW ROUTINES BELOW ----- !

subroutine mpas_dmpar_get_exch_list(haloLayer, ownedListField, neededListField, offsetListField, ownedLimitField)!{{{

      implicit none

      integer, intent(in) :: haloLayer
      type (field1dInteger), pointer :: ownedListField, neededListField
      type (field0dInteger), pointer, optional :: offsetListField
      type (field0dInteger), pointer, optional :: ownedLimitField

      type (dm_info), pointer :: dminfo

      integer :: i, j, k, kk, iBlock
      integer :: totalSize, nMesgRecv, nMesgSend, recvNeighbor, sendNeighbor, currentProc, offset
      integer :: totalSent, totalRecv
      integer, allocatable, dimension(:) :: numToSend, numToRecv
      integer, allocatable, dimension(:) :: ownedList, ownedListIndex, ownedBlock, neededList, neededListIndex, neededBlock
      integer, allocatable, dimension(:) :: offsetList, ownedLimitList
      integer, allocatable, dimension(:,:) :: ownedListSorted, ownedBlockSorted, recipientList
      integer, allocatable, dimension(:) :: ownerListIn, ownerListOut
      integer, allocatable, dimension(:) :: packingOrder
      type (mpas_exchange_list), pointer :: exchListPtr, exchListPtr2
      type (field1dInteger), pointer :: fieldCursor, fieldCursor2
      type (field0dInteger), pointer :: offsetCursor, ownedLimitCursor
      integer :: nOwnedBlocks, nNeededBlocks
      integer :: nOwnedList, nNeededList
      integer :: mpi_ierr, mpi_rreq, mpi_sreq

      type (hashtable) :: neededHash
      integer :: nUniqueNeededList
      integer, dimension(:,:), pointer :: uniqueSortedNeededList


      !
      ! *** NB: This code assumes that block % blockID values are local block IDs and are in the range [1, numBlocks]
      !         where numBlocks is the number of blocks owned by each task
      !


      ! For the ownedListField:
      !    - ownedList contains a list of the global indices owned by all blocks
      !    - ownedListIndex contains a list of the block-local indices of the global indices owned by all blocks
      !    - ownedBlock contains the local block ID associated with each index 
      !
      ! Example:
      !    ownedList      := ( 21 13 15 01 05 06 33 42 44 45 )     ! Global indices from all blocks on this task
      !    ownedListIndex := (  1  2  3  4  1  2  3  4  5  6 )     ! Local  indices of global indices on each block
      !    ownedBlock     := (  1  1  1  1  2  2  2  2  2  2 )     ! Local  indices of global indices on each block
      !
    
      ! For the neededListField:
      !    similar to the owneListField...

      dminfo => ownedListField % block % domain % dminfo

      ! 
      ! Determine total number of owned blocks on this task
      ! 
      nOwnedBlocks = 0
      fieldCursor => ownedListField
      do while (associated(fieldCursor))
        nOwnedBlocks = nOwnedBlocks + 1
        if(associated(fieldCursor % sendList % halos(haloLayer) % exchList)) then
          call mpas_dmpar_destroy_exchange_list(fieldCursor % sendList % halos(haloLayer) % exchList)
        end if

        if(associated(fieldCursor % copyList % halos(haloLayer) % exchList)) then
          call mpas_dmpar_destroy_exchange_list(fieldCursor % copyList % halos(haloLayer) % exchList)
        end if
        fieldCursor => fieldCursor % next
      end do

      !
      ! Determine total number of needed indices on this task
      !
      nNeededList = 0
      nNeededBlocks = 0
      fieldCursor => neededListField
      do while (associated(fieldCursor))
        nNeededBlocks = nNeededBlocks + 1
        nNeededList = nNeededList + fieldCursor % dimSizes(1)
        if(associated(fieldCursor % recvList % halos(haloLayer) % exchList)) then
          call mpas_dmpar_destroy_exchange_list(fieldCursor % recvList % halos(haloLayer) % exchList)
        end if

        fieldCursor => fieldCursor % next
      end do

      !
      ! Determine unique list of needed elements.
      !
      nUniqueNeededList = 0
      call mpas_hash_init(neededHash)
      fieldCursor => neededListField
      do while (associated(fieldCursor))
        do i = 1, fieldCursor % dimSizes(1)
          if(.not. mpas_hash_search(neededHash, fieldCursor % array(i))) then
            nUniqueNeededList = nUniqueNeededList + 1
            call mpas_hash_insert(neededHash, fieldCursor % array(i))
          end if
        end do
        fieldCursor => fieldCursor % next
      end do

      kk = mpas_hash_size(neededHash)

      nUniqueNeededList = mpas_hash_size(neededHash)
      allocate(uniqueSortedNeededList(2,nUniqueNeededList))
      allocate(packingOrder(nUniqueNeededList))
      call mpas_hash_destroy(neededHash)
      call mpas_hash_init(neededHash)

      j = 0
      fieldCursor => neededListField
      do while (associated(fieldCursor))
        do i = 1, fieldCursor % dimSizes(1)
          if(.not. mpas_hash_search(neededHash, fieldCursor % array(i))) then
            j = j +1
            uniqueSortedNeededList(1, j) = fieldCursor % array(i)
            uniqueSortedNeededList(2, j) = fieldCursor % block % localBlockID
            call mpas_hash_insert(neededHash, fieldCursor % array(i))
          end if
        end do
        fieldCursor => fieldCursor % next
      end do

      kk = mpas_hash_size(neededHash)

      call mpas_hash_destroy(neededHash)
      call mpas_quicksort(nUniqueNeededList, uniqueSortedNeededList)

      !
      ! Get list of index offsets for all blocks
      !
      allocate(offsetList(nNeededBlocks))
      if (present(offsetListField)) then
        offsetCursor => offsetListField
        do while (associated(offsetCursor))
          offsetList(offsetCursor % block % localBlockID+1) = offsetCursor % scalar
          offsetCursor => offsetCursor % next
        end do
      else
        offsetList(:) = 0
      end if 

      !
      ! Get list of bounds limit for owned elements
      ! 
      allocate(ownedLimitList(nOwnedBlocks))
      if(present(ownedLimitField)) then
        ownedLimitCursor => ownedLimitField
        do while(associated(ownedLimitCursor))
          ownedLimitList(ownedLimitCursor % block % localBlockID+1) = ownedLimitCursor % scalar
          ownedLimitCursor => ownedLimitCursor % next
        end do
      else
        fieldCursor => ownedListField
        do while(associated(fieldCursor))
          ownedLimitList(fieldCursor % block % localBlockID+1) = fieldCursor % dimSizes(1)
          fieldCursor => fieldCursor % next
        end do
      end if

      ! 
      ! Determine total number of owned indices on this task, and 
      !   initialize output send and recv lists for ownedListField
      ! 
      nOwnedList = 0
      fieldCursor => ownedListField
      do while (associated(fieldCursor))
        iBlock = fieldcursor % block % localBlockID + 1
        nOwnedList = nOwnedList + ownedLimitList(iBlock)
        fieldCursor => fieldCursor % next
      end do

#ifdef _MPI
      !
      ! Gather list of all owned indices and their associated blocks on this task
      !
      allocate(ownedList(nOwnedList))
      allocate(ownedBlock(nOwnedList))
      ownedBlock = -1
      ownedList = -1
      fieldCursor => ownedListField
      i = 1
      do while (associated(fieldCursor))
        iBlock = fieldCursor % block % localBlockID + 1
        ownedList(i:i+ownedLimitList(iBlock)-1) = fieldCursor % array(1:ownedLimitList(iBlock))
        ownedBlock(i:i+ownedLimitList(iBlock)-1) = fieldCursor % block % localBlockID
        i = i + ownedLimitList(iBlock)
        fieldCursor => fieldCursor % next
      end do

      !
      ! Gather list of all needed indices and their associated blocks on this task
      !
      allocate(neededList(nNeededList))
      allocate(neededBlock(nNeededList))
      fieldCursor => neededListField
      i = 1
      do while (associated(fieldCursor))
        neededList(i:i+fieldCursor % dimSizes(1)-1) = fieldCursor % array(:)
        neededBlock(i:i+fieldCursor % dimSizes(1)-1) = fieldCursor % block % localBlockID
        i = i + fieldCursor % dimSizes(1)
        fieldCursor => fieldCursor % next
      end do

      !
      ! Obtain sorted list of global indices owned by this task and the associated local indices and block IDs
      !
      allocate(ownedListIndex(nOwnedList))
      allocate(ownedListSorted(2,nOwnedList))
      allocate(recipientList(2,nOwnedList))
      j = 1
      k = 1
      do i=1,nOwnedList
        ownedListSorted(1,i) = ownedList(i)
        if (i > 1) then
          if(ownedBlock(i) /= ownedBlock(i-1)) k = 1
        end if
        ownedListIndex(i) = k
        ownedListSorted(2,i) = j
        j = j + 1
        k = k + 1
      end do
      call mpas_quicksort(nOwnedList, ownedListSorted)

      allocate(ownedBlockSorted(2,nOwnedList))
      do i=1,nOwnedList
        ownedBlockSorted(1,i) = ownedList(i)
        ownedBlockSorted(2,i) = ownedBlock(i)
      end do
      call mpas_quicksort(nOwnedList, ownedBlockSorted)


      allocate(neededListIndex(nNeededList))
      j = 1
      do i=1,nNeededList
        if (i > 1) then 
          if(neededBlock(i) /= neededBlock(i-1)) j = 1
        end if
        neededListIndex(i) = j
        j = j + 1
      end do

      !
      ! Set totalSize to the maximum number of items in any task's needed list
      !
      call MPI_Allreduce(nUniqueNeededList, totalSize, 1, MPI_INTEGERKIND, MPI_MAX, dminfo % comm, mpi_ierr)

      allocate(ownerListIn(totalSize))
      allocate(ownerListOut(totalSize))

      nMesgSend = nUniqueNeededList
      nMesgRecv = nUniqueNeededList
      ownerListOut(1:nUniqueNeededList) = uniqueSortedNeededList(1,1:nUniqueNeededList)

      recvNeighbor = mod(dminfo % my_proc_id + dminfo % nprocs - 1, dminfo % nprocs)
      sendNeighbor = mod(dminfo % my_proc_id + 1, dminfo % nprocs)

      allocate(numToSend(nOwnedBlocks))
      allocate(numToRecv(nNeededBlocks))

      ! Initial send of data to neighbors.
      if(dminfo % nProcs == 1) then
        ownerListIn = ownerListOut
      else
        call MPI_Irecv(nMesgRecv, 1, MPI_INTEGERKIND, recvNeighbor, recvNeighbor, dminfo % comm, mpi_rreq, mpi_ierr)
        call MPI_Isend(nMesgSend, 1, MPI_INTEGERKIND, sendNeighbor, dminfo % my_proc_id, dminfo % comm, mpi_sreq, mpi_ierr)
        call MPI_Wait(mpi_rreq, MPI_STATUS_IGNORE, mpi_ierr)
        call MPI_Wait(mpi_sreq, MPI_STATUS_IGNORE, mpi_ierr)
        call MPI_Irecv(ownerListIn, nMesgRecv, MPI_INTEGERKIND, recvNeighbor, recvNeighbor, dminfo % comm, mpi_rreq, mpi_ierr)
        call MPI_Isend(ownerListOut, nMesgSend, MPI_INTEGERKIND, sendNeighbor, dminfo % my_proc_id, dminfo % comm, mpi_sreq, mpi_ierr)
        call MPI_Wait(mpi_rreq, MPI_STATUS_IGNORE, mpi_ierr)
        call MPI_Wait(mpi_sreq, MPI_STATUS_IGNORE, mpi_ierr)
      end if

      ! 
      ! For each processor (not including ourself), mark the indices that we will provide to
      !    that processor in ownerListOut, and build a send list for that processor if we
      !    do need to send any indices
      ! 
      do i=2, dminfo % nprocs
        recipientList = -1
        numToSend(:) = 0
        totalSent = 0

        currentProc = mod(dminfo % my_proc_id + dminfo % nprocs - i + 1, dminfo % nprocs)
        do j=1,nMesgRecv
          if (ownerListIn(j) > 0) then
            k = mpas_binary_search(ownedListSorted, 2, 1, nOwnedList, ownerListIn(j))
            if (k <= nOwnedList) then
              iBlock = ownedBlock(ownedListSorted(2,k)) + 1
              numToSend(iBlock) = numToSend(iBlock) + 1
              totalSent = totalSent + 1

              ! recipientList(1,:) represents the index in the srcList to place this data
              recipientList(1,ownedListSorted(2,k)) = numToSend(iBlock)
              ! recipientList(2,:) represnets the index in the buffer to place this data
              recipientList(2,ownedListSorted(2,k)) = totalSent

              ownerListOut(j) = -1 * dminfo % my_proc_id
            else
              ownerListOut(j) = ownerListIn(j)
            end if
          else
            ownerListOut(j) = ownerListIn(j)
          end if
        end do

        fieldCursor => ownedListField
        do while (associated(fieldCursor))
          iBlock = fieldCursor % block % localBlockID + 1

          if (numToSend(iBlock) > 0) then
            ! Find end of send list
            if(.not.associated(fieldCursor % sendList % halos(haloLayer) % exchList)) then
              allocate(fieldCursor % sendList % halos(haloLayer) % exchList)
              exchListPtr => fieldCursor % sendList % halos(haloLayer) % exchList
              nullify(exchListPtr % next)
            else
              exchListPtr => fieldCursor % sendList % halos(haloLayer) % exchList
              exchListPtr2 => fieldCursor % sendList % halos(haloLayer) % exchList % next
              do while(associated(exchListPtr2))
                exchListPtr => exchListPtr % next
                exchListPtr2 => exchListPtr % next
              end do

              allocate(exchListPtr % next)
              exchListPtr => exchListPtr % next
              nullify(exchListPtr % next)
            end if

            exchListPtr % endPointID = currentProc
            exchListPtr % nlist = numToSend(iBlock)
            allocate(exchListPtr % srcList(numToSend(iBlock)))
            allocate(exchListPtr % destList(numToSend(iBlock)))
            exchListPtr % srcList = -1
            exchListPtr % destList = -1

            kk = 1
            do j=1,nOwnedList
              if (recipientList(1,j) /= -1) then
                if(ownedBlock(j) == fieldCursor % block % localBlockID) then
                  exchListPtr % srcList(recipientList(1,j)) = ownedListIndex(j)
                  exchListPtr % destList(recipientList(1,j)) = recipientList(2,j)
                  kk = kk + 1
                end if
              end if
            end do
          end if

          fieldCursor => fieldCursor % next
        end do

        nMesgSend = nMesgRecv
        call MPI_Irecv(nMesgRecv, 1, MPI_INTEGERKIND, recvNeighbor, recvNeighbor, dminfo % comm, mpi_rreq, mpi_ierr)
        call MPI_Isend(nMesgSend, 1, MPI_INTEGERKIND, sendNeighbor, dminfo % my_proc_id, dminfo % comm, mpi_sreq, mpi_ierr)
        call MPI_Wait(mpi_rreq, MPI_STATUS_IGNORE, mpi_ierr)
        call MPI_Wait(mpi_sreq, MPI_STATUS_IGNORE, mpi_ierr)
        call MPI_Irecv(ownerListIn, nMesgRecv, MPI_INTEGERKIND, recvNeighbor, recvNeighbor, dminfo % comm, mpi_rreq, mpi_ierr)
        call MPI_Isend(ownerListOut, nMesgSend, MPI_INTEGERKIND, sendNeighbor, dminfo % my_proc_id, dminfo % comm, mpi_sreq, mpi_ierr)
        call MPI_Wait(mpi_rreq, MPI_STATUS_IGNORE, mpi_ierr)
        call MPI_Wait(mpi_sreq, MPI_STATUS_IGNORE, mpi_ierr)
      end do

      !
      ! With our needed list returned to us, build receive lists based on which indices were
      !    marked by other tasks
      !
      do i=0, dminfo % nprocs - 1
        if(i == dminfo % my_proc_id) cycle

        numToRecv(:) = 0
        packingOrder = 0

        k = 0
        do j=1,nUniqueNeededList
          if (ownerListIn(j) == -i) then
            k = k + 1
            packingOrder(j) = k
          end if
        end do

        fieldCursor => neededListField
        do while (associated(fieldCursor))
          do j = 1, fieldCursor % dimSizes(1)
            k = mpas_binary_search(uniqueSortedNeededList, 2, 1, nUniqueNeededList, fieldCursor % array(j))
            if(k <= nUniqueNeededList) then
              if(ownerListIn(k) == -i) then
                iBlock = fieldCursor % block % localBlockID + 1
                numToRecv(iBlock) = numToRecv(iBlock) + 1
              end if
            end if
          end do
          fieldCursor => fieldCursor % next
        end do

        fieldCursor => neededListField
        do while (associated(fieldCursor))
          iBlock = fieldCursor % block % localBlockID + 1

          if (numToRecv(iBlock) > 0) then
            if(.not.associated(fieldCursor % recvList % halos(haloLayer) % exchList)) then
              allocate(fieldCursor % recvList % halos(haloLayer) % exchList)
              exchListPtr => fieldCursor % recvList % halos(haloLayer) % exchList
              nullify(exchListPtr % next)
            else
              ! Find end of recv list
              exchListPtr => fieldCursor % recvList % halos(haloLayer) % exchList
              exchListPtr2 => fieldCursor % recvList % halos(haloLayer) % exchList % next
              do while(associated(exchListPtr2))
                exchListPtr => exchListPtr % next
                exchListPtr2 => exchListPtr % next
              end do

              allocate(exchListPtr % next)
              exchListPtr => exchListPtr % next
              nullify(exchListPtr % next)
            end if

            exchListPtr % endPointID = i
            exchListPtr % nlist = numToRecv(iBlock)
            allocate(exchListPtr % srcList(exchListPtr % nList))
            allocate(exchListPtr % destList(exchListPtr % nList))
            exchListPtr % srcList = -1
            exchListPtr % destList = -1

            kk = 0
            do j=1,fieldCursor % dimSizes(1)
              k = mpas_binary_search(uniqueSortedNeededList, 2, 1, nUniqueNeededList, fieldCursor % array(j))
              if(k <= nUniqueNeededList) then
                if (ownerListIn(k) == -i) then
                  kk = kk + 1
                  exchListPtr % srcList(kk) = packingOrder(k)
                  exchListPtr % destList(kk) = j + offsetList(iBlock)
                end if
              end if
            end do
          end if

          fieldCursor => fieldCursor % next
        end do
      end do

      !
      ! Free up memory
      !
      deallocate(numToSend)
      deallocate(numToRecv)
      deallocate(ownedList)
      deallocate(ownedListIndex)
      deallocate(ownedBlock)
      deallocate(neededList)
      deallocate(neededListIndex)
      deallocate(neededBlock)
      deallocate(ownedListSorted)
      deallocate(ownedBlockSorted)
      deallocate(recipientList)
      deallocate(ownerListIn)
      deallocate(ownerListOut)
      deallocate(uniqueSortedNeededList)
      deallocate(packingOrder)
#endif

      ! Build Copy Lists
      allocate(numToSend(1))
      fieldCursor => ownedListField
      do while (associated(fieldCursor))
        iBlock = fieldCursor % block % localBlockID + 1
        nOwnedList = ownedLimitList(iBlock)
        allocate(ownedListSorted(2, nOwnedList))
        allocate(recipientList(2, nOwnedList))

        do i = 1, nOwnedList
          ownedListSorted(1, i) = fieldCursor % array(i)
          ownedListSorted(2, i) = i
        end do

        call mpas_quicksort(nOwnedList, ownedListSorted)

        fieldCursor2 => neededListField
        do while(associated(fieldCursor2))
          if(associated(fieldCursor, fieldCursor2)) then
            fieldCursor2 => fieldCursor2 % next
            cycle
          end if

          numToSend = 0
          recipientList = -1

          do i = 1, fieldCursor2 % dimSizes(1)
            k = mpas_binary_search(ownedListSorted, 2, 1, nOwnedList, fieldCursor2 % array(i))
            if (k <= nOwnedList) then
              numToSend(1) = numToSend(1) + 1
              ! recipientList(1,:) represents the needed block id
              recipientList(1,ownedListSorted(2,k)) = fieldCursor2 % block % localBlockID
              ! recipientList(2,:) represnets the index in the buffer to place this data
              recipientList(2,ownedListSorted(2,k)) = i
            end if
          end do
          
          if(numToSend(1) > 0) then
            if(.not.associated(fieldCursor % copyList % halos(haloLayer) % exchList)) then
              allocate(fieldCursor % copyList % halos(haloLayer) % exchList)
              exchListPtr => fieldCursor % copyList % halos(haloLayer) % exchList
              nullify(exchListPtr % next)
            else
              ! Find end of copy list
              exchListPtr => fieldCursor % copyList % halos(haloLayer) % exchList
              exchListPtr2 => fieldCursor % copyList % halos(haloLayer) % exchList % next
              do while(associated(exchListPtr2))
                exchListPtr => exchListPtr % next
                exchListPtr2 => exchListPtr % next
              end do

              allocate(exchListPtr % next)
              exchListPtr => exchListPtr % next
              nullify(exchListPtr % next)
            end if
    
            exchListPtr % endPointID = fieldCursor2 % block % localBlockID
            exchListPtr % nlist = numToSend(1)
            allocate(exchListPtr % srcList(numToSend(1)))
            allocate(exchListPtr % destList(numToSend(1)))
            exchListPtr % srcList = -1
            exchListPtr % destList = -1

            kk = 1
            do j=1,nOwnedList
             if(recipientList(1,j) == fieldCursor2 % block % localBlockID) then
               exchListPtr % srcList(kk) = j
               exchListPtr % destList(kk) = recipientList(2,j) + offSetList(fieldCursor2 % block % localBlockID+1)
               kk = kk + 1
             end if
            end do
          end if
          fieldCursor2 => fieldCursor2 % next
        end do

        deallocate(recipientList)
        deallocate(ownedListSorted)
        fieldCursor => fieldCursor % next
      end do
      deallocate(numToSend)
      deallocate(offSetList)

   end subroutine mpas_dmpar_get_exch_list!}}}


   subroutine mpas_dmpar_alltoall_field1d_integer(fieldIn, fieldout, haloLayersIn)!{{{

     implicit none

     type (field1dInteger), pointer :: fieldIn
     type (field1dInteger), pointer :: fieldOut
     integer, dimension(:), pointer, optional :: haloLayersIn

     type (field1dInteger), pointer :: fieldInPtr, fieldOutPtr
     type (mpas_exchange_list), pointer :: exchListPtr, exchListPtr2
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i
     integer :: nHaloLayers
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo

     if(present(haloLayersIn)) then
       nHaloLayers = size(haloLayersIn)
       allocate(haloLayers(nHaloLayers))
       do iHalo = 1, nHaloLayers
         haloLayers(iHalo) = haloLayersIn(iHalo)
       end do
     else
       nHaloLayers = size(fieldIn % sendList % halos)
       allocate(haloLayers(nHaloLayers))
       do iHalo = 1, nHaloLayers
         haloLayers(iHalo) = iHalo
       end do
     end if

#ifdef _MPI
     nullify(sendList)
     nullify(recvList)

     ! Setup recieve lists.
     do iHalo = 1, nHaloLayers
       fieldOutPtr => fieldOut
       do while(associated(fieldOutPtr))
         exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           comm_list_found = .false.
  
           ! Search for an already created commList to this processor.
           commListPtr => recvList
           do while(associated(commListPtr))
             if(commListPtr % procID == exchListPtr % endPointID) then
               comm_list_found = .true.
               exit
             end if
  
             commListPtr => commListPtr % next
           end do
  
           ! If no comm list exists, create a new one.
           if(.not. comm_list_found) then
             if(.not.associated(recvList)) then
               allocate(recvList)
               nullify(recvList % next)
               commListPtr => recvList
             else
               commListPtr => recvList
               commListPtr2 => commListPtr % next
               do while(associated(commListPtr2))
                 commListPtr => commListPtr % next
                 commListPtr2 => commListPtr % next
               end do

               allocate(commListPtr % next)
               commListPtr => commListPtr % next
               nullify(commListPtr % next)
             end if
  
             commListPtr % procID = exchListPtr % endPointID
             commListPtr % nList = 0
           end if

           exchListPtr => exchListPtr % next
         end do
  
         fieldOutPtr => fieldOutPtr % next
       end do
     end do

     ! Determine size of receive list buffers.
     commListPtr => recvList
     do while(associated(commListPtr))
       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               nAdded = max(nAdded, maxval(exchListPtr % srcList))
             end if
             exchListPtr => exchListPtr % next
           end do
  
           fieldOutPtr => fieldOutPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do
       commListPtr % nList = nAdded

       commListPtr => commListPtr % next
     end do

     ! Allocate buffers for recieves, and initiate mpi_irecv calls.
     commListPtr => recvList
     do while(associated(commListPtr))
       allocate(commListPtr % ibuffer(commListPtr % nList))
       nullify(commListPtr % rbuffer)
       commListPtr % ibuffer = 0
       call MPI_Irecv(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)
       commListPtr => commListPtr % next
     end do

     ! Setup send lists, and determine the size of their buffers.
     do iHalo = 1, nHaloLayers
       fieldInPtr => fieldIn
       do while(associated(fieldInPtr))
         exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           comm_list_found = .false.
  
           ! Search for an already created commList to this processor.
           commListPtr => sendList
           do while(associated(commListPtr))
             if(commListPtr % procID == exchListPtr % endPointID) then
               commListPtr % nList = commListPtr % nList + exchListPtr % nList 
               comm_list_found = .true.
               exit
             end if
  
             commListPtr => commListPtr % next
           end do
  
           ! If no comm list exists, create a new one.
           if(.not. comm_list_found) then
             if(.not.associated(sendList)) then
               allocate(sendList)
               nullify(sendList % next)
               commListPtr => sendList
             else
               commListPtr => sendList
               commListPtr2 => commListPtr % next
               do while(associated(commListPtr2))
                 commListPtr => commListPtr % next
                 commListPtr2 => commListPtr % next
               end do
    
               allocate(commListPtr % next)
               commListPtr => commListPtr % next
               nullify(commListPtr % next)
             end if
             commListPtr % procID = exchListPtr % endPointID
             commListPtr % nList = exchListPtr % nList
           end if
  
           exchListPtr => exchListPtr % next
         end do
  
         fieldInPtr => fieldInPtr % next
       end do
     end do

     ! Allocate sendLists, copy data into buffer, and initiate mpi_isends
     commListPtr => sendList
     do while(associated(commListPtr))
       allocate(commListPtr % ibuffer(commListPtr % nList))
       nullify(commListPtr % rbuffer)
       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               do i = 1, exchListPtr % nList
                 iBuffer = exchListPtr % destList(i) + bufferOffset
                 commListPtr % ibuffer(iBuffer) = fieldInPtr % array(exchListPtr % srcList(i))
                 nAdded = nAdded + 1
               end do
             end if
  
             exchListPtr => exchListPtr % next
           end do
  
           fieldInPtr => fieldInPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do

       call MPI_Isend(commListPtr % ibuffer, commListPtr % nlist, MPI_INTEGERKIND, &
                      commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)

       commListPtr => commListPtr % next
     end do

#endif     

     ! Handle Local Copies. Only local copies if no MPI
     do iHalo = 1, nHaloLayers
       fieldInPtr => fieldIn
       do while(associated(fieldInPtr))
         exchListPtr => fieldInPtr % copyList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             if(exchListPtr % endPointID == fieldOutPtr % block % localBlockID) then
               do i = 1, exchListPtr % nList
                 fieldOutPtr % array(exchListPtr % destList(i)) = fieldInPtr % array(exchListPtr % srcList(i))
               end do
             end if
             fieldOutPtr => fieldOutPtr % next
           end do
  
           exchListPtr => exchListPtr % next
         end do
         fieldInPtr => fieldInPtr % next
       end do
     end do

#ifdef _MPI
     ! Wait for MPI_Irecv's to finish, and unpack data.
     commListPtr => recvList
     do while(associated(commListPtr))
       call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               do i = 1, exchListPtr % nList
                 iBuffer = exchListPtr % srcList(i) + bufferOffset
                 fieldOutPtr % array(exchListPtr % destList(i)) = commListPtr % ibuffer(iBuffer)
               end do
               nAdded = max(nAdded, maxval(exchListPtr % srcList))
             end if
             exchListPtr => exchListPtr % next
           end do
  
           fieldOutPtr => fieldOutPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do

       commListPtr => commListPtr % next
     end do

     ! Wait for MPI_Isend's to finish.
     commListPtr => sendList
     do while(associated(commListPtr))
       call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
       commListPtr => commListPtr % next
     end do

     ! Destroy commLists.
     call mpas_dmpar_destroy_communication_list(sendList)
     call mpas_dmpar_destroy_communication_list(recvList)
#endif

     deallocate(haloLayers)

   end subroutine mpas_dmpar_alltoall_field1d_integer!}}}

   subroutine mpas_dmpar_alltoall_field2d_integer(fieldIn, fieldout, haloLayersIn)!{{{

     implicit none

     type (field2dInteger), pointer :: fieldIn
     type (field2dInteger), pointer :: fieldOut
     integer, dimension(:), pointer, optional :: haloLayersIn

     type (field2dInteger), pointer :: fieldInPtr, fieldOutPtr
     type (mpas_exchange_list), pointer :: exchListPtr, exchListPtr2
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i, j
     integer :: nHaloLayers
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo

     if(present(haloLayersIn)) then
       nHaloLayers = size(haloLayersIn)
       allocate(haloLayers(nHaloLayers))
       do iHalo = 1, nHaloLayers
         haloLayers(iHalo) = haloLayersIn(iHalo)
       end do
     else
       nHaloLayers = size(fieldIn % sendList % halos)
       allocate(haloLayers(nHaloLayers))
       do iHalo = 1, nHaloLayers
         haloLayers(iHalo) = iHalo
       end do
     end if

#ifdef _MPI
     nullify(sendList)
     nullify(recvList)

     ! Setup recieve lists
     do iHalo = 1, nHaloLayers
       fieldOutPtr => fieldOut
       do while(associated(fieldOutPtr))
         exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           comm_list_found = .false.
  
           ! Search for an already created commList to this processor.
           commListPtr => recvList
           do while(associated(commListPtr))
             if(commListPtr % procID == exchListPtr % endPointID) then
               comm_list_found = .true.
               exit
             end if
  
             commListPtr => commListPtr % next
           end do
  
           ! If no comm list exists, create a new one.
           if(.not. comm_list_found) then
             if(.not.associated(recvList)) then
               allocate(recvList)
               nullify(recvList % next)
               commListPtr => recvList
             else
               commListPtr => recvList
               commListPtr2 => commListPtr % next
               do while(associated(commListPtr2))
                 commListPtr => commListPtr % next
                 commListPtr2 => commListPtr % next
               end do

               allocate(commListPtr % next)
               commListPtr => commListPtr % next
               nullify(commListPtr % next)
             end if
  
             commListPtr % procID = exchListPtr % endPointID
           end if

           exchListPtr => exchListPtr % next
         end do
  
         fieldOutPtr => fieldOutPtr % next
       end do
     end do

     ! Determine size of receive list buffers.
     commListPtr => recvList
     do while(associated(commListPtr))
       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1))
             end if
             exchListPtr => exchListPtr % next
           end do
  
           fieldOutPtr => fieldOutPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do
       commListPtr % nList = bufferOffset

       commListPtr => commListPtr % next
     end do

     ! Allocate buffers for recieves, and initiate mpi_irecv calls.
     commListPtr => recvList
     do while(associated(commListPtr))
       allocate(commListPtr % ibuffer(commListPtr % nList))
       nullify(commListPtr % rbuffer)
       call MPI_Irecv(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)
       commListPtr => commListPtr % next
     end do

     ! Setup send lists, and determine the size of their buffers.
     do iHalo = 1, nHaloLayers
       fieldInPtr => fieldIn
       do while(associated(fieldInPtr))
         exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           comm_list_found = .false.
  
           ! Search for an already created commList to this processor.
           commListPtr => sendList
           do while(associated(commListPtr))
             if(commListPtr % procID == exchListPtr % endPointID) then
               commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldInPtr % dimSizes(1)
               comm_list_found = .true.
               exit
             end if
  
             commListPtr => commListPtr % next
           end do
  
           ! If no comm list exists, create a new one.
           if(.not. comm_list_found) then
             if(.not.associated(sendList)) then
               allocate(sendList)
               nullify(sendList % next)
               commListPtr => sendList
             else
               commListPtr => sendList
               commListPtr2 => commListPtr % next
               do while(associated(commListPtr2))
                 commListPtr => commListPtr % next
                 commListPtr2 => commListPtr % next
               end do
    
               allocate(commListPtr % next)
               commListPtr => commListPtr % next
               nullify(commListPtr % next)
             end if
             commListPtr % procID = exchListPtr % endPointID
             commListPtr % nList = exchListPtr % nList * fieldInPtr % dimSizes(1)
           end if
  
           exchListPtr => exchListPtr % next
         end do
  
         fieldInPtr => fieldInPtr % next
       end do
     end do

     ! Allocate sendLists, copy data into buffer, and initiate mpi_isends
     commListPtr => sendList
     do while(associated(commListPtr))
       allocate(commListPtr % ibuffer(commListPtr % nList))
       nullify(commListPtr % rbuffer)
       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               do i = 1, exchListPtr % nList
                 do j = 1, fieldInPtr % dimSizes(1)
                   iBuffer = (exchListPtr % destList(i)-1) * fieldInPtr % dimSizes(1) + j + bufferOffset
                   commListPtr % ibuffer(iBuffer) = fieldInPtr % array(j, exchListPtr % srcList(i))
                   nAdded = nAdded + 1
                 end do
               end do
             end if
  
             exchListPtr => exchListPtr % next
           end do
  
           fieldInPtr => fieldInPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do

       call MPI_Isend(commListPtr % ibuffer, commListPtr % nlist, MPI_INTEGERKIND, &
                      commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)

       commListPtr => commListPtr % next
     end do
#endif     

     ! Handle Local Copies. Only local copies if no MPI
     do iHalo = 1, nHaloLayers
       fieldInPtr => fieldIn
       do while(associated(fieldInPtr))
         exchListPtr => fieldInPtr % copyList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             if(exchListPtr % endPointID == fieldOutPtr % block % localBlockID) then
               do i = 1, exchListPtr % nList
                 fieldOutPtr % array(:, exchListPtr % destList(i)) = fieldInPtr % array(:, exchListPtr % srcList(i))
               end do
             end if
             fieldOutPtr => fieldOutPtr % next
           end do
  
           exchListPtr => exchListPtr % next
         end do
         fieldInPtr => fieldInPtr % next
       end do
     end do

#ifdef _MPI
     ! Wait for MPI_Irecv's to finish, and unpack data.
     commListPtr => recvList
     do while(associated(commListPtr))
       call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               do i = 1, exchListPtr % nList
                 do j = 1, fieldOutPtr % dimSizes(1)
                   iBuffer = (exchListPtr % srcList(i)-1) * fieldOutPtr % dimSizes(1) + j + bufferOffset
                   fieldOutPtr % array(j, exchListPtr % destList(i)) = commListPtr % ibuffer(iBuffer)
                 end do
               end do
               nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1))
             end if
             exchListPtr => exchListPtr % next
           end do
  
           fieldOutPtr => fieldOutPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do

       commListPtr => commListPtr % next
     end do

     ! Wait for MPI_Isend's to finish.
     commListPtr => sendList
     do while(associated(commListPtr))
       call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
       commListPtr => commListPtr % next
     end do

     ! Destroy commLists.
     call mpas_dmpar_destroy_communication_list(sendList)
     call mpas_dmpar_destroy_communication_list(recvList)
#endif

     deallocate(haloLayers)

   end subroutine mpas_dmpar_alltoall_field2d_integer!}}}

   subroutine mpas_dmpar_alltoall_field3d_integer(fieldIn, fieldout, haloLayersIn)!{{{

     implicit none

     type (field3dInteger), pointer :: fieldIn
     type (field3dInteger), pointer :: fieldOut
     integer, dimension(:), pointer, optional :: haloLayersIn

     type (field3dInteger), pointer :: fieldInPtr, fieldOutPtr
     type (mpas_exchange_list), pointer :: exchListPtr, exchListPtr2
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i, j, k
     integer :: nHaloLayers
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo

     if(present(haloLayersIn)) then
       nHaloLayers = size(haloLayersIn)
       allocate(haloLayers(nHaloLayers))
       do iHalo = 1, nHaloLayers
         haloLayers(iHalo) = haloLayersIn(iHalo)
       end do
     else
       nHaloLayers = size(fieldIn % sendList % halos)
       allocate(haloLayers(nHaloLayers))
       do iHalo = 1, nHaloLayers
         haloLayers(iHalo) = iHalo
       end do
     end if

#ifdef _MPI
     nullify(sendList)
     nullify(recvList)

     ! Setup recieve lists.
     do iHalo = 1, nHaloLayers
       fieldOutPtr => fieldOut
       do while(associated(fieldOutPtr))
         exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           comm_list_found = .false.
  
           ! Search for an already created commList to this processor.
           commListPtr => recvList
           do while(associated(commListPtr))
             if(commListPtr % procID == exchListPtr % endPointID) then
               commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2)
               comm_list_found = .true.
               exit
             end if
  
             commListPtr => commListPtr % next
           end do
  
           ! If no comm list exists, create a new one.
           if(.not. comm_list_found) then
             if(.not.associated(recvList)) then
               allocate(recvList)
               nullify(recvList % next)
               commListPtr => recvList
             else
               commListPtr => recvList
               commListPtr2 => commListPtr % next
               do while(associated(commListPtr2))
                 commListPtr => commListPtr % next
                 commListPtr2 => commListPtr % next
               end do

               allocate(commListPtr % next)
               commListPtr => commListPtr % next
               nullify(commListPtr % next)
             end if
  
             commListPtr % procID = exchListPtr % endPointID
             commListPtr % nList = exchListPtr % nList * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2)
           end if

           exchListPtr => exchListPtr % next
         end do
  
         fieldOutPtr => fieldOutPtr % next
       end do
     end do

     ! Determine size of receive list buffers
     commListPtr => recvList
     do while(associated(commListPtr))
       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2))
             end if
             exchListPtr => exchListPtr % next
           end do
  
           fieldOutPtr => fieldOutPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do
       commListPtr % nList = nAdded

       commListPtr => commListPtr % next
     end do

     ! Allocate buffers for recieves, and initiate mpi_irecv calls.
     commListPtr => recvList
     do while(associated(commListPtr))
       allocate(commListPtr % ibuffer(commListPtr % nList))
       nullify(commListPtr % rbuffer)
       call MPI_Irecv(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)
       commListPtr => commListPtr % next
     end do

     ! Setup send lists, and determine the size of their buffers.
     do iHalo = 1, nHaloLayers
       fieldInPtr => fieldIn
       do while(associated(fieldInPtr))
         exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           comm_list_found = .false.
  
           ! Search for an already created commList to this processor.
           commListPtr => sendList
           do while(associated(commListPtr))
             if(commListPtr % procID == exchListPtr % endPointID) then
               commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2)
               comm_list_found = .true.
               exit
             end if
  
             commListPtr => commListPtr % next
           end do
  
           ! If no comm list exists, create a new one.
           if(.not. comm_list_found) then
             if(.not.associated(sendList)) then
               allocate(sendList)
               nullify(sendList % next)
               commListPtr => sendList
             else
               commListPtr => sendList
               commListPtr2 => commListPtr % next
               do while(associated(commListPtr2))
                 commListPtr => commListPtr % next
                 commListPtr2 => commListPtr % next
               end do
    
               allocate(commListPtr % next)
               commListPtr => commListPtr % next
               nullify(commListPtr % next)
             end if
             commListPtr % procID = exchListPtr % endPointID
             commListPtr % nList = exchListPtr % nList * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2)
           end if
  
           exchListPtr => exchListPtr % next
         end do
  
         fieldInPtr => fieldInPtr % next
       end do
     end do

     ! Allocate sendLists, copy data into buffer, and initiate mpi_isends
     commListPtr => sendList
     do while(associated(commListPtr))
       allocate(commListPtr % ibuffer(commListPtr % nList))
       nullify(commListPtr % rbuffer)
       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               do i = 1, exchListPtr % nList
                 do j = 1, fieldInPtr % dimSizes(2)
                   do k = 1, fieldInPtr % dimSizes(1)
                     iBuffer = (exchListPtr % destList(i)-1) * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2) + (j-1) * fieldInPtr % dimSizes(1) + k + bufferOffset
                     commListPtr % ibuffer(iBuffer) = fieldInPtr % array(k, j, exchListPtr % srcList(i))
                     nAdded = nAdded + 1
                   end do
                 end do
               end do
             end if
  
             exchListPtr => exchListPtr % next
           end do
  
           fieldInPtr => fieldInPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do

       call MPI_Isend(commListPtr % ibuffer, commListPtr % nlist, MPI_INTEGERKIND, &
                      commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)

       commListPtr => commListPtr % next
     end do

#endif     

     ! Handle Local Copies. Only local copies if no MPI
     do iHalo = 1, nHaloLayers
       fieldInPtr => fieldIn
       do while(associated(fieldInPtr))
         exchListPtr => fieldInPtr % copyList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             if(exchListPtr % endPointID == fieldOutPtr % block % localBlockID) then
               do i = 1, exchListPtr % nList
                 fieldOutPtr % array(:, :, exchListPtr % destList(i)) = fieldInPtr % array(:, :, exchListPtr % srcList(i))
               end do
             end if
             fieldOutPtr => fieldOutPtr % next
           end do
  
           exchListPtr => exchListPtr % next
         end do
         fieldInPtr => fieldInPtr % next
       end do
     end do

#ifdef _MPI
     ! Wait for MPI_Irecv's to finish, and unpack data.
     commListPtr => recvList
     do while(associated(commListPtr))
       call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               do i = 1, exchListPtr % nList
                 do j = 1, fieldOutPtr % dimSizes(2)
                   do k = 1, fieldOutPtr % dimSizes(1)
                     iBuffer = (exchListPtr % srcList(i)-1) * fieldOutPtr % dimSizes(2) * fieldOutPtr % dimSizes(1) + (j-1) * fieldOutPtr % dimSizes(1) + k + bufferOffset
                     fieldOutPtr % array(k, j, exchListPtr % destList(i)) = commListPtr % ibuffer(iBuffer)
                   end do
                 end do
               end do
               nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2))
             end if
             exchListPtr => exchListPtr % next
           end do
  
           fieldOutPtr => fieldOutPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do

       commListPtr => commListPtr % next
     end do

     ! Wait for MPI_Isend's to finish.
     commListPtr => sendList
     do while(associated(commListPtr))
       call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
       commListPtr => commListPtr % next
     end do

     ! Destroy commLists.
     call mpas_dmpar_destroy_communication_list(sendList)
     call mpas_dmpar_destroy_communication_list(recvList)
#endif

     deallocate(haloLayers)

   end subroutine mpas_dmpar_alltoall_field3d_integer!}}}

   subroutine mpas_dmpar_alltoall_field1d_real(fieldIn, fieldout, haloLayersIn)!{{{

     implicit none

     type (field1dReal), pointer :: fieldIn
     type (field1dReal), pointer :: fieldOut
     integer, dimension(:), pointer, optional :: haloLayersIn

     type (field1dReal), pointer :: fieldInPtr, fieldOutPtr
     type (mpas_exchange_list), pointer :: exchListPtr, exchListPtr2
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i
     integer :: nHaloLayers
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo

     if(present(haloLayersIn)) then
       nHaloLayers = size(haloLayersIn)
       allocate(haloLayers(nHaloLayers))
       do iHalo = 1, nHaloLayers
         haloLayers(iHalo) = haloLayersIn(iHalo)
       end do
     else
       nHaloLayers = size(fieldIn % sendList % halos)
       allocate(haloLayers(nHaloLayers))
       do iHalo = 1, nHaloLayers
         haloLayers(iHalo) = iHalo
       end do
     end if

#ifdef _MPI
     nullify(sendList)
     nullify(recvList)

     ! Setup recieve lists.
     do iHalo = 1, nHaloLayers
       fieldOutPtr => fieldOut
       do while(associated(fieldOutPtr))
         exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           comm_list_found = .false.
  
           ! Search for an already created commList to this processor.
           commListPtr => recvList
           do while(associated(commListPtr))
             if(commListPtr % procID == exchListPtr % endPointID) then
               comm_list_found = .true.
               exit
             end if
  
             commListPtr => commListPtr % next
           end do
  
           ! If no comm list exists, create a new one.
           if(.not. comm_list_found) then
             if(.not.associated(recvList)) then
               allocate(recvList)
               nullify(recvList % next)
               commListPtr => recvList
             else
               commListPtr => recvList
               commListPtr2 => commListPtr % next
               do while(associated(commListPtr2))
                 commListPtr => commListPtr % next
                 commListPtr2 => commListPtr % next
               end do

               allocate(commListPtr % next)
               commListPtr => commListPtr % next
               nullify(commListPtr % next)
             end if
  
             commListPtr % procID = exchListPtr % endPointID
           end if

           exchListPtr => exchListPtr % next
         end do
  
         fieldOutPtr => fieldOutPtr % next
       end do
     end do

     ! Determine size of receive list buffers
     commListPtr => recvList
     do while(associated(commListPtr))
       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               nAdded = max(nAdded, maxval(exchListPtr % srcList))
             end if
             exchListPtr => exchListPtr % next
           end do
  
           fieldOutPtr => fieldOutPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do
       commListPtr % nList = nAdded

       commListPtr => commListPtr % next
     end do

     ! Allocate buffers for recieves, and initiate mpi_irecv calls.
     commListPtr => recvList
     do while(associated(commListPtr))
       allocate(commListPtr % rbuffer(commListPtr % nList))
       nullify(commListPtr % ibuffer)
       call MPI_Irecv(commListPtr % rbuffer, commListPtr % nList, MPI_realKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)
       commListPtr => commListPtr % next
     end do

     ! Setup send lists, and determine the size of their buffers.
     do iHalo = 1, nHaloLayers
       fieldInPtr => fieldIn
       do while(associated(fieldInPtr))
         exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           comm_list_found = .false.
  
           ! Search for an already created commList to this processor.
           commListPtr => sendList
           do while(associated(commListPtr))
             if(commListPtr % procID == exchListPtr % endPointID) then
               commListPtr % nList = commListPtr % nList + exchListPtr % nList 
               comm_list_found = .true.
               exit
             end if
  
             commListPtr => commListPtr % next
           end do
  
           ! If no comm list exists, create a new one.
           if(.not. comm_list_found) then
             if(.not.associated(sendList)) then
               allocate(sendList)
               nullify(sendList % next)
               commListPtr => sendList
             else
               commListPtr => sendList
               commListPtr2 => commListPtr % next
               do while(associated(commListPtr2))
                 commListPtr => commListPtr % next
                 commListPtr2 => commListPtr % next
               end do
    
               allocate(commListPtr % next)
               commListPtr => commListPtr % next
               nullify(commListPtr % next)
             end if
             commListPtr % procID = exchListPtr % endPointID
             commListPtr % nList = exchListPtr % nList
           end if
  
           exchListPtr => exchListPtr % next
         end do
  
         fieldInPtr => fieldInPtr % next
       end do
     end do

     ! Allocate sendLists, copy data into buffer, and initiate mpi_isends
     commListPtr => sendList
     do while(associated(commListPtr))
       allocate(commListPtr % rbuffer(commListPtr % nList))
       nullify(commListPtr % ibuffer)
       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               do i = 1, exchListPtr % nList
                 iBuffer = exchListPtr % destList(i) + bufferOffset
                 commListPtr % rbuffer(iBuffer) = fieldInPtr % array(exchListPtr % srcList(i))
                 nAdded = nAdded + 1
               end do
             end if
  
             exchListPtr => exchListPtr % next
           end do
  
           fieldInPtr => fieldInPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do

       call MPI_Isend(commListPtr % rbuffer, commListPtr % nlist, MPI_realKIND, &
                      commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)

       commListPtr => commListPtr % next
     end do

#endif     

     ! Handle Local Copies. Only local copies if no MPI
     do iHalo = 1, nHaloLayers
       fieldInPtr => fieldIn
       do while(associated(fieldInPtr))
         exchListPtr => fieldInPtr % copyList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             if(exchListPtr % endPointID == fieldOutPtr % block % localBlockID) then
               do i = 1, exchListPtr % nList
                 fieldOutPtr % array(exchListPtr % destList(i)) = fieldInPtr % array(exchListPtr % srcList(i))
               end do
             end if
             fieldOutPtr => fieldOutPtr % next
           end do
  
           exchListPtr => exchListPtr % next
         end do
         fieldInPtr => fieldInPtr % next
       end do
     end do

#ifdef _MPI
     ! Wait for MPI_Irecv's to finish, and unpack data.
     commListPtr => recvList
     do while(associated(commListPtr))
       call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               do i = 1, exchListPtr % nList
                 iBuffer = exchListPtr % srcList(i) + bufferOffset
                 fieldOutPtr % array(exchListPtr % destList(i)) = commListPtr % rbuffer(iBuffer)
               end do
               nAdded = max(nAdded, maxval(exchListPtr % srcList))
             end if
             exchListPtr => exchListPtr % next
           end do
  
           fieldOutPtr => fieldOutPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do

       commListPtr => commListPtr % next
     end do

     ! Wait for MPI_Isend's to finish.
     commListPtr => sendList
     do while(associated(commListPtr))
       call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
       commListPtr => commListPtr % next
     end do

     ! Destroy commLists.
     call mpas_dmpar_destroy_communication_list(sendList)
     call mpas_dmpar_destroy_communication_list(recvList)
#endif

     deallocate(haloLayers)

   end subroutine mpas_dmpar_alltoall_field1d_real!}}}

   subroutine mpas_dmpar_alltoall_field2d_real(fieldIn, fieldout, haloLayersIn)!{{{

     implicit none

     type (field2dReal), pointer :: fieldIn
     type (field2dReal), pointer :: fieldOut
     integer, dimension(:), pointer, optional :: haloLayersIn

     type (field2dReal), pointer :: fieldInPtr, fieldOutPtr
     type (mpas_exchange_list), pointer :: exchListPtr, exchListPtr2
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i, j
     integer :: nHaloLayers
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo

     if(present(haloLayersIn)) then
       nHaloLayers = size(haloLayersIn)
       allocate(haloLayers(nHaloLayers))
       do iHalo = 1, nHaloLayers
         haloLayers(iHalo) = haloLayersIn(iHalo)
       end do
     else
       nHaloLayers = size(fieldIn % sendList % halos)
       allocate(haloLayers(nHaloLayers))
       do iHalo = 1, nHaloLayers
         haloLayers(iHalo) = iHalo
       end do
     end if

#ifdef _MPI
     nullify(sendList)
     nullify(recvList)

     ! Setup recieve lists, and determine the size of their buffers.
     do iHalo = 1, nHaloLayers
       fieldOutPtr => fieldOut
       do while(associated(fieldOutPtr))
         exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           comm_list_found = .false.
  
           ! Search for an already created commList to this processor.
           commListPtr => recvList
           do while(associated(commListPtr))
             if(commListPtr % procID == exchListPtr % endPointID) then
               comm_list_found = .true.
               exit
             end if
  
             commListPtr => commListPtr % next
           end do
  
           ! If no comm list exists, create a new one.
           if(.not. comm_list_found) then
             if(.not.associated(recvList)) then
               allocate(recvList)
               nullify(recvList % next)
               commListPtr => recvList
             else
               commListPtr => recvList
               commListPtr2 => commListPtr % next
               do while(associated(commListPtr2))
                 commListPtr => commListPtr % next
                 commListPtr2 => commListPtr % next
               end do

               allocate(commListPtr % next)
               commListPtr => commListPtr % next
               nullify(commListPtr % next)
             end if
  
             commListPtr % procID = exchListPtr % endPointID
           end if

           exchListPtr => exchListPtr % next
         end do
  
         fieldOutPtr => fieldOutPtr % next
       end do
     end do

     ! Determine size of receive list buffers.
     commListPtr => recvList
     do while(associated(commListPtr))
       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1))
             end if
             exchListPtr => exchListPtr % next
           end do
  
           fieldOutPtr => fieldOutPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do
       commListPtr % nList = nAdded

       commListPtr => commListPtr % next
     end do

     ! Allocate buffers for recieves, and initiate mpi_irecv calls.
     commListPtr => recvList
     do while(associated(commListPtr))
       allocate(commListPtr % rbuffer(commListPtr % nList))
       nullify(commListPtr % ibuffer)
       call MPI_Irecv(commListPtr % rbuffer, commListPtr % nList, MPI_realKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)
       commListPtr => commListPtr % next
     end do

     ! Setup send lists, and determine the size of their buffers.
     do iHalo = 1, nHaloLayers
       fieldInPtr => fieldIn
       do while(associated(fieldInPtr))
         exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           comm_list_found = .false.
  
           ! Search for an already created commList to this processor.
           commListPtr => sendList
           do while(associated(commListPtr))
             if(commListPtr % procID == exchListPtr % endPointID) then
               commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldInPtr % dimSizes(1)
               comm_list_found = .true.
               exit
             end if
  
             commListPtr => commListPtr % next
           end do
  
           ! If no comm list exists, create a new one.
           if(.not. comm_list_found) then
             if(.not.associated(sendList)) then
               allocate(sendList)
               nullify(sendList % next)
               commListPtr => sendList
             else
               commListPtr => sendList
               commListPtr2 => commListPtr % next
               do while(associated(commListPtr2))
                 commListPtr => commListPtr % next
                 commListPtr2 => commListPtr % next
               end do
    
               allocate(commListPtr % next)
               commListPtr => commListPtr % next
               nullify(commListPtr % next)
             end if
             commListPtr % procID = exchListPtr % endPointID
             commListPtr % nList = exchListPtr % nList * fieldInPtr % dimSizes(1)
           end if
  
           exchListPtr => exchListPtr % next
         end do
  
         fieldInPtr => fieldInPtr % next
       end do
     end do

     ! Allocate sendLists, copy data into buffer, and initiate mpi_isends
     commListPtr => sendList
     do while(associated(commListPtr))
       allocate(commListPtr % rbuffer(commListPtr % nList))
       nullify(commListPtr % ibuffer)
       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               do i = 1, exchListPtr % nList
                 do j = 1, fieldInPtr % dimSizes(1)
                   iBuffer = (exchListPtr % destList(i)-1) * fieldInPtr % dimSizes(1) + j + bufferOffset
                   commListPtr % rbuffer(iBuffer) = fieldInPtr % array(j, exchListPtr % srcList(i))
                   nAdded = nAdded + 1
                 end do
               end do
             end if
  
             exchListPtr => exchListPtr % next
           end do
  
           fieldInPtr => fieldInPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do

       call MPI_Isend(commListPtr % rbuffer, commListPtr % nlist, MPI_realKIND, &
                      commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)

       commListPtr => commListPtr % next
     end do

#endif     

     ! Handle Local Copies. Only local copies if no MPI
     do iHalo = 1, nHaloLayers
       fieldInPtr => fieldIn
       do while(associated(fieldInPtr))
         exchListPtr => fieldInPtr % copyList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             if(exchListPtr % endPointID == fieldOutPtr % block % localBlockID) then
               do i = 1, exchListPtr % nList
                 fieldOutPtr % array(:, exchListPtr % destList(i)) = fieldInPtr % array(:, exchListPtr % srcList(i))
               end do
             end if
             fieldOutPtr => fieldOutPtr % next
           end do
  
           exchListPtr => exchListPtr % next
         end do
         fieldInPtr => fieldInPtr % next
       end do
     end do

#ifdef _MPI
     ! Wait for MPI_Irecv's to finish, and unpack data.
     commListPtr => recvList
     do while(associated(commListPtr))
       call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               do i = 1, exchListPtr % nList
                 do j = 1, fieldOutPtr % dimSizes(1)
                   iBuffer = (exchListPtr % srcList(i)-1) * fieldOutPtr % dimSizes(1) + j + bufferOffset
                   fieldOutPtr % array(j, exchListPtr % destList(i)) = commListPtr % rbuffer(iBuffer)
                 end do
               end do
               nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1))
             end if
             exchListPtr => exchListPtr % next
           end do
  
           fieldOutPtr => fieldOutPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do

       commListPtr => commListPtr % next
     end do

     ! Wait for MPI_Isend's to finish.
     commListPtr => sendList
     do while(associated(commListPtr))
       call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
       commListPtr => commListPtr % next
     end do

     ! Destroy commLists.
     call mpas_dmpar_destroy_communication_list(sendList)
     call mpas_dmpar_destroy_communication_list(recvList)
#endif

     deallocate(haloLayers)

   end subroutine mpas_dmpar_alltoall_field2d_real!}}}

   subroutine mpas_dmpar_alltoall_field3d_real(fieldIn, fieldout, haloLayersIn)!{{{

     implicit none

     type (field3dReal), pointer :: fieldIn
     type (field3dReal), pointer :: fieldOut
     integer, dimension(:), pointer, optional :: haloLayersIn

     type (field3dReal), pointer :: fieldInPtr, fieldOutPtr
     type (mpas_exchange_list), pointer :: exchListPtr, exchListPtr2
     type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
     type (dm_info), pointer :: dminfo

     logical :: comm_list_found

     integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
     integer :: nAdded, bufferOffset
     integer :: mpi_ierr
     integer :: iHalo, iBuffer, i, j, k
     integer :: nHaloLayers
     integer, dimension(:), pointer :: haloLayers

     dminfo => fieldIn % block % domain % dminfo

     if(present(haloLayersIn)) then
       nHaloLayers = size(haloLayersIn)
       allocate(haloLayers(nHaloLayers))
       do iHalo = 1, nHaloLayers
         haloLayers(iHalo) = haloLayersIn(iHalo)
       end do
     else
       nHaloLayers = size(fieldIn % sendList % halos)
       allocate(haloLayers(nHaloLayers))
       do iHalo = 1, nHaloLayers
         haloLayers(iHalo) = iHalo
       end do
     end if

#ifdef _MPI
     nullify(sendList)
     nullify(recvList)

     ! Setup recieve lists.
     do iHalo = 1, nHaloLayers
       fieldOutPtr => fieldOut
       do while(associated(fieldOutPtr))
         exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           comm_list_found = .false.
  
           ! Search for an already created commList to this processor.
           commListPtr => recvList
           do while(associated(commListPtr))
             if(commListPtr % procID == exchListPtr % endPointID) then
               commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2)
               comm_list_found = .true.
               exit
             end if
  
             commListPtr => commListPtr % next
           end do
  
           ! If no comm list exists, create a new one.
           if(.not. comm_list_found) then
             if(.not.associated(recvList)) then
               allocate(recvList)
               nullify(recvList % next)
               commListPtr => recvList
             else
               commListPtr => recvList
               commListPtr2 => commListPtr % next
               do while(associated(commListPtr2))
                 commListPtr => commListPtr % next
                 commListPtr2 => commListPtr % next
               end do

               allocate(commListPtr % next)
               commListPtr => commListPtr % next
               nullify(commListPtr % next)
             end if
  
             commListPtr % procID = exchListPtr % endPointID
             commListPtr % nList = exchListPtr % nList * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2)
           end if

           exchListPtr => exchListPtr % next
         end do
  
         fieldOutPtr => fieldOutPtr % next
       end do
     end do

     ! Determine size of receive list buffers.
     commListPtr => recvList
     do while(associated(commListPtr))
       call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2))
             end if
             exchListPtr => exchListPtr % next
           end do
  
           fieldOutPtr => fieldOutPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do
       commListPtr % nList = nAdded

       commListPtr => commListPtr % next
     end do

     ! Allocate buffers for recieves, and initiate mpi_irecv calls.
     commListPtr => recvList
     do while(associated(commListPtr))
       allocate(commListPtr % rbuffer(commListPtr % nList))
       nullify(commListPtr % ibuffer)
       call MPI_Irecv(commListPtr % rbuffer, commListPtr % nList, MPI_realKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)
       commListPtr => commListPtr % next
     end do

     ! Setup send lists, and determine the size of their buffers.
     do iHalo = 1, nHaloLayers
       fieldInPtr => fieldIn
       do while(associated(fieldInPtr))
         exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           comm_list_found = .false.
  
           ! Search for an already created commList to this processor.
           commListPtr => sendList
           do while(associated(commListPtr))
             if(commListPtr % procID == exchListPtr % endPointID) then
               commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2)
               comm_list_found = .true.
               exit
             end if
  
             commListPtr => commListPtr % next
           end do
  
           ! If no comm list exists, create a new one.
           if(.not. comm_list_found) then
             if(.not.associated(sendList)) then
               allocate(sendList)
               nullify(sendList % next)
               commListPtr => sendList
             else
               commListPtr => sendList
               commListPtr2 => commListPtr % next
               do while(associated(commListPtr2))
                 commListPtr => commListPtr % next
                 commListPtr2 => commListPtr % next
               end do
    
               allocate(commListPtr % next)
               commListPtr => commListPtr % next
               nullify(commListPtr % next)
             end if
             commListPtr % procID = exchListPtr % endPointID
             commListPtr % nList = exchListPtr % nList * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2)
           end if
  
           exchListPtr => exchListPtr % next
         end do
  
         fieldInPtr => fieldInPtr % next
       end do
     end do

     ! Allocate sendLists, copy data into buffer, and initiate mpi_isends
     commListPtr => sendList
     do while(associated(commListPtr))
       allocate(commListPtr % rbuffer(commListPtr % nList))
       nullify(commListPtr % ibuffer)
       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldInPtr => fieldIn
         do while(associated(fieldInPtr))
           exchListPtr => fieldInPtr % sendList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               do i = 1, exchListPtr % nList
                 do j = 1, fieldInPtr % dimSizes(2)
                   do k = 1, fieldInPtr % dimSizes(1)
                     iBuffer = (exchListPtr % destList(i)-1) * fieldInPtr % dimSizes(1) * fieldInPtr % dimSizes(2) + (j-1) * fieldInPtr % dimSizes(1) + k + bufferOffset
                     commListPtr % rbuffer(iBuffer) = fieldInPtr % array(k, j, exchListPtr % srcList(i))
                     nAdded = nAdded + 1
                   end do
                 end do
               end do
             end if
  
             exchListPtr => exchListPtr % next
           end do
  
           fieldInPtr => fieldInPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do

       call MPI_Isend(commListPtr % rbuffer, commListPtr % nlist, MPI_realKIND, &
                      commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)

       commListPtr => commListPtr % next
     end do

#endif     

     ! Handle Local Copies. Only local copies if no MPI
     do iHalo = 1, nHaloLayers
       fieldInPtr => fieldIn
       do while(associated(fieldInPtr))
         exchListPtr => fieldInPtr % copyList % halos(haloLayers(iHalo)) % exchList
         do while(associated(exchListPtr))
           fieldOutPtr => fieldOut
           do while(associated(fieldOutPtr))
             if(exchListPtr % endPointID == fieldOutPtr % block % localBlockID) then
               do i = 1, exchListPtr % nList
                 fieldOutPtr % array(:, :, exchListPtr % destList(i)) = fieldInPtr % array(:, :, exchListPtr % srcList(i))
               end do
             end if
             fieldOutPtr => fieldOutPtr % next
           end do
  
           exchListPtr => exchListPtr % next
         end do
         fieldInPtr => fieldInPtr % next
       end do
     end do

#ifdef _MPI
     ! Wait for MPI_Irecv's to finish, and unpack data.
     commListPtr => recvList
     do while(associated(commListPtr))
       call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)

       bufferOffset = 0
       do iHalo = 1, nHaloLayers
         nAdded = 0
         fieldOutPtr => fieldOut
         do while(associated(fieldOutPtr))
           exchListPtr => fieldOutPtr % recvList % halos(haloLayers(iHalo)) % exchList
           do while(associated(exchListPtr))
             if(exchListPtr % endPointID == commListPtr % procID) then
               do i = 1, exchListPtr % nList
                 do j = 1, fieldOutPtr % dimSizes(2)
                   do k = 1, fieldOutPtr % dimSizes(1)
                     iBuffer = (exchListPtr % srcList(i)-1) * fieldOutPtr % dimSizes(2) * fieldOutPtr % dimSizes(1) + (j-1) * fieldOutPtr % dimSizes(1) + k + bufferOffset
                     fieldOutPtr % array(k, j, exchListPtr % destList(i)) = commListPtr % rbuffer(iBuffer)
                   end do
                 end do
               end do
               nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldOutPtr % dimSizes(1) * fieldOutPtr % dimSizes(2))
             end if
             exchListPtr => exchListPtr % next
           end do
  
           fieldOutPtr => fieldOutPtr % next
         end do
         bufferOffset = bufferOffset + nAdded
       end do

       commListPtr => commListPtr % next
     end do

     ! Wait for MPI_Isend's to finish.
     commListPtr => sendList
     do while(associated(commListPtr))
       call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
       commListPtr => commListPtr % next
     end do

     ! Destroy commLists.
     call mpas_dmpar_destroy_communication_list(sendList)
     call mpas_dmpar_destroy_communication_list(recvList)
#endif

     deallocate(haloLayers)

   end subroutine mpas_dmpar_alltoall_field3d_real!}}}


   subroutine mpas_dmpar_exch_halo_field1d_integer(field, haloLayersIn)!{{{

      implicit none

      type (field1DInteger), pointer :: field
      integer, dimension(:), intent(in), optional :: haloLayersIn

      type (dm_info), pointer :: dminfo
      type (field1DInteger), pointer :: fieldCursor, fieldCursor2
      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
      integer :: mpi_ierr
      integer :: nHaloLayers, iHalo, i
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      logical :: comm_list_found

      do i = 1, 1
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo

      if(present(haloLayersIn)) then
        nHaloLayers = size(haloLayersIn)
        allocate(haloLayers(nHaloLayers))
        do iHalo = 1, nHaloLayers
          haloLayers(iHalo) = haloLayersIn(iHalo)
        end do
      else
        nHaloLayers = size(field % sendList % halos)
        allocate(haloLayers(nHaloLayers))
        do iHalo = 1, nHaloLayers
          haloLayers(iHalo) = iHalo
        end do
      end if

#ifdef _MPI
      ! Allocate communication lists, and setup dead header node.
      allocate(sendList)
      nullify(sendList % next)
      sendList % procID = -1
      sendList % nList = 0

      allocate(recvList)
      nullify(recvList % next)
      recvList % procID = -1
      recvList % nList = 0

      dminfo   = field % block % domain % dminfo

      ! Determine size of buffers for communication lists
      fieldCursor => field
      do while(associated(fieldCursor))

        ! Need to aggregate across halo layers
        do iHalo = 1, nHaloLayers
          
          ! Determine size from send lists
          exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
          do while(associated(exchListPtr))
            comm_list_found = .false.

            commListPtr => sendList
            do while(associated(commListPtr))
              if(commListPtr % procID == exchListPtr % endPointId) then
                comm_list_found = .true.
                commListPtr % nList = commListPtr % nList + exchListPtr % nList
                exit
              end if

              commListPtr => commListPtr % next
            end do

            if(.not. comm_list_found) then
              commListPtr => sendList
              commListPtr2 => commListPtr % next
              do while(associated(commListPtr2))
                commListPtr => commListPtr % next
                commListPtr2 => commListPtr % next
              end do

              allocate(commListPtr % next)
              commListPtr => commListPtr % next
              nullify(commListPtr % next)
              commListPtr % procID = exchListPtr % endPointID
              commListPtr % nList = exchListPtr % nList
            end if

            exchListPtr => exchListPtr % next
          end do

          ! Setup recv lists
          exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
          do while(associated(exchListPtr))
            comm_list_found = .false.

            commListPtr => recvList
            do while(associated(commListPtr))
              if(commListPtr % procID == exchListPtr % endPointId) then
                comm_list_found = .true.
                exit
              end if

              commListPtr => commListPtr % next
            end do

            if(.not. comm_list_found) then
              commListPtr => recvList
              commListPtr2 => commListPtr % next
              do while(associated(commListPtr2))
                commListPtr => commListPtr % next
                commListPtr2 => commListPtr % next
              end do

              allocate(commListPtr % next)
              commListPtr => commListPtr % next
              nullify(commListPtr % next)
              commListPtr % procID = exchListPtr % endPointID
            end if

            exchListPtr => exchListPtr % next
          end do
        end do

        fieldCursor => fieldCursor % next
      end do

      ! Remove the dead head pointer on send and recv list
      commListPtr => sendList
      sendList => sendList % next
      deallocate(commListPtr)

      commListPtr => recvList
      recvList => recvList % next
      deallocate(commListPtr)

      ! Determine size of recieve lists
      commListPtr => recvList
      do while(associated(commListPtr))
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                nAdded = max(nAdded, maxval(exchListPtr % srcList))
              end if
              exchListPtr => exchListPtr % next
            end do
            
            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do
        commListPtr % nList = bufferOffset

        commListPtr => commListPtr % next
      end do


      ! Allocate space in recv lists, and initiate mpi_irecv calls
      commListPtr => recvList
      do while(associated(commListPtr))
        allocate(commListPtr % ibuffer(commListPtr % nList))
        nullify(commListPtr % rbuffer)
        call MPI_Irecv(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)

        commListPtr => commListPtr % next
      end do

      ! Allocate space in send lists, copy data into buffer, and initiate mpi_isend calls
      commListPtr => sendList
      do while(associated(commListPtr))
        allocate(commListPtr % ibuffer(commListPtr % nList))
        nullify(commListPtr % rbuffer)
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                do  i = 1, exchListPtr % nList
                  commListPtr % ibuffer(exchListPtr % destList(i) + bufferOffset) = fieldCursor % array(exchListPtr % srcList(i))
                  nAdded = nAdded + 1

                end do
              end if

              exchListPtr => exchListPtr % next
            end do

            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do

        call MPI_Isend(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)

        commListPtr => commListPtr % next
      end do
#endif

      ! Handle local copy. If MPI is off, then only local copies are performed.
      fieldCursor => field
      do while(associated(fieldCursor))
        do iHalo = 1, nHaloLayers
          exchListPtr => fieldCursor % copyList % halos(haloLayers(iHalo)) % exchList

          do while(associated(exchListPtr))
            fieldCursor2 => field
            do while(associated(fieldCursor2))
              if(exchListPtr % endPointID == fieldCursor2 % block % localBlockID) then
                do i = 1, exchListPtr % nList
                  fieldCursor2 % array(exchListPtr % destList(i)) = fieldCursor % array(exchListPtr % srcList(i))
                end do
              end if
              
              fieldCursor2 => fieldCursor2 % next
            end do

            exchListPtr => exchListPtr % next
          end do
        end do

        fieldCursor => fieldCursor % next
      end do

#ifdef _MPI

      ! Wait for mpi_irecv to finish, and unpack data from buffer
      commListPtr => recvList
      do while(associated(commListPtr))
        call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                do i = 1, exchListPtr % nList
                  fieldCursor % array(exchListPtr % destList(i)) = commListPtr % ibuffer(exchListPtr % srcList(i) + bufferOffset)
                end do
                nAdded = max(nAdded, maxval(exchListPtr % srcList))
              end if
              exchListPtr => exchListPtr % next
            end do
            
            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do
        commListPtr => commListPtr % next
      end do

      ! wait for mpi_isend to finish.
      commListPtr => sendList
      do while(associated(commListPtr))
        call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
        commListPtr => commListPtr % next
      end do

     ! Destroy commLists.
     call mpas_dmpar_destroy_communication_list(sendList)
     call mpas_dmpar_destroy_communication_list(recvList)
#endif

     deallocate(haloLayers)

   end subroutine mpas_dmpar_exch_halo_field1d_integer!}}}

   subroutine mpas_dmpar_exch_halo_field2d_integer(field, haloLayersIn)!{{{

      implicit none

      type (field2DInteger), pointer :: field
      integer, dimension(:), intent(in), optional :: haloLayersIn

      type (dm_info), pointer :: dminfo
      type (field2DInteger), pointer :: fieldCursor, fieldCursor2
      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
      integer :: mpi_ierr
      integer :: nHaloLayers, iHalo, i, j
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      logical :: comm_list_found

      do i = 1, 2
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo

      if(present(haloLayersIn)) then
        nHaloLayers = size(haloLayersIn)
        allocate(haloLayers(nHaloLayers))
        do iHalo = 1, nHaloLayers
          haloLayers(iHalo) = haloLayersIn(iHalo)
        end do
      else
        nHaloLayers = size(field % sendList % halos)
        allocate(haloLayers(nHaloLayers))
        do iHalo = 1, nHaloLayers
          haloLayers(iHalo) = iHalo
        end do
      end if

#ifdef _MPI
      ! Allocate communication lists, and setup dead header nodes
      allocate(sendList)
      nullify(sendList % next)
      sendList % procID = -1
      sendList % nList = 0

      allocate(recvList)
      nullify(recvList % next)
      recvList % procID = -1
      recvList % nList = 0

      dminfo   = field % block % domain % dminfo

      ! Determine size of buffers for communication lists
      fieldCursor => field
      do while(associated(fieldCursor))

        ! Need to aggregate across halo layers
        do iHalo = 1, nHaloLayers
          
          ! Determine size from send lists
          exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
          do while(associated(exchListPtr))
            comm_list_found = .false.

            commListPtr => sendList
            do while(associated(commListPtr))
              if(commListPtr % procID == exchListPtr % endPointId) then
                comm_list_found = .true.
                commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldCursor % dimSizes(1)
                exit
              end if

              commListPtr => commListPtr % next
            end do

            if(.not. comm_list_found) then
              commListPtr => sendList
              commListPtr2 => commListPtr % next
              do while(associated(commListPtr2))
                commListPtr => commListPtr % next
                commListPtr2 => commListPtr % next
              end do

              allocate(commListPtr % next)
              commListPtr => commListPtr % next
              nullify(commListPtr % next)
              commListPtr % procID = exchListPtr % endPointID
              commListPtr % nList = exchListPtr % nList * fieldCursor % dimSizes(1)
            end if

            exchListPtr => exchListPtr % next
          end do

          ! Setup recv lists
          exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
          do while(associated(exchListPtr))
            comm_list_found = .false.

            commListPtr => recvList
            do while(associated(commListPtr))
              if(commListPtr % procID == exchListPtr % endPointId) then
                comm_list_found = .true.
                exit
              end if

              commListPtr => commListPtr % next
            end do

            if(.not. comm_list_found) then
              commListPtr => recvList
              commListPtr2 => commListPtr % next
              do while(associated(commListPtr2))
                commListPtr => commListPtr % next
                commListPtr2 => commListPtr % next
              end do

              allocate(commListPtr % next)
              commListPtr => commListPtr % next
              nullify(commListPtr % next)
              commListPtr % procID = exchListPtr % endPointID
            end if

            exchListPtr => exchListPtr % next
          end do
        end do

        fieldCursor => fieldCursor % next
      end do

      ! Remove the dead head pointer on send and recv list
      commListPtr => sendList
      sendList => sendList % next
      deallocate(commListPtr)

      commListPtr => recvList
      recvList => recvList % next
      deallocate(commListPtr)

      ! Determine size of recieve list buffers
      commListPtr => recvList
      do while(associated(commListPtr))
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1))
              end if
              exchListPtr => exchListPtr % next
            end do
            
            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do
        commListPtr % nList = bufferOffset

        commListPtr => commListPtr % next
      end do

      ! Allocate space in recv lists, and initiate mpi_irecv calls
      commListPtr => recvList
      do while(associated(commListPtr))
        allocate(commListPtr % ibuffer(commListPtr % nList))
        nullify(commListPtr % rbuffer)
        call MPI_Irecv(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)

        commListPtr => commListPtr % next
      end do

      ! Allocate space in send lists, copy data into buffer, and initiate mpi_isend calls
      commListPtr => sendList
      do while(associated(commListPtr))
        allocate(commListPtr % ibuffer(commListPtr % nList))
        nullify(commListPtr % rbuffer)
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                do  i = 1, exchListPtr % nList
                  do j = 1, fieldCursor % dimSizes(1)
                    commListPtr % ibuffer((exchListPtr % destList(i)-1) * fieldCursor % dimSizes(1) + j + bufferOffset) = fieldCursor % array(j, exchListPtr % srcList(i))
                    nAdded = nAdded + 1
                  end do
                end do
              end if

              exchListPtr => exchListPtr % next
            end do

            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do

        call MPI_Isend(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)
        commListPtr => commListPtr % next
      end do
#endif

      ! Handle local copy. If MPI is off, then only local copies are performed.
      fieldCursor => field
      do while(associated(fieldCursor))
        do iHalo = 1, nHaloLayers
          exchListPtr => fieldCursor % copyList % halos(haloLayers(iHalo)) % exchList

          do while(associated(exchListPtr))
            fieldCursor2 => field
            do while(associated(fieldCursor2))
              if(exchListPtr % endPointID == fieldCursor2 % block % localBlockID) then
                do i = 1, exchListPtr % nList
                  fieldCursor2 % array(:, exchListPtr % destList(i)) = fieldCursor % array(:, exchListPtr % srcList(i))
                end do
              end if
              
              fieldCursor2 => fieldCursor2 % next
            end do

            exchListPtr => exchListPtr % next
          end do
        end do

        fieldCursor => fieldCursor % next
      end do

#ifdef _MPI

      ! Wait for mpi_irecv to finish, and unpack data from buffer
      commListPtr => recvList
      do while(associated(commListPtr))
        call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                do i = 1, exchListPtr % nList
                  do j = 1, fieldCursor % dimSizes(1)
                    fieldCursor % array(j, exchListPtr % destList(i)) = commListPtr % ibuffer((exchListPtr % srcList(i)-1)*fieldCursor % dimSizes(1) + j + bufferOffset)
                  end do
                end do
                nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1))
              end if
              exchListPtr => exchListPtr % next
            end do
            
            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do
        commListPtr => commListPtr % next
      end do

      ! wait for mpi_isend to finish.
      commListPtr => sendList
      do while(associated(commListPtr))
        call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
        commListPtr => commListPtr % next
      end do

     ! Destroy commLists.
     call mpas_dmpar_destroy_communication_list(sendList)
     call mpas_dmpar_destroy_communication_list(recvList)
#endif

     deallocate(haloLayers)

   end subroutine mpas_dmpar_exch_halo_field2d_integer!}}}

   subroutine mpas_dmpar_exch_halo_field3d_integer(field, haloLayersIn)!{{{

      implicit none

      type (field3DInteger), pointer :: field
      integer, dimension(:), intent(in), optional :: haloLayersIn

      type (dm_info), pointer :: dminfo
      type (field3DInteger), pointer :: fieldCursor, fieldCursor2
      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
      integer :: mpi_ierr
      integer :: nHaloLayers, iHalo, i, j, k
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      logical :: comm_list_found

      do i = 1, 3
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo

      if(present(haloLayersIn)) then
        nHaloLayers = size(haloLayersIn)
        allocate(haloLayers(nHaloLayers))
        do iHalo = 1, nHaloLayers
          haloLayers(iHalo) = haloLayersIn(iHalo)
        end do
      else
        nHaloLayers = size(field % sendList % halos)
        allocate(haloLayers(nHaloLayers))
        do iHalo = 1, nHaloLayers
          haloLayers(iHalo) = iHalo
        end do
      end if

#ifdef _MPI
      ! Allocate communication lists, and setup dead header nodes
      allocate(sendList)
      nullify(sendList % next)
      sendList % procID = -1
      sendList % nList = 0

      allocate(recvList)
      nullify(recvList % next)
      recvList % procID = -1
      recvList % nList = 0

      dminfo   = field % block % domain % dminfo

      ! Determine size of buffers for communication lists
      fieldCursor => field
      do while(associated(fieldCursor))

        ! Need to aggregate across halo layers
        do iHalo = 1, nHaloLayers
          
          ! Determine size from send lists
          exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
          do while(associated(exchListPtr))
            comm_list_found = .false.

            commListPtr => sendList
            do while(associated(commListPtr))
              if(commListPtr % procID == exchListPtr % endPointId) then
                comm_list_found = .true.
                commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2)
                exit
              end if

              commListPtr => commListPtr % next
            end do

            if(.not. comm_list_found) then
              commListPtr => sendList
              commListPtr2 => commListPtr % next
              do while(associated(commListPtr2))
                commListPtr => commListPtr % next
                commListPtr2 => commListPtr % next
              end do

              allocate(commListPtr % next)
              commListPtr => commListPtr % next
              nullify(commListPtr % next)
              commListPtr % procID = exchListPtr % endPointID
              commListPtr % nList = exchListPtr % nList * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2)
            end if

            exchListPtr => exchListPtr % next
          end do

          ! Setup recv lists
          exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
          do while(associated(exchListPtr))
            comm_list_found = .false.

            commListPtr => recvList
            do while(associated(commListPtr))
              if(commListPtr % procID == exchListPtr % endPointId) then
                comm_list_found = .true.
                exit
              end if

              commListPtr => commListPtr % next
            end do

            if(.not. comm_list_found) then
              commListPtr => recvList
              commListPtr2 => commListPtr % next
              do while(associated(commListPtr2))
                commListPtr => commListPtr % next
                commListPtr2 => commListPtr % next
              end do

              allocate(commListPtr % next)
              commListPtr => commListPtr % next
              nullify(commListPtr % next)
              commListPtr % procID = exchListPtr % endPointID
            end if

            exchListPtr => exchListPtr % next
          end do
        end do

        fieldCursor => fieldCursor % next
      end do

      ! Remove the dead head pointer on send and recv list
      commListPtr => sendList
      sendList => sendList % next
      deallocate(commListPtr)

      commListPtr => recvList
      recvList => recvList % next
      deallocate(commListPtr)

      ! Determine size of recv lists
      commListPtr => recvList
      do while(associated(commListPtr))
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2))
              end if
              exchListPtr => exchListPtr % next
            end do
            
            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do
        commListPtr % nList = bufferOffset

        commListPtr => commListPtr % next
      end do

      ! Allocate space in recv lists, and initiate mpi_irecv calls
      commListPtr => recvList
      do while(associated(commListPtr))
        allocate(commListPtr % ibuffer(commListPtr % nList))
        nullify(commListPtr % rbuffer)
        call MPI_Irecv(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)

        commListPtr => commListPtr % next
      end do

      ! Allocate space in send lists, copy data into buffer, and initiate mpi_isend calls
      commListPtr => sendList
      do while(associated(commListPtr))
        allocate(commListPtr % ibuffer(commListPtr % nList))
        nullify(commListPtr % rbuffer)
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                do  i = 1, exchListPtr % nList
                  do j = 1, fieldCursor % dimSizes(2)
                    do k = 1, fieldCursor % dimSizes(1)
                      commListPtr % ibuffer((exchListPtr % destList(i)-1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                          + (j-1) * fieldCursor % dimSizes(1) + k  + bufferOffset) = fieldCursor % array(k, j, exchListPtr % srcList(i))
                      nAdded = nAdded + 1
                    end do
                  end do
                end do
              end if

              exchListPtr => exchListPtr % next
            end do

            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do

        call MPI_Isend(commListPtr % ibuffer, commListPtr % nList, MPI_INTEGERKIND, commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)
        commListPtr => commListPtr % next
      end do
#endif

      ! Handle local copy. If MPI is off, then only local copies are performed.
      fieldCursor => field
      do while(associated(fieldCursor))
        do iHalo = 1, nHaloLayers
          exchListPtr => fieldCursor % copyList % halos(haloLayers(iHalo)) % exchList

          do while(associated(exchListPtr))
            fieldCursor2 => field
            do while(associated(fieldCursor2))
              if(exchListPtr % endPointID == fieldCursor2 % block % localBlockID) then
                do i = 1, exchListPtr % nList
                  fieldCursor2 % array(:, :, exchListPtr % destList(i)) = fieldCursor % array(:, :, exchListPtr % srcList(i))
                end do
              end if
              
              fieldCursor2 => fieldCursor2 % next
            end do

            exchListPtr => exchListPtr % next
          end do
        end do

        fieldCursor => fieldCursor % next
      end do

#ifdef _MPI

      ! Wait for mpi_irecv to finish, and unpack data from buffer
      commListPtr => recvList
      do while(associated(commListPtr))
        call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                do i = 1, exchListPtr % nList
                  do j = 1, fieldCursor % dimSizes(2)
                    do k = 1, fieldCursor % dimSizes(1)
                      fieldCursor % array(k, j, exchListPtr % destList(i)) = commListPtr % ibuffer((exchListPtr % srcList(i)-1)*fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                                                           + (j-1)*fieldCursor % dimSizes(1) + k + bufferOffset)
                    end do
                  end do
                end do
                nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2))
              end if
              exchListPtr => exchListPtr % next
            end do
            
            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do
        commListPtr => commListPtr % next
      end do

      ! wait for mpi_isend to finish.
      commListPtr => sendList
      do while(associated(commListPtr))
        call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
        commListPtr => commListPtr % next
      end do

     ! Destroy commLists.
     call mpas_dmpar_destroy_communication_list(sendList)
     call mpas_dmpar_destroy_communication_list(recvList)
#endif

     deallocate(haloLayers)

   end subroutine mpas_dmpar_exch_halo_field3d_integer!}}}

   subroutine mpas_dmpar_exch_halo_field1d_real(field, haloLayersIn)!{{{

      implicit none

      type (field1dReal), pointer :: field
      integer, dimension(:), intent(in), optional :: haloLayersIn

      type (dm_info), pointer :: dminfo
      type (field1dReal), pointer :: fieldCursor, fieldCursor2
      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
      integer :: mpi_ierr
      integer :: nHaloLayers, iHalo, i
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      logical :: comm_list_found

      do i = 1, 1
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo

      if(present(haloLayersIn)) then
        nHaloLayers = size(haloLayersIn)
        allocate(haloLayers(nHaloLayers))
        do iHalo = 1, nHaloLayers
          haloLayers(iHalo) = haloLayersIn(iHalo)
        end do
      else
        nHaloLayers = size(field % sendList % halos)
        allocate(haloLayers(nHaloLayers))
        do iHalo = 1, nHaloLayers
          haloLayers(iHalo) = iHalo
        end do
      end if

#ifdef _MPI
      ! Allocate communication lists, and setup dead header nodes
      allocate(sendList)
      nullify(sendList % next)
      sendList % procID = -1
      sendList % nList = 0

      allocate(recvList)
      nullify(recvList % next)
      recvList % procID = -1
      recvList % nList = 0

      dminfo   = field % block % domain % dminfo

      ! Determine size of buffers for communication lists
      fieldCursor => field
      do while(associated(fieldCursor))

        ! Need to aggregate across halo layers
        do iHalo = 1, nHaloLayers
          
          ! Determine size from send lists
          exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
          do while(associated(exchListPtr))
            comm_list_found = .false.

            commListPtr => sendList
            do while(associated(commListPtr))
              if(commListPtr % procID == exchListPtr % endPointId) then
                comm_list_found = .true.
                commListPtr % nList = commListPtr % nList + exchListPtr % nList
                exit
              end if

              commListPtr => commListPtr % next
            end do

            if(.not. comm_list_found) then
              commListPtr => sendList
              commListPtr2 => commListPtr % next
              do while(associated(commListPtr2))
                commListPtr => commListPtr % next
                commListPtr2 => commListPtr % next
              end do

              allocate(commListPtr % next)
              commListPtr => commListPtr % next
              nullify(commListPtr % next)
              commListPtr % procID = exchListPtr % endPointID
              commListPtr % nList = exchListPtr % nList
            end if

            exchListPtr => exchListPtr % next
          end do

          ! Setup recv lists
          exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
          do while(associated(exchListPtr))
            comm_list_found = .false.

            commListPtr => recvList
            do while(associated(commListPtr))
              if(commListPtr % procID == exchListPtr % endPointId) then
                comm_list_found = .true.
                exit
              end if

              commListPtr => commListPtr % next
            end do

            if(.not. comm_list_found) then
              commListPtr => recvList
              commListPtr2 => commListPtr % next
              do while(associated(commListPtr2))
                commListPtr => commListPtr % next
                commListPtr2 => commListPtr % next
              end do

              allocate(commListPtr % next)
              commListPtr => commListPtr % next
              nullify(commListPtr % next)
              commListPtr % procID = exchListPtr % endPointID
            end if

            exchListPtr => exchListPtr % next
          end do
        end do

        fieldCursor => fieldCursor % next
      end do

      ! Remove the dead head pointer on send and recv list
      commListPtr => sendList
      sendList => sendList % next
      deallocate(commListPtr)

      commListPtr => recvList
      recvList => recvList % next
      deallocate(commListPtr)

      ! Determine size of recv lists
      commListPtr => recvList
      do while(associated(commListPtr))
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                nAdded = max(nAdded, maxval(exchListPtr % srcList))
              end if
              exchListPtr => exchListPtr % next
            end do
            
            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do
        commListPtr % nList = bufferOffset

        commListPtr => commListPtr % next
      end do


      ! Allocate space in recv lists, and initiate mpi_irecv calls
      commListPtr => recvList
      do while(associated(commListPtr))
        allocate(commListPtr % rbuffer(commListPtr % nList))
        nullify(commListPtr % ibuffer)
        call MPI_Irecv(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)

        commListPtr => commListPtr % next
      end do

      ! Allocate space in send lists, copy data into buffer, and initiate mpi_isend calls
      commListPtr => sendList
      do while(associated(commListPtr))
        allocate(commListPtr % rbuffer(commListPtr % nList))
        nullify(commListPtr % ibuffer)
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                do  i = 1, exchListPtr % nList
                  commListPtr % rbuffer(exchListPtr % destList(i) + bufferOffset) = fieldCursor % array(exchListPtr % srcList(i))
                  nAdded = nAdded + 1
                end do
              end if

              exchListPtr => exchListPtr % next
            end do

            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do

        call MPI_Isend(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)
        commListPtr => commListPtr % next
      end do
#endif

      ! Handle local copy. If MPI is off, then only local copies are performed.
      fieldCursor => field
      do while(associated(fieldCursor))
        do iHalo = 1, nHaloLayers
          exchListPtr => fieldCursor % copyList % halos(haloLayers(iHalo)) % exchList

          do while(associated(exchListPtr))
            fieldCursor2 => field
            do while(associated(fieldCursor2))
              if(exchListPtr % endPointID == fieldCursor2 % block % localBlockID) then
                do i = 1, exchListPtr % nList
                  fieldCursor2 % array(exchListPtr % destList(i)) = fieldCursor % array(exchListPtr % srcList(i))
                end do
              end if
              
              fieldCursor2 => fieldCursor2 % next
            end do

            exchListPtr => exchListPtr % next
          end do
        end do

        fieldCursor => fieldCursor % next
      end do

#ifdef _MPI

      ! Wait for mpi_irecv to finish, and unpack data from buffer
      commListPtr => recvList
      do while(associated(commListPtr))
        call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                do i = 1, exchListPtr % nList
                  fieldCursor % array(exchListPtr % destList(i)) = commListPtr % rbuffer(exchListPtr % srcList(i) + bufferOffset)
                end do
                nAdded = max(nAdded, maxval(exchListPtr % srcList))
              end if
              exchListPtr => exchListPtr % next
            end do
            
            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do
        commListPtr => commListPtr % next
      end do

      ! wait for mpi_isend to finish.
      commListPtr => sendList
      do while(associated(commListPtr))
        call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
        commListPtr => commListPtr % next
      end do

     ! Destroy commLists.
     call mpas_dmpar_destroy_communication_list(sendList)
     call mpas_dmpar_destroy_communication_list(recvList)
#endif

     deallocate(haloLayers)

   end subroutine mpas_dmpar_exch_halo_field1d_real!}}}

   subroutine mpas_dmpar_exch_halo_field2d_real(field, haloLayersIn)!{{{

      implicit none

      type (field2dReal), pointer :: field
      integer, dimension(:), intent(in), optional :: haloLayersIn

      type (dm_info), pointer :: dminfo
      type (field2dReal), pointer :: fieldCursor, fieldCursor2
      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
      integer :: mpi_ierr
      integer :: nHaloLayers, iHalo, i, j
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      logical :: comm_list_found

      do i = 1, 2
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo

      if(present(haloLayersIn)) then
        nHaloLayers = size(haloLayersIn)
        allocate(haloLayers(nHaloLayers))
        do iHalo = 1, nHaloLayers
          haloLayers(iHalo) = haloLayersIn(iHalo)
        end do
      else
        nHaloLayers = size(field % sendList % halos)
        allocate(haloLayers(nHaloLayers))
        do iHalo = 1, nHaloLayers
          haloLayers(iHalo) = iHalo
        end do
      end if

#ifdef _MPI
      ! Allocate communication lists, and setup dead header nodes
      allocate(sendList)
      nullify(sendList % next)
      sendList % procID = -1
      sendList % nList = 0

      allocate(recvList)
      nullify(recvList % next)
      recvList % procID = -1
      recvList % nList = 0

      dminfo   = field % block % domain % dminfo

      ! Determine size of buffers for communication lists
      fieldCursor => field
      do while(associated(fieldCursor))

        ! Need to aggregate across halo layers
        do iHalo = 1, nHaloLayers
          
          ! Determine size from send lists
          exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
          do while(associated(exchListPtr))
            comm_list_found = .false.

            commListPtr => sendList
            do while(associated(commListPtr))
              if(commListPtr % procID == exchListPtr % endPointId) then
                comm_list_found = .true.
                commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldCursor % dimSizes(1)
                exit
              end if

              commListPtr => commListPtr % next
            end do

            if(.not. comm_list_found) then
              commListPtr => sendList
              commListPtr2 => commListPtr % next
              do while(associated(commListPtr2))
                commListPtr => commListPtr % next
                commListPtr2 => commListPtr % next
              end do

              allocate(commListPtr % next)
              commListPtr => commListPtr % next
              nullify(commListPtr % next)
              commListPtr % procID = exchListPtr % endPointID
              commListPtr % nList = exchListPtr % nList * fieldCursor % dimSizes(1)
            end if

            exchListPtr => exchListPtr % next
          end do

          ! Setup recv lists
          exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
          do while(associated(exchListPtr))
            comm_list_found = .false.

            commListPtr => recvList
            do while(associated(commListPtr))
              if(commListPtr % procID == exchListPtr % endPointId) then
                comm_list_found = .true.
                commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldCursor % dimSizes(1)
                exit
              end if

              commListPtr => commListPtr % next
            end do

            if(.not. comm_list_found) then
              commListPtr => recvList
              commListPtr2 => commListPtr % next
              do while(associated(commListPtr2))
                commListPtr => commListPtr % next
                commListPtr2 => commListPtr % next
              end do

              allocate(commListPtr % next)
              commListPtr => commListPtr % next
              nullify(commListPtr % next)
              commListPtr % procID = exchListPtr % endPointID
              commListPtr % nList = exchListPtr % nList * fieldCursor % dimSizes(1)
            end if

            exchListPtr => exchListPtr % next
          end do
        end do

        fieldCursor => fieldCursor % next
      end do

      ! Remove the dead head pointer on send and recv list
      commListPtr => sendList
      sendList => sendList % next
      deallocate(commListPtr)

      commListPtr => recvList
      recvList => recvList % next
      deallocate(commListPtr)

      ! Determine size of recv lists
      commListPtr => recvList
      do while(associated(commListPtr))
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1))
              end if
              exchListPtr => exchListPtr % next
            end do
            
            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do
        commListPtr % nList = bufferOffset

        commListPtr => commListPtr % next
      end do

      ! Allocate space in recv lists, and initiate mpi_irecv calls
      commListPtr => recvList
      do while(associated(commListPtr))
        allocate(commListPtr % rbuffer(commListPtr % nList))
        nullify(commListPtr % ibuffer)
        call MPI_Irecv(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)

        commListPtr => commListPtr % next
      end do

      ! Allocate space in send lists, copy data into buffer, and initiate mpi_isend calls
      commListPtr => sendList
      do while(associated(commListPtr))
        allocate(commListPtr % rbuffer(commListPtr % nList))
        nullify(commListPtr % ibuffer)
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                do  i = 1, exchListPtr % nList
                  do j = 1, fieldCursor % dimSizes(1)
                    commListPtr % rbuffer((exchListPtr % destList(i)-1) * fieldCursor % dimSizes(1) + j + bufferOffset) = fieldCursor % array(j, exchListPtr % srcList(i))
                    nAdded = nAdded + 1
                  end do
                end do
              end if

              exchListPtr => exchListPtr % next
            end do

            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do

        call MPI_Isend(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)
        commListPtr => commListPtr % next
      end do
#endif

      ! Handle local copy. If MPI is off, then only local copies are performed.
      fieldCursor => field
      do while(associated(fieldCursor))
        do iHalo = 1, nHaloLayers
          exchListPtr => fieldCursor % copyList % halos(haloLayers(iHalo)) % exchList

          do while(associated(exchListPtr))
            fieldCursor2 => field
            do while(associated(fieldCursor2))
              if(exchListPtr % endPointID == fieldCursor2 % block % localBlockID) then
                do i = 1, exchListPtr % nList
                  fieldCursor2 % array(:, exchListPtr % destList(i)) = fieldCursor % array(:, exchListPtr % srcList(i))
                end do
              end if
              
              fieldCursor2 => fieldCursor2 % next
            end do

            exchListPtr => exchListPtr % next
          end do
        end do

        fieldCursor => fieldCursor % next
      end do

#ifdef _MPI

      ! Wait for mpi_irecv to finish, and unpack data from buffer
      commListPtr => recvList
      do while(associated(commListPtr))
        call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                do i = 1, exchListPtr % nList
                  do j = 1, fieldCursor % dimSizes(1)
                    fieldCursor % array(j, exchListPtr % destList(i)) = commListPtr % rbuffer((exchListPtr % srcList(i)-1)*fieldCursor % dimSizeS(1) + j + bufferOffset)
                  end do
                end do
                nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1))
              end if
              exchListPtr => exchListPtr % next
            end do
            
            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do
        commListPtr => commListPtr % next
      end do

      ! wait for mpi_isend to finish.
      commListPtr => sendList
      do while(associated(commListPtr))
        call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
        commListPtr => commListPtr % next
      end do

     ! Destroy commLists.
     call mpas_dmpar_destroy_communication_list(sendList)
     call mpas_dmpar_destroy_communication_list(recvList)
#endif

     deallocate(haloLayers)

   end subroutine mpas_dmpar_exch_halo_field2d_real!}}}

   subroutine mpas_dmpar_exch_halo_field3d_real(field, haloLayersIn)!{{{

      implicit none

      type (field3dReal), pointer :: field
      integer, dimension(:), intent(in), optional :: haloLayersIn

      type (dm_info), pointer :: dminfo
      type (field3dReal), pointer :: fieldCursor, fieldCursor2
      type (mpas_exchange_list), pointer :: exchListPtr
      type (mpas_communication_list), pointer :: sendList, recvList, commListPtr, commListPtr2
      integer :: mpi_ierr
      integer :: nHaloLayers, iHalo, i, j, k
      integer :: bufferOffset, nAdded
      integer, dimension(:), pointer :: haloLayers

      logical :: comm_list_found

      do i = 1, 3
        if(field % dimSizes(i) <= 0) then
          return
        end if
      end do

      dminfo => field % block % domain % dminfo

      if(present(haloLayersIn)) then
        nHaloLayers = size(haloLayersIn)
        allocate(haloLayers(nHaloLayers))
        do iHalo = 1, nHaloLayers
          haloLayers(iHalo) = haloLayersIn(iHalo)
        end do
      else
        nHaloLayers = size(field % sendList % halos)
        allocate(haloLayers(nHaloLayers))
        do iHalo = 1, nHaloLayers
          haloLayers(iHalo) = iHalo
        end do
      end if

#ifdef _MPI
      ! Allocate communication lists, and setup dead header nodes
      allocate(sendList)
      nullify(sendList % next)
      sendList % procID = -1
      sendList % nList = 0

      allocate(recvList)
      nullify(recvList % next)
      recvList % procID = -1
      recvList % nList = 0

      dminfo   = field % block % domain % dminfo

      ! Determine size of buffers for communication lists
      fieldCursor => field
      do while(associated(fieldCursor))

        ! Need to aggregate across halo layers
        do iHalo = 1, nHaloLayers
          
          ! Determine size from send lists
          exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
          do while(associated(exchListPtr))
            comm_list_found = .false.

            commListPtr => sendList
            do while(associated(commListPtr))
              if(commListPtr % procID == exchListPtr % endPointId) then
                comm_list_found = .true.
                commListPtr % nList = commListPtr % nList + exchListPtr % nList * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2)
                exit
              end if

              commListPtr => commListPtr % next
            end do

            if(.not. comm_list_found) then
              commListPtr => sendList
              commListPtr2 => commListPtr % next
              do while(associated(commListPtr2))
                commListPtr => commListPtr % next
                commListPtr2 => commListPtr % next
              end do

              allocate(commListPtr % next)
              commListPtr => commListPtr % next
              nullify(commListPtr % next)
              commListPtr % procID = exchListPtr % endPointID
              commListPtr % nList = exchListPtr % nList * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2)
            end if

            exchListPtr => exchListPtr % next
          end do

          ! Setup recv lists
          exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
          do while(associated(exchListPtr))
            comm_list_found = .false.

            commListPtr => recvList
            do while(associated(commListPtr))
              if(commListPtr % procID == exchListPtr % endPointId) then
                comm_list_found = .true.
                exit
              end if

              commListPtr => commListPtr % next
            end do

            if(.not. comm_list_found) then
              commListPtr => recvList
              commListPtr2 => commListPtr % next
              do while(associated(commListPtr2))
                commListPtr => commListPtr % next
                commListPtr2 => commListPtr % next
              end do

              allocate(commListPtr % next)
              commListPtr => commListPtr % next
              nullify(commListPtr % next)
              commListPtr % procID = exchListPtr % endPointID
            end if

            exchListPtr => exchListPtr % next
          end do
        end do

        fieldCursor => fieldCursor % next
      end do

      ! Remove the dead head pointer on send and recv list
      commListPtr => sendList
      sendList => sendList % next
      deallocate(commListPtr)

      commListPtr => recvList
      recvList => recvList % next
      deallocate(commListPtr)

      ! Determine size of recv lists
      commListPtr => recvList
      do while(associated(commListPtr))
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2))
              end if
              exchListPtr => exchListPtr % next
            end do
            
            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do
        commListPtr % nList = bufferOffset

        commListPtr => commListPtr % next
      end do

      ! Allocate space in recv lists, and initiate mpi_irecv calls
      commListPtr => recvList
      do while(associated(commListPtr))
        allocate(commListPtr % rbuffer(commListPtr % nList))
        nullify(commListPtr % ibuffer)
        call MPI_Irecv(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, commListPtr % procID, dminfo % comm, commListPtr % reqID, mpi_ierr)

        commListPtr => commListPtr % next
      end do

      ! Allocate space in send lists, copy data into buffer, and initiate mpi_isend calls
      commListPtr => sendList
      do while(associated(commListPtr))
        allocate(commListPtr % rbuffer(commListPtr % nList))
        nullify(commListPtr % ibuffer)
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                do  i = 1, exchListPtr % nList
                  do j = 1, fieldCursor % dimSizes(2)
                    do k = 1, fieldCursor % dimSizes(1)
                      commListPtr % rbuffer((exchListPtr % destList(i)-1) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                          + (j-1) * fieldCursor % dimSizes(1) + k  + bufferOffset) = fieldCursor % array(k, j, exchListPtr % srcList(i))
                      nAdded = nAdded + 1
                    end do
                  end do
                end do
              end if

              exchListPtr => exchListPtr % next
            end do

            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do

        call MPI_Isend(commListPtr % rbuffer, commListPtr % nList, MPI_REALKIND, commListPtr % procID, dminfo % my_proc_id, dminfo % comm, commListPtr % reqID, mpi_ierr)
        commListPtr => commListPtr % next
      end do
#endif

      ! Handle local copy. If MPI is off, then only local copies are performed.
      fieldCursor => field
      do while(associated(fieldCursor))
        do iHalo = 1, nHaloLayers
          exchListPtr => fieldCursor % copyList % halos(haloLayers(iHalo)) % exchList

          do while(associated(exchListPtr))
            fieldCursor2 => field
            do while(associated(fieldCursor2))
              if(exchListPtr % endPointID == fieldCursor2 % block % localBlockID) then
                do i = 1, exchListPtr % nList
                  fieldCursor2 % array(:, :, exchListPtr % destList(i)) = fieldCursor % array(:, :, exchListPtr % srcList(i))
                end do
              end if
              
              fieldCursor2 => fieldCursor2 % next
            end do

            exchListPtr => exchListPtr % next
          end do
        end do

        fieldCursor => fieldCursor % next
      end do

#ifdef _MPI

      ! Wait for mpi_irecv to finish, and unpack data from buffer
      commListPtr => recvList
      do while(associated(commListPtr))
        call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))
            exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if(exchListPtr % endPointID == commListPtr % procID) then
                do i = 1, exchListPtr % nList
                  do j = 1, fieldCursor % dimSizes(2)
                    do k = 1, fieldCursor % dimSizes(1)
                      fieldCursor % array(k, j, exchListPtr % destList(i)) = commListPtr % rbuffer((exchListPtr % srcList(i)-1)*fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2) &
                                                                           + (j-1)*fieldCursor % dimSizes(1) + k + bufferOffset)
                    end do
                  end do
                end do
                nAdded = max(nAdded, maxval(exchListPtr % srcList) * fieldCursor % dimSizes(1) * fieldCursor % dimSizes(2))
              end if
              exchListPtr => exchListPtr % next
            end do
            
            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do
        commListPtr => commListPtr % next
      end do

      ! wait for mpi_isend to finish.
      commListPtr => sendList
      do while(associated(commListPtr))
        call MPI_Wait(commListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
        commListPtr => commListPtr % next
      end do

     ! Destroy commLists.
     call mpas_dmpar_destroy_communication_list(sendList)
     call mpas_dmpar_destroy_communication_list(recvList)
#endif

     deallocate(haloLayers)

   end subroutine mpas_dmpar_exch_halo_field3d_real!}}}

   subroutine mpas_dmpar_init_mulithalo_exchange_list(exchList, nHalos)!{{{
     type (mpas_multihalo_exchange_list), pointer :: exchList
     integer, intent(in) :: nHalos

     integer :: i

     allocate(exchList)
     allocate(exchList % halos(nHalos))
     do i = 1, nHalos
       nullify(exchList % halos(i) % exchList)
     end do
   end subroutine mpas_dmpar_init_mulithalo_exchange_list!}}}

   subroutine mpas_dmpar_destroy_mulithalo_exchange_list(exchList)!{{{
     type (mpas_multihalo_exchange_list), pointer :: exchList

     integer :: nHalos
     integer :: i

     nHalos = size(exchList % halos)

     do i = 1, nHalos
       call mpas_dmpar_destroy_exchange_list(exchList % halos(i) % exchList)
     end do

     deallocate(exchList % halos)
     deallocate(exchList)
     nullify(exchList)
   end subroutine mpas_dmpar_destroy_mulithalo_exchange_list!}}}

   subroutine mpas_dmpar_destroy_communication_list(commList)!{{{
     type (mpas_communication_list), pointer :: commList
     type (mpas_communication_list), pointer :: commListPtr

     commListPtr => commList
     do while(associated(commListPtr))
       if(associated(commList % next)) then
         commList => commList % next
       else
         nullify(commList)
       end if

       if(associated(commListPtr % ibuffer)) then
         deallocate(commListPtr % ibuffer)
       end if

       if(associated(commListPtr % rbuffer)) then
         deallocate(commListPtr % rbuffer)
       end if

       deallocate(commListPtr)
       commListPtr => commList
     end do

   end subroutine mpas_dmpar_destroy_communication_list!}}}

   subroutine mpas_dmpar_destroy_exchange_list(exchList)!{{{
     type (mpas_exchange_list), pointer :: exchList
     type (mpas_exchange_list), pointer :: exchListPtr

     exchListPtr => exchList
     do while(associated(exchList))
       if(associated(exchList % next)) then
         exchList => exchList % next
       else
         nullify(exchList)
       end if

       if(associated(exchListPtr % srcList)) then
         deallocate(exchListPtr % srcList)
       end if

       if(associated(exchListPtr % destList)) then
         deallocate(exchListPtr % destList)
       end if

       deallocate(exchListPtr)
       exchListPtr => exchList
     end do

   end subroutine mpas_dmpar_destroy_exchange_list!}}}

   subroutine mpas_dmpar_copy_field1d_integer(field)!{{{
       type (field1dInteger), pointer :: field
       type (field1dInteger), pointer :: fieldCursor

       if(associated(field % next)) then
         fieldCursor => field % next
         do while(associated(fieldCursor))
           fieldCursor % array = field % array
           fieldCursor => fieldCursor % next
         end do
       end if
   end subroutine mpas_dmpar_copy_field1d_integer!}}}

   subroutine mpas_dmpar_copy_field2d_integer(field)!{{{
       type (field2dInteger), pointer :: field
       type (field2dInteger), pointer :: fieldCursor

       if(associated(field % next)) then
         fieldCursor => field % next
         do while(associated(fieldCursor))
           fieldCursor % array = field % array
           fieldCursor => fieldCursor % next
         end do
       end if
   end subroutine mpas_dmpar_copy_field2d_integer!}}}

   subroutine mpas_dmpar_copy_field3d_integer(field)!{{{
       type (field3dInteger), pointer :: field
       type (field3dInteger), pointer :: fieldCursor

       if(associated(field % next)) then
         fieldCursor => field % next
         do while(associated(fieldCursor))
           fieldCursor % array = field % array
           fieldCursor => fieldCursor % next
         end do
       end if
   end subroutine mpas_dmpar_copy_field3d_integer!}}}

   subroutine mpas_dmpar_copy_field1d_real(field)!{{{
       type (field1dReal), pointer :: field
       type (field1dReal), pointer :: fieldCursor


       if(associated(field % next)) then
         fieldCursor => field
         do while(associated(fieldCursor))
           fieldCursor % array(:) = field % array(:)
           fieldCursor => fieldCursor % next
         end do
       end if
   end subroutine mpas_dmpar_copy_field1d_real!}}}

   subroutine mpas_dmpar_copy_field2d_real(field)!{{{
       type (field2dReal), pointer :: field
       type (field2dReal), pointer :: fieldCursor

       if(associated(field % next)) then
         fieldCursor => field % next
         do while(associated(fieldCursor))
           fieldCursor % array = field % array
           fieldCursor => fieldCursor % next
         end do
       end if
   end subroutine mpas_dmpar_copy_field2d_real!}}}

   subroutine mpas_dmpar_copy_field3d_real(field)!{{{
       type (field3dReal), pointer :: field
       type (field3dReal), pointer :: fieldCursor

       if(associated(field % next)) then
         fieldCursor => field % next
         do while(associated(fieldCursor))
           fieldCursor % array = field % array
           fieldCursor => fieldCursor % next
         end do
       end if
   end subroutine mpas_dmpar_copy_field3d_real!}}}

end module mpas_dmpar
