module mpas_dmpar

   use mpas_kind_types
   use mpas_sort

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


   type dm_info
      integer :: nprocs, my_proc_id, comm, info
      logical :: using_external_comm
   end type dm_info


   type exchange_list
      integer :: procID
      integer :: nlist
      integer, dimension(:), pointer :: list
      type (exchange_list), pointer :: next
      real (kind=RKIND), dimension(:), pointer :: rbuffer
      integer, dimension(:), pointer           :: ibuffer
      integer :: reqID
   end type exchange_list


   interface mpas_dmpar_alltoall_field
      module procedure mpas_dmpar_alltoall_field1d_integer
      module procedure mpas_dmpar_alltoall_field2d_integer
      module procedure mpas_dmpar_alltoall_field1d_real
      module procedure mpas_dmpar_alltoall_field2d_real
      module procedure mpas_dmpar_alltoall_field3d_real
   end interface


   contains


   subroutine mpas_dmpar_init(dminfo, mpi_comm)

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

      call open_streams(dminfo % my_proc_id)

      dminfo % info = MPI_INFO_NULL
#else
      dminfo % comm = 0
      dminfo % my_proc_id = IO_NODE
      dminfo % nprocs = 1
      dminfo % using_external_comm = .false.
#endif

   end subroutine mpas_dmpar_init


   subroutine mpas_dmpar_finalize(dminfo)

      implicit none

      type (dm_info), intent(inout) :: dminfo

#ifdef _MPI
      integer :: mpi_ierr

      if (.not. dminfo % using_external_comm) then
         call MPI_Finalize(mpi_ierr)
      end if
#endif

   end subroutine mpas_dmpar_finalize


   subroutine mpas_dmpar_abort(dminfo)

      implicit none

      type (dm_info), intent(in) :: dminfo

#ifdef _MPI
      integer :: mpi_ierr, mpi_errcode

      call MPI_Abort(dminfo % comm, mpi_errcode, mpi_ierr)
#endif

      stop

   end subroutine mpas_dmpar_abort


   subroutine mpas_dmpar_global_abort(mesg)

      implicit none

      character (len=*), intent(in) :: mesg

#ifdef _MPI
      integer :: mpi_ierr, mpi_errcode

      write(0,*) trim(mesg)
      call MPI_Abort(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
#endif

      write(0,*) trim(mesg)
      stop

   end subroutine mpas_dmpar_global_abort


   subroutine mpas_dmpar_bcast_int(dminfo, i)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(inout) :: i

#ifdef _MPI
      integer :: mpi_ierr

      call MPI_Bcast(i, 1, MPI_INTEGERKIND, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_bcast_int


   subroutine mpas_dmpar_bcast_ints(dminfo, n, iarray)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: n
      integer, dimension(n), intent(inout) :: iarray

#ifdef _MPI
      integer :: mpi_ierr

      call MPI_Bcast(iarray, n, MPI_INTEGERKIND, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_bcast_ints


   subroutine mpas_dmpar_bcast_real(dminfo, r)

      implicit none

      type (dm_info), intent(in) :: dminfo
      real (kind=RKIND), intent(inout) :: r

#ifdef _MPI
      integer :: mpi_ierr

      call MPI_Bcast(r, 1, MPI_REALKIND, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_bcast_real


   subroutine mpas_dmpar_bcast_reals(dminfo, n, rarray)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: n
      real (kind=RKIND), dimension(n), intent(inout) :: rarray

#ifdef _MPI
      integer :: mpi_ierr

      call MPI_Bcast(rarray, n, MPI_REALKIND, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_bcast_reals


   subroutine mpas_dmpar_bcast_logical(dminfo, l)

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

   end subroutine mpas_dmpar_bcast_logical


   subroutine mpas_dmpar_bcast_char(dminfo, c)

      implicit none

      type (dm_info), intent(in) :: dminfo
      character (len=*), intent(inout) :: c

#ifdef _MPI
      integer :: mpi_ierr

      call MPI_Bcast(c, len(c), MPI_CHARACTER, IO_NODE, dminfo % comm, mpi_ierr)
#endif

   end subroutine mpas_dmpar_bcast_char


   subroutine mpas_dmpar_sum_int(dminfo, i, isum)

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

   end subroutine mpas_dmpar_sum_int


   subroutine mpas_dmpar_sum_real(dminfo, r, rsum)

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

   end subroutine mpas_dmpar_sum_real


   subroutine mpas_dmpar_min_int(dminfo, i, imin)

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

   end subroutine mpas_dmpar_min_int


   subroutine mpas_dmpar_min_real(dminfo, r, rmin)

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

   end subroutine mpas_dmpar_min_real


   subroutine mpas_dmpar_max_int(dminfo, i, imax)

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

   end subroutine mpas_dmpar_max_int


   subroutine mpas_dmpar_max_real(dminfo, r, rmax)

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

   end subroutine mpas_dmpar_max_real


   subroutine mpas_dmpar_sum_int_array(dminfo, nElements, inArray, outArray)

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

   end subroutine mpas_dmpar_sum_int_array


   subroutine mpas_dmpar_min_int_array(dminfo, nElements, inArray, outArray)
   
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

   end subroutine mpas_dmpar_min_int_array


   subroutine mpas_dmpar_max_int_array(dminfo, nElements, inArray, outArray)

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

   end subroutine mpas_dmpar_max_int_array


   subroutine mpas_dmpar_sum_real_array(dminfo, nElements, inArray, outArray)

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

   end subroutine mpas_dmpar_sum_real_array


   subroutine mpas_dmpar_min_real_array(dminfo, nElements, inArray, outArray)

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

   end subroutine mpas_dmpar_min_real_array


   subroutine mpas_dmpar_max_real_array(dminfo, nElements, inArray, outArray)

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

   end subroutine mpas_dmpar_max_real_array


   subroutine mpas_dmpar_scatter_ints(dminfo, nprocs, noutlist, displs, counts, inlist, outlist)

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

   end subroutine mpas_dmpar_scatter_ints


   subroutine mpas_dmpar_get_index_range(dminfo, &
                                    global_start, global_end, &
                                    local_start, local_end)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: global_start, global_end
      integer, intent(out) :: local_start, local_end

      local_start = nint(real(dminfo % my_proc_id) * real(global_end - global_start + 1) / real(dminfo % nprocs)) + 1
      local_end   = nint(real(dminfo % my_proc_id + 1) * real(global_end - global_start + 1) / real(dminfo % nprocs)) 

   end subroutine mpas_dmpar_get_index_range

  
   subroutine mpas_dmpar_compute_index_range(dminfo, &
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
      
   
   end subroutine mpas_dmpar_compute_index_range


   subroutine mpas_dmpar_get_owner_list(dminfo, &
                                   nOwnedList, nNeededList, &
                                   ownedList, neededList, &
                                   sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nOwnedList, nNeededList
      integer, dimension(nOwnedList), intent(in) :: ownedList
      integer, dimension(nNeededList), intent(in) :: neededList
      type (exchange_list), pointer :: sendList
      type (exchange_list), pointer :: recvList

      integer :: i, j, k, kk
      integer :: totalSize, nMesgRecv, nMesgSend, recvNeighbor, sendNeighbor, currentProc
      integer :: numToSend, numToRecv
      integer, dimension(nOwnedList) :: recipientList
      integer, dimension(2,nOwnedList) :: ownedListSorted
      integer, allocatable, dimension(:) :: ownerListIn, ownerListOut
      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: mpi_ierr, mpi_rreq, mpi_sreq

#ifdef _MPI
      allocate(sendList)
      allocate(recvList)
      nullify(sendList % next)
      nullify(recvList % next)
      sendListPtr => sendList
      recvListPtr => recvList

      do i=1,nOwnedList
         ownedListSorted(1,i) = ownedList(i)
         ownedListSorted(2,i) = i
      end do
      call quicksort(nOwnedList, ownedListSorted)

      call MPI_Allreduce(nNeededList, totalSize, 1, MPI_INTEGERKIND, MPI_MAX, dminfo % comm, mpi_ierr)

      allocate(ownerListIn(totalSize))
      allocate(ownerListOut(totalSize))

      nMesgRecv = nNeededList
      ownerListIn(1:nNeededList) = neededList(1:nNeededList)

      recvNeighbor = mod(dminfo % my_proc_id + dminfo % nprocs - 1, dminfo % nprocs)
      sendNeighbor = mod(dminfo % my_proc_id + 1, dminfo % nprocs)

      do i=1, dminfo % nprocs

         recipientList(:) = -1
         numToSend = 0

         currentProc = mod(dminfo % my_proc_id + dminfo % nprocs - i + 1, dminfo % nprocs)
         do j=1,nMesgRecv
            if (ownerListIn(j) > 0) then
               k = mpas_binary_search(ownedListSorted, 2, 1, nOwnedList, ownerListIn(j))
               if (k <= nOwnedList) then
                  ownerListOut(j) = -1 * dminfo % my_proc_id
                  numToSend = numToSend + 1
                  recipientList(ownedListSorted(2,k)) = numToSend
               else
                  ownerListOut(j) = ownerListIn(j)
               end if
            else
               ownerListOut(j) = ownerListIn(j)
            end if
         end do

         if (numToSend > 0) then
            allocate(sendListPtr % next)
            sendListPtr => sendListPtr % next
            sendListPtr % procID = currentProc
            sendListPtr % nlist = numToSend
            allocate(sendListPtr % list(numToSend))
            nullify(sendListPtr % next)
            kk = 1
            do j=1,nOwnedList
               if (recipientList(j) /= -1) then
                  sendListPtr % list(recipientList(j)) = j
                  kk = kk + 1
               end if
            end do
         end if

         nMesgSend = nMesgRecv
         call MPI_Irecv(nMesgRecv, 1, MPI_INTEGERKIND, recvNeighbor, i, dminfo % comm, mpi_rreq, mpi_ierr)
         call MPI_Isend(nMesgSend, 1, MPI_INTEGERKIND, sendNeighbor, i, dminfo % comm, mpi_sreq, mpi_ierr)
         call MPI_Wait(mpi_rreq, MPI_STATUS_IGNORE, mpi_ierr)
         call MPI_Wait(mpi_sreq, MPI_STATUS_IGNORE, mpi_ierr)
         call MPI_Irecv(ownerListIn, nMesgRecv, MPI_INTEGERKIND, recvNeighbor, i, dminfo % comm, mpi_rreq, mpi_ierr)
         call MPI_Isend(ownerListOut, nMesgSend, MPI_INTEGERKIND, sendNeighbor, i, dminfo % comm, mpi_sreq, mpi_ierr)
         call MPI_Wait(mpi_rreq, MPI_STATUS_IGNORE, mpi_ierr)
         call MPI_Wait(mpi_sreq, MPI_STATUS_IGNORE, mpi_ierr)
      end do

      do i=0, dminfo % nprocs - 1

         numToRecv = 0
         do j=1,nNeededList
            if (ownerListIn(j) == -i) numToRecv = numToRecv + 1
         end do
         if (numToRecv > 0) then
            allocate(recvListPtr % next)
            recvListPtr => recvListPtr % next
            recvListPtr % procID = i
            recvListPtr % nlist = numToRecv
            allocate(recvListPtr % list(numToRecv))
            nullify(recvListPtr % next)
            kk = 1
            do j=1,nNeededList
               if (ownerListIn(j) == -i) then
                  recvListPtr % list(kk) = j
                  kk = kk + 1
               end if
            end do
         end if

      end do

      deallocate(ownerListIn)
      deallocate(ownerListOut)

      sendListPtr => sendList
      sendList => sendList % next
      deallocate(sendListPtr)

      recvListPtr => recvList
      recvList => recvList % next
      deallocate(recvListPtr)

#else
      allocate(recvList)
      recvList % procID = dminfo % my_proc_id
      recvList % nlist = nNeededList
      allocate(recvList % list(nNeededList))
      nullify(recvList % next)
      do j=1,nNeededList
         recvList % list(j) = j
      end do

      allocate(sendList)
      sendList % procID = dminfo % my_proc_id
      sendList % nlist = nOwnedList
      allocate(sendList % list(nOwnedList))
      nullify(sendList % next)
      do j=1,nOwnedList
         sendList % list(j) = j
      end do
#endif

   end subroutine mpas_dmpar_get_owner_list


   subroutine mpas_dmpar_alltoall_field1d_integer(dminfo, arrayIn, arrayOut, nOwnedList, nNeededList, sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, dimension(*), intent(in) :: arrayIn
      integer, dimension(*), intent(inout) :: arrayOut
      integer, intent(in) :: nOwnedList, nNeededList
      type (exchange_list), pointer :: sendList, recvList

      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
      integer :: mpi_ierr
      integer :: i

#ifdef _MPI

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID == dminfo % my_proc_id) exit
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID == dminfo % my_proc_id) exit
         recvListPtr => recvListPtr % next
      end do

      if (associated(recvListPtr) .and. associated(sendListPtr)) then
         do i=1,recvListPtr % nlist
            arrayOut(recvListPtr % list(i)) = arrayIn(sendListPtr % list(i))
         end do
      end if

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            allocate(recvListPtr % ibuffer(recvListPtr % nlist))
            call MPI_Irecv(recvListPtr % ibuffer, recvListPtr % nlist, MPI_INTEGERKIND, &
                           recvListPtr % procID, recvListPtr % procID, dminfo % comm, recvListPtr % reqID, mpi_ierr)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            allocate(sendListPtr % ibuffer(sendListPtr % nlist))
            call mpas_pack_send_buf1d_integer(nOwnedList, arrayIn, sendListPtr, 1, sendListPtr % nlist, &
                                   sendListPtr % ibuffer, nPacked, lastPackedIdx)
            call MPI_Isend(sendListPtr % ibuffer, sendListPtr % nlist, MPI_INTEGERKIND, &
                           sendListPtr % procID, dminfo % my_proc_id, dminfo % comm, sendListPtr % reqID, mpi_ierr)
         end if
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(recvListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            call mpas_unpack_recv_buf1d_integer(nNeededList, arrayOut, recvListPtr, 1, recvListPtr % nlist, &
                                     recvListPtr % ibuffer, nUnpacked, lastUnpackedIdx)
            deallocate(recvListPtr % ibuffer)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(sendListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            deallocate(sendListPtr % ibuffer)
         end if
         sendListPtr => sendListPtr % next
      end do

#else
      if (nOwnedList /= nNeededList) then
         write(0,*) 'Error in dmpar_alltoall_field: For non-dmpar, ', &
           'arrayIn and arrayOut dims must match.'
         call mpas_dmpar_abort(dminfo)
      else
         arrayOut(1:nNeededList) = arrayIn(1:nOwnedList)
      end if
#endif

   end subroutine mpas_dmpar_alltoall_field1d_integer


   subroutine mpas_dmpar_alltoall_field2d_integer(dminfo, arrayIn, arrayOut, dim1, nOwnedList, nNeededList, sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: dim1, nOwnedList, nNeededList
      integer, dimension(dim1,*), intent(in) :: arrayIn
      integer, dimension(dim1,*), intent(inout) :: arrayOut
      type (exchange_list), pointer :: sendList, recvList

      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
      integer :: mpi_ierr
      integer :: i, d2

#ifdef _MPI

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID == dminfo % my_proc_id) exit
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID == dminfo % my_proc_id) exit
         recvListPtr => recvListPtr % next
      end do

      if (associated(recvListPtr) .and. associated(sendListPtr)) then
         do i=1,recvListPtr % nlist
            arrayOut(:,recvListPtr % list(i)) = arrayIn(:,sendListPtr % list(i))
         end do
      end if

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            d2 = dim1 * recvListPtr % nlist
            allocate(recvListPtr % ibuffer(d2))
            call MPI_Irecv(recvListPtr % ibuffer, d2, MPI_INTEGERKIND, &
                           recvListPtr % procID, recvListPtr % procID, dminfo % comm, recvListPtr % reqID, mpi_ierr)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            d2 = dim1 * sendListPtr % nlist
            allocate(sendListPtr % ibuffer(d2))
            call mpas_pack_send_buf2d_integer(1, dim1, nOwnedList, arrayIn, sendListPtr, 1, d2, &
                                   sendListPtr % ibuffer, nPacked, lastPackedIdx)
            call MPI_Isend(sendListPtr % ibuffer, d2, MPI_INTEGERKIND, &
                           sendListPtr % procID, dminfo % my_proc_id, dminfo % comm, sendListPtr % reqID, mpi_ierr)
         end if
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(recvListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            d2 = dim1 * recvListPtr % nlist
            call mpas_unpack_recv_buf2d_integer(1, dim1, nNeededList, arrayOut, recvListPtr, 1, d2, &
                                     recvListPtr % ibuffer, nUnpacked, lastUnpackedIdx)
            deallocate(recvListPtr % ibuffer)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(sendListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            deallocate(sendListPtr % ibuffer)
         end if
         sendListPtr => sendListPtr % next
      end do

#else
      if (nOwnedList /= nNeededList) then
         write(0,*) 'Error in dmpar_alltoall_field: For non-dmpar, ', &
           'arrayIn and arrayOut dims must match.'
         call mpas_dmpar_abort(dminfo)
      else
         arrayOut(:,1:nNeededList) = arrayIn(:,1:nOwnedList)
      end if
#endif

   end subroutine mpas_dmpar_alltoall_field2d_integer


   subroutine mpas_dmpar_alltoall_field1d_real(dminfo, arrayIn, arrayOut, nOwnedList, nNeededList, sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      real (kind=RKIND), dimension(*), intent(in) :: arrayIn
      real (kind=RKIND), dimension(*), intent(inout) :: arrayOut
      integer, intent(in) :: nOwnedList, nNeededList
      type (exchange_list), pointer :: sendList, recvList

      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
      integer :: mpi_ierr
      integer :: i

#ifdef _MPI

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID == dminfo % my_proc_id) exit
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID == dminfo % my_proc_id) exit
         recvListPtr => recvListPtr % next
      end do

      if (associated(recvListPtr) .and. associated(sendListPtr)) then
         do i=1,recvListPtr % nlist
            arrayOut(recvListPtr % list(i)) = arrayIn(sendListPtr % list(i))
         end do
      end if

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            allocate(recvListPtr % rbuffer(recvListPtr % nlist))
            call MPI_Irecv(recvListPtr % rbuffer, recvListPtr % nlist, MPI_REALKIND, &
                           recvListPtr % procID, recvListPtr % procID, dminfo % comm, recvListPtr % reqID, mpi_ierr)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            allocate(sendListPtr % rbuffer(sendListPtr % nlist))
            call mpas_pack_send_buf1d_real(nOwnedList, arrayIn, sendListPtr, 1, sendListPtr % nlist, &
                                   sendListPtr % rbuffer, nPacked, lastPackedIdx)
            call MPI_Isend(sendListPtr % rbuffer, sendListPtr % nlist, MPI_REALKIND, &
                           sendListPtr % procID, dminfo % my_proc_id, dminfo % comm, sendListPtr % reqID, mpi_ierr)
         end if
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(recvListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            call mpas_unpack_recv_buf1d_real(nNeededList, arrayOut, recvListPtr, 1, recvListPtr % nlist, &
                                     recvListPtr % rbuffer, nUnpacked, lastUnpackedIdx)
            deallocate(recvListPtr % rbuffer)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(sendListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            deallocate(sendListPtr % rbuffer)
         end if
         sendListPtr => sendListPtr % next
      end do

#else
      if (nOwnedList /= nNeededList) then
         write(0,*) 'Error in dmpar_alltoall_field: For non-dmpar, ', &
           'arrayIn and arrayOut dims must match.'
         call mpas_dmpar_abort(dminfo)
      else
         arrayOut(1:nNeededList) = arrayIn(1:nOwnedList)
      end if
#endif

   end subroutine mpas_dmpar_alltoall_field1d_real


   subroutine mpas_dmpar_alltoall_field2d_real(dminfo, arrayIn, arrayOut, dim1, nOwnedList, nNeededList, sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: dim1, nOwnedList, nNeededList
      real (kind=RKIND), dimension(dim1,*), intent(in) :: arrayIn
      real (kind=RKIND), dimension(dim1,*), intent(inout) :: arrayOut
      type (exchange_list), pointer :: sendList, recvList

      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
      integer :: mpi_ierr
      integer :: i, d2

#ifdef _MPI

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID == dminfo % my_proc_id) exit
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID == dminfo % my_proc_id) exit
         recvListPtr => recvListPtr % next
      end do

      if (associated(recvListPtr) .and. associated(sendListPtr)) then
         do i=1,recvListPtr % nlist
            arrayOut(:,recvListPtr % list(i)) = arrayIn(:,sendListPtr % list(i))
         end do
      end if

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            d2 = dim1 * recvListPtr % nlist
            allocate(recvListPtr % rbuffer(d2))
            call MPI_Irecv(recvListPtr % rbuffer, d2, MPI_REALKIND, &
                           recvListPtr % procID, recvListPtr % procID, dminfo % comm, recvListPtr % reqID, mpi_ierr)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            d2 = dim1 * sendListPtr % nlist
            allocate(sendListPtr % rbuffer(d2))
            call mpas_pack_send_buf2d_real(1, dim1, nOwnedList, arrayIn, sendListPtr, 1, d2, &
                                   sendListPtr % rbuffer, nPacked, lastPackedIdx)
            call MPI_Isend(sendListPtr % rbuffer, d2, MPI_REALKIND, &
                           sendListPtr % procID, dminfo % my_proc_id, dminfo % comm, sendListPtr % reqID, mpi_ierr)
         end if
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(recvListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            d2 = dim1 * recvListPtr % nlist
            call mpas_unpack_recv_buf2d_real(1, dim1, nNeededList, arrayOut, recvListPtr, 1, d2, &
                                     recvListPtr % rbuffer, nUnpacked, lastUnpackedIdx)
            deallocate(recvListPtr % rbuffer)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(sendListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            deallocate(sendListPtr % rbuffer)
         end if
         sendListPtr => sendListPtr % next
      end do

#else
      if (nOwnedList /= nNeededList) then
         write(0,*) 'Error in dmpar_alltoall_field: For non-dmpar, ', &
           'arrayIn and arrayOut dims must match.'
         call mpas_dmpar_abort(dminfo)
      else
         arrayOut(:,1:nNeededList) = arrayIn(:,1:nOwnedList)
      end if
#endif

   end subroutine mpas_dmpar_alltoall_field2d_real

  
   subroutine mpas_dmpar_alltoall_field3d_real(dminfo, arrayIn, arrayOut, dim1, dim2, nOwnedList, nNeededList, sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: dim1, dim2, nOwnedList, nNeededList
      real (kind=RKIND), dimension(dim1,dim2,*), intent(in) :: arrayIn
      real (kind=RKIND), dimension(dim1,dim2,*), intent(inout) :: arrayOut
      type (exchange_list), pointer :: sendList, recvList

      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
      integer :: mpi_ierr
      integer :: i, d3

#ifdef _MPI

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID == dminfo % my_proc_id) exit
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID == dminfo % my_proc_id) exit
         recvListPtr => recvListPtr % next
      end do

      if (associated(recvListPtr) .and. associated(sendListPtr)) then
         do i=1,recvListPtr % nlist
            arrayOut(:,:,recvListPtr % list(i)) = arrayIn(:,:,sendListPtr % list(i))
         end do
      end if

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            d3 = dim1 * dim2 * recvListPtr % nlist
            allocate(recvListPtr % rbuffer(d3))
            call MPI_Irecv(recvListPtr % rbuffer, d3, MPI_REALKIND, &
                           recvListPtr % procID, recvListPtr % procID, dminfo % comm, recvListPtr % reqID, mpi_ierr)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            d3 = dim1 * dim2 * sendListPtr % nlist
            allocate(sendListPtr % rbuffer(d3))
            call mpas_pack_send_buf3d_real(1, dim1, 1, dim2, nOwnedList, arrayIn, sendListPtr, 1, d3, &
                                   sendListPtr % rbuffer, nPacked, lastPackedIdx)
            call MPI_Isend(sendListPtr % rbuffer, d3, MPI_REALKIND, &
                           sendListPtr % procID, dminfo % my_proc_id, dminfo % comm, sendListPtr % reqID, mpi_ierr)
         end if
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(recvListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            d3 = dim1 * dim2 * recvListPtr % nlist
            call mpas_unpack_recv_buf3d_real(1, dim1, 1, dim2, nNeededList, arrayOut, recvListPtr, 1, d3, &
                                     recvListPtr % rbuffer, nUnpacked, lastUnpackedIdx)
            deallocate(recvListPtr % rbuffer)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(sendListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            deallocate(sendListPtr % rbuffer)
         end if
         sendListPtr => sendListPtr % next
      end do

#else
      if (nOwnedList /= nNeededList) then
         write(0,*) 'Error in dmpar_alltoall_field: For non-dmpar, ', &
           'arrayIn and arrayOut dims must match.'
         call mpas_dmpar_abort(dminfo)
      else
         arrayOut(:,:,1:nNeededList) = arrayIn(:,:,1:nOwnedList)
      end if
#endif

   end subroutine mpas_dmpar_alltoall_field3d_real

  
   subroutine mpas_pack_send_buf1d_integer(nField, field, sendList, startPackIdx, nBuffer, buffer, nPacked, lastPackedIdx)

      implicit none

      integer, intent(in) :: nField, nBuffer, startPackIdx
      integer, dimension(*), intent(in) :: field
      type (exchange_list), intent(in) :: sendList
      integer, dimension(nBuffer), intent(out) :: buffer
      integer, intent(inout) :: nPacked, lastPackedIdx

      integer :: i

      nPacked = 0
      do i=startPackIdx, sendList % nlist
         nPacked = nPacked + 1
         if (nPacked > nBuffer) then
            nPacked = nPacked - 1
            lastPackedIdx = i - 1
            return
         end if
         buffer(nPacked) = field(sendList % list(i))
      end do
      lastPackedIdx = sendList % nlist

   end subroutine mpas_pack_send_buf1d_integer


   subroutine mpas_pack_send_buf2d_integer(ds, de, nField, field, sendList, startPackIdx, nBuffer, buffer, nPacked, lastPackedIdx)

      implicit none

      integer, intent(in) :: ds, de, nField, nBuffer, startPackIdx
      integer, dimension(ds:de,*), intent(in) :: field
      type (exchange_list), intent(in) :: sendList
      integer, dimension(nBuffer), intent(out) :: buffer
      integer, intent(inout) :: nPacked, lastPackedIdx

      integer :: i, n

      n = de-ds+1

      if (n > nBuffer) then
         write(0,*) 'packSendBuf2dInteger: Not enough space in buffer', &
          ' to fit a single slice.'
         return
      end if

      nPacked = 0
      do i=startPackIdx, sendList % nlist
         nPacked = nPacked + n
         if (nPacked > nBuffer) then
            nPacked = nPacked - n
            lastPackedIdx = i - 1
            return
         end if
         buffer(nPacked-n+1:nPacked) = field(ds:de,sendList % list(i))
      end do
      lastPackedIdx = sendList % nlist

   end subroutine mpas_pack_send_buf2d_integer


   subroutine mpas_pack_send_buf3d_integer(d1s, d1e, d2s, d2e, nField, field, sendList, startPackIdx, nBuffer, buffer, nPacked, lastPackedIdx)

      implicit none

      integer, intent(in) :: d1s, d1e, d2s, d2e, nField, nBuffer, startPackIdx
      integer, dimension(d1s:d1e,d2s:d2e,*), intent(in) :: field
      type (exchange_list), intent(in) :: sendList
      integer, dimension(nBuffer), intent(out) :: buffer
      integer, intent(inout) :: nPacked, lastPackedIdx

      integer :: i, j, k, n

      n = (d1e-d1s+1) * (d2e-d2s+1)

      if (n > nBuffer) then
         write(0,*) 'packSendBuf3dInteger: Not enough space in buffer', &
          ' to fit a single slice.'
         return
      end if

      nPacked = 0
      do i=startPackIdx, sendList % nlist
         nPacked = nPacked + n
         if (nPacked > nBuffer) then
            nPacked = nPacked - n
            lastPackedIdx = i - 1
            return
         end if
         k = nPacked-n+1
         do j=d2s,d2e
            buffer(k:k+d1e-d1s) = field(d1s:d1e,j,sendList % list(i))
            k = k + d1e-d1s+1
         end do
      end do
      lastPackedIdx = sendList % nlist

   end subroutine mpas_pack_send_buf3d_integer


   subroutine mpas_pack_send_buf1d_real(nField, field, sendList, startPackIdx, nBuffer, buffer, nPacked, lastPackedIdx)

      implicit none

      integer, intent(in) :: nField, nBuffer, startPackIdx
      real (kind=RKIND), dimension(*), intent(in) :: field
      type (exchange_list), intent(in) :: sendList
      real (kind=RKIND), dimension(nBuffer), intent(out) :: buffer
      integer, intent(inout) :: nPacked, lastPackedIdx

      integer :: i

      nPacked = 0
      do i=startPackIdx, sendList % nlist
         nPacked = nPacked + 1
         if (nPacked > nBuffer) then
            nPacked = nPacked - 1
            lastPackedIdx = i - 1
            return
         end if
         buffer(nPacked) = field(sendList % list(i))
      end do
      lastPackedIdx = sendList % nlist

   end subroutine mpas_pack_send_buf1d_real


   subroutine mpas_pack_send_buf2d_real(ds, de, nField, field, sendList, startPackIdx, nBuffer, buffer, nPacked, lastPackedIdx)

      implicit none

      integer, intent(in) :: ds, de, nField, nBuffer, startPackIdx
      real (kind=RKIND), dimension(ds:de,*), intent(in) :: field
      type (exchange_list), intent(in) :: sendList
      real (kind=RKIND), dimension(nBuffer), intent(out) :: buffer
      integer, intent(inout) :: nPacked, lastPackedIdx

      integer :: i, n

      n = de-ds+1

      if (n > nBuffer) then
         write(0,*) 'packSendBuf2dReal: Not enough space in buffer', &
          ' to fit a single slice.'
         return
      end if

      nPacked = 0
      do i=startPackIdx, sendList % nlist
         nPacked = nPacked + n
         if (nPacked > nBuffer) then
            nPacked = nPacked - n
            lastPackedIdx = i - 1
            return
         end if
         buffer(nPacked-n+1:nPacked) = field(ds:de,sendList % list(i))
      end do
      lastPackedIdx = sendList % nlist

   end subroutine mpas_pack_send_buf2d_real


   subroutine mpas_pack_send_buf3d_real(d1s, d1e, d2s, d2e, nField, field, sendList, startPackIdx, nBuffer, buffer, nPacked, lastPackedIdx)

      implicit none

      integer, intent(in) :: d1s, d1e, d2s, d2e, nField, nBuffer, startPackIdx
      real (kind=RKIND), dimension(d1s:d1e,d2s:d2e,*), intent(in) :: field
      type (exchange_list), intent(in) :: sendList
      real (kind=RKIND), dimension(nBuffer), intent(out) :: buffer
      integer, intent(inout) :: nPacked, lastPackedIdx

      integer :: i, j, k, n

      n = (d1e-d1s+1) * (d2e-d2s+1)

      if (n > nBuffer) then
         write(0,*) 'packSendBuf3dReal: Not enough space in buffer', &
          ' to fit a single slice.'
         return
      end if

      nPacked = 0
      do i=startPackIdx, sendList % nlist
         nPacked = nPacked + n
         if (nPacked > nBuffer) then
            nPacked = nPacked - n
            lastPackedIdx = i - 1
            return
         end if
         k = nPacked-n+1
         do j=d2s,d2e
            buffer(k:k+d1e-d1s) = field(d1s:d1e,j,sendList % list(i))
            k = k + d1e-d1s+1
         end do
      end do
      lastPackedIdx = sendList % nlist

   end subroutine mpas_pack_send_buf3d_real


   subroutine mpas_unpack_recv_buf1d_integer(nField, field, recvList, startUnpackIdx, nBuffer, buffer, nUnpacked, lastUnpackedIdx)

      implicit none

      integer, intent(in) :: nField, nBuffer, startUnpackIdx
      integer, dimension(*), intent(inout) :: field
      type (exchange_list), intent(in) :: recvList
      integer, dimension(nBuffer), intent(in) :: buffer
      integer, intent(inout) :: nUnpacked, lastUnpackedIdx

      integer :: i

      nUnpacked = 0
      do i=startUnpackIdx, recvList % nlist
         nUnpacked = nUnpacked + 1
         if (nUnpacked > nBuffer) then
            nUnpacked = nUnpacked - 1
            lastUnpackedIdx = i - 1
            return
         end if
         field(recvList % list(i)) = buffer(nUnpacked)
      end do
      lastUnpackedIdx = recvList % nlist

   end subroutine mpas_unpack_recv_buf1d_integer


   subroutine mpas_unpack_recv_buf2d_integer(ds, de, nField, field, recvList, startUnpackIdx, nBuffer, buffer, nUnpacked, lastUnpackedIdx)

      implicit none

      integer, intent(in) :: ds, de, nField, nBuffer, startUnpackIdx
      integer, dimension(ds:de,*), intent(inout) :: field
      type (exchange_list), intent(in) :: recvList
      integer, dimension(nBuffer), intent(in) :: buffer
      integer, intent(inout) :: nUnpacked, lastUnpackedIdx

      integer :: i, n

      n = de-ds+1

      nUnpacked = 0
      do i=startUnpackIdx, recvList % nlist
         nUnpacked = nUnpacked + n
         if (nUnpacked > nBuffer) then
            nUnpacked = nUnpacked - n
            lastUnpackedIdx = i - 1
            return
         end if
         field(ds:de,recvList % list(i)) = buffer(nUnpacked-n+1:nUnpacked)
      end do
      lastUnpackedIdx = recvList % nlist

   end subroutine mpas_unpack_recv_buf2d_integer


   subroutine mpas_unpack_recv_buf3d_integer(d1s, d1e, d2s, d2e, nField, field, recvList, startUnpackIdx, nBuffer, buffer, &
                                  nUnpacked, lastUnpackedIdx)

      implicit none

      integer, intent(in) :: d1s, d1e, d2s, d2e, nField, nBuffer, startUnpackIdx
      integer, dimension(d1s:d1e,d2s:d2e,*), intent(inout) :: field
      type (exchange_list), intent(in) :: recvList
      integer, dimension(nBuffer), intent(in) :: buffer
      integer, intent(inout) :: nUnpacked, lastUnpackedIdx

      integer :: i, j, k, n

      n = (d1e-d1s+1) * (d2e-d2s+1)

      nUnpacked = 0
      do i=startUnpackIdx, recvList % nlist
         nUnpacked = nUnpacked + n
         if (nUnpacked > nBuffer) then
            nUnpacked = nUnpacked - n
            lastUnpackedIdx = i - 1
            return
         end if
         k = nUnpacked-n+1
         do j=d2s,d2e
            field(d1s:d1e,j,recvList % list(i)) = buffer(k:k+d1e-d1s)
            k = k + d1e-d1s+1
         end do
      end do
      lastUnpackedIdx = recvList % nlist

   end subroutine mpas_unpack_recv_buf3d_integer


   subroutine mpas_dmpar_exch_halo_field1d_integer(dminfo, array, dim1, sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: dim1
      integer, dimension(*), intent(inout) :: array
      type (exchange_list), pointer :: sendList, recvList

      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
      integer :: mpi_ierr

#ifdef _MPI

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            allocate(recvListPtr % ibuffer(recvListPtr % nlist))
            call MPI_Irecv(recvListPtr % ibuffer, recvListPtr % nlist, MPI_INTEGERKIND, &
                           recvListPtr % procID, recvListPtr % procID, dminfo % comm, recvListPtr % reqID, mpi_ierr)
         end if
         recvListPtr => recvListPtr % next
      end do
      
      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            allocate(sendListPtr % ibuffer(sendListPtr % nlist))
            call mpas_pack_send_buf1d_integer(dim1, array, sendListPtr, 1, sendListPtr % nlist, sendListPtr % ibuffer, nPacked, lastPackedIdx)
            call MPI_Isend(sendListPtr % ibuffer, sendListPtr % nlist, MPI_INTEGERKIND, &
                           sendListPtr % procID, dminfo % my_proc_id, dminfo % comm, sendListPtr % reqID, mpi_ierr)
         end if
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(recvListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            call mpas_unpack_recv_buf1d_integer(dim1, array, recvListPtr, 1, recvListPtr % nlist, recvListPtr % ibuffer, nUnpacked, lastUnpackedIdx)
            deallocate(recvListPtr % ibuffer)
         end if
         recvListPtr => recvListPtr % next
      end do
      
      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(sendListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            deallocate(sendListPtr % ibuffer)
         end if
         sendListPtr => sendListPtr % next
      end do

#endif

   end subroutine mpas_dmpar_exch_halo_field1d_integer


   subroutine mpas_dmpar_exch_halo_field2d_integer(dminfo, array, dim1, dim2, sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: dim1, dim2
      integer, dimension(dim1,*), intent(inout) :: array
      type (exchange_list), pointer :: sendList, recvList

      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
      integer :: mpi_ierr
      integer :: d2

#ifdef _MPI

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            d2 = dim1 * recvListPtr % nlist
            allocate(recvListPtr % ibuffer(d2))
            call MPI_Irecv(recvListPtr % ibuffer, d2, MPI_INTEGERKIND, &
                           recvListPtr % procID, recvListPtr % procID, dminfo % comm, recvListPtr % reqID, mpi_ierr)
         end if
         recvListPtr => recvListPtr % next
      end do
      
      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            d2 = dim1 * sendListPtr % nlist
            allocate(sendListPtr % ibuffer(d2))
            call mpas_pack_send_buf2d_integer(1, dim1, dim2, array, sendListPtr, 1, d2, sendListPtr % ibuffer, nPacked, lastPackedIdx)
            call MPI_Isend(sendListPtr % ibuffer, d2, MPI_INTEGERKIND, &
                           sendListPtr % procID, dminfo % my_proc_id, dminfo % comm, sendListPtr % reqID, mpi_ierr)
         end if
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(recvListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            d2 = dim1 * recvListPtr % nlist
            call mpas_unpack_recv_buf2d_integer(1, dim1, dim2, array, recvListPtr, 1, d2, recvListPtr % ibuffer, nUnpacked, lastUnpackedIdx)
            deallocate(recvListPtr % ibuffer)
         end if
         recvListPtr => recvListPtr % next
      end do
      
      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(sendListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            deallocate(sendListPtr % ibuffer)
         end if
         sendListPtr => sendListPtr % next
      end do

#endif

   end subroutine mpas_dmpar_exch_halo_field2d_integer


   subroutine mpas_dmpar_exch_halo_field3d_integer(dminfo, array, dim1, dim2, dim3, sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: dim1, dim2, dim3
      integer, dimension(dim1,dim2,*), intent(inout) :: array
      type (exchange_list), pointer :: sendList, recvList

      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
      integer :: mpi_ierr
      integer :: d3

#ifdef _MPI

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            d3 = dim1 * dim2 * recvListPtr % nlist
            allocate(recvListPtr % ibuffer(d3))
            call MPI_Irecv(recvListPtr % ibuffer, d3, MPI_INTEGERKIND, &
                           recvListPtr % procID, recvListPtr % procID, dminfo % comm, recvListPtr % reqID, mpi_ierr)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            d3 = dim1 * dim2 * sendListPtr % nlist
            allocate(sendListPtr % ibuffer(d3))
            call mpas_pack_send_buf3d_integer(1, dim1, 1, dim2, dim3, array, sendListPtr, 1, d3, &
                                   sendListPtr % ibuffer, nPacked, lastPackedIdx)
            call MPI_Isend(sendListPtr % ibuffer, d3, MPI_INTEGERKIND, &
                           sendListPtr % procID, dminfo % my_proc_id, dminfo % comm, sendListPtr % reqID, mpi_ierr)
         end if
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(recvListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            d3 = dim1 * dim2 * recvListPtr % nlist
            call mpas_unpack_recv_buf3d_integer(1, dim1, 1, dim2, dim3, array, recvListPtr, 1, d3, &
                                     recvListPtr % ibuffer, nUnpacked, lastUnpackedIdx)
            deallocate(recvListPtr % ibuffer)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(sendListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            deallocate(sendListPtr % ibuffer)
         end if
         sendListPtr => sendListPtr % next
      end do

#endif

   end subroutine mpas_dmpar_exch_halo_field3d_integer

  
   subroutine mpas_unpack_recv_buf1d_real(nField, field, recvList, startUnpackIdx, nBuffer, buffer, nUnpacked, lastUnpackedIdx)

      implicit none

      integer, intent(in) :: nField, nBuffer, startUnpackIdx
      real (kind=RKIND), dimension(*), intent(inout) :: field
      type (exchange_list), intent(in) :: recvList
      real (kind=RKIND), dimension(nBuffer), intent(in) :: buffer
      integer, intent(inout) :: nUnpacked, lastUnpackedIdx

      integer :: i

      nUnpacked = 0
      do i=startUnpackIdx, recvList % nlist
         nUnpacked = nUnpacked + 1
         if (nUnpacked > nBuffer) then
            nUnpacked = nUnpacked - 1
            lastUnpackedIdx = i - 1
            return
         end if
         field(recvList % list(i)) = buffer(nUnpacked)
      end do
      lastUnpackedIdx = recvList % nlist

   end subroutine mpas_unpack_recv_buf1d_real


   subroutine mpas_unpack_recv_buf2d_real(ds, de, nField, field, recvList, startUnpackIdx, nBuffer, buffer, nUnpacked, lastUnpackedIdx)

      implicit none

      integer, intent(in) :: ds, de, nField, nBuffer, startUnpackIdx
      real (kind=RKIND), dimension(ds:de,*), intent(inout) :: field
      type (exchange_list), intent(in) :: recvList
      real (kind=RKIND), dimension(nBuffer), intent(in) :: buffer
      integer, intent(inout) :: nUnpacked, lastUnpackedIdx

      integer :: i, n

      n = de-ds+1

      nUnpacked = 0
      do i=startUnpackIdx, recvList % nlist
         nUnpacked = nUnpacked + n
         if (nUnpacked > nBuffer) then
            nUnpacked = nUnpacked - n
            lastUnpackedIdx = i - 1
            return
         end if
         field(ds:de,recvList % list(i)) = buffer(nUnpacked-n+1:nUnpacked)
      end do
      lastUnpackedIdx = recvList % nlist

   end subroutine mpas_unpack_recv_buf2d_real


   subroutine mpas_unpack_recv_buf3d_real(d1s, d1e, d2s, d2e, nField, field, recvList, startUnpackIdx, nBuffer, buffer, &
                                  nUnpacked, lastUnpackedIdx)

      implicit none

      integer, intent(in) :: d1s, d1e, d2s, d2e, nField, nBuffer, startUnpackIdx
      real (kind=RKIND), dimension(d1s:d1e,d2s:d2e,*), intent(inout) :: field
      type (exchange_list), intent(in) :: recvList
      real (kind=RKIND), dimension(nBuffer), intent(in) :: buffer
      integer, intent(inout) :: nUnpacked, lastUnpackedIdx

      integer :: i, j, k, n

      n = (d1e-d1s+1) * (d2e-d2s+1)

      nUnpacked = 0
      do i=startUnpackIdx, recvList % nlist
         nUnpacked = nUnpacked + n
         if (nUnpacked > nBuffer) then
            nUnpacked = nUnpacked - n
            lastUnpackedIdx = i - 1
            return
         end if
         k = nUnpacked-n+1
         do j=d2s,d2e
            field(d1s:d1e,j,recvList % list(i)) = buffer(k:k+d1e-d1s)
            k = k + d1e-d1s+1
         end do
      end do
      lastUnpackedIdx = recvList % nlist

   end subroutine mpas_unpack_recv_buf3d_real


   subroutine mpas_dmpar_exch_halo_field1d_real(dminfo, array, dim1, sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: dim1
      real (kind=RKIND), dimension(*), intent(inout) :: array
      type (exchange_list), pointer :: sendList, recvList

      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
      integer :: mpi_ierr

#ifdef _MPI

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            allocate(recvListPtr % rbuffer(recvListPtr % nlist))
            call MPI_Irecv(recvListPtr % rbuffer, recvListPtr % nlist, MPI_REALKIND, &
                           recvListPtr % procID, recvListPtr % procID, dminfo % comm, recvListPtr % reqID, mpi_ierr)
         end if
         recvListPtr => recvListPtr % next
      end do
      
      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            allocate(sendListPtr % rbuffer(sendListPtr % nlist))
            call mpas_pack_send_buf1d_real(dim1, array, sendListPtr, 1, sendListPtr % nlist, sendListPtr % rbuffer, nPacked, lastPackedIdx)
            call MPI_Isend(sendListPtr % rbuffer, sendListPtr % nlist, MPI_REALKIND, &
                           sendListPtr % procID, dminfo % my_proc_id, dminfo % comm, sendListPtr % reqID, mpi_ierr)
         end if
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(recvListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            call mpas_unpack_recv_buf1d_real(dim1, array, recvListPtr, 1, recvListPtr % nlist, recvListPtr % rbuffer, nUnpacked, lastUnpackedIdx)
            deallocate(recvListPtr % rbuffer)
         end if
         recvListPtr => recvListPtr % next
      end do
      
      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(sendListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            deallocate(sendListPtr % rbuffer)
         end if
         sendListPtr => sendListPtr % next
      end do

#endif

   end subroutine mpas_dmpar_exch_halo_field1d_real


   subroutine mpas_dmpar_exch_halo_field2d_real(dminfo, array, dim1, dim2, sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: dim1, dim2
      real (kind=RKIND), dimension(dim1,*), intent(inout) :: array
      type (exchange_list), pointer :: sendList, recvList

      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
      integer :: mpi_ierr
      integer :: d2

#ifdef _MPI

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            d2 = dim1 * recvListPtr % nlist
            allocate(recvListPtr % rbuffer(d2))
            call MPI_Irecv(recvListPtr % rbuffer, d2, MPI_REALKIND, &
                           recvListPtr % procID, recvListPtr % procID, dminfo % comm, recvListPtr % reqID, mpi_ierr)
         end if
         recvListPtr => recvListPtr % next
      end do
      
      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            d2 = dim1 * sendListPtr % nlist
            allocate(sendListPtr % rbuffer(d2))
            call mpas_pack_send_buf2d_real(1, dim1, dim2, array, sendListPtr, 1, d2, sendListPtr % rbuffer, nPacked, lastPackedIdx)
            call MPI_Isend(sendListPtr % rbuffer, d2, MPI_REALKIND, &
                           sendListPtr % procID, dminfo % my_proc_id, dminfo % comm, sendListPtr % reqID, mpi_ierr)
         end if
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(recvListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            d2 = dim1 * recvListPtr % nlist
            call mpas_unpack_recv_buf2d_real(1, dim1, dim2, array, recvListPtr, 1, d2, recvListPtr % rbuffer, nUnpacked, lastUnpackedIdx)
            deallocate(recvListPtr % rbuffer)
         end if
         recvListPtr => recvListPtr % next
      end do
      
      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(sendListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            deallocate(sendListPtr % rbuffer)
         end if
         sendListPtr => sendListPtr % next
      end do

#endif

   end subroutine mpas_dmpar_exch_halo_field2d_real


   subroutine mpas_dmpar_exch_halo_field3d_real(dminfo, array, dim1, dim2, dim3, sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: dim1, dim2, dim3
      real (kind=RKIND), dimension(dim1,dim2,*), intent(inout) :: array
      type (exchange_list), pointer :: sendList, recvList

      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
      integer :: mpi_ierr
      integer :: d3

#ifdef _MPI

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            d3 = dim1 * dim2 * recvListPtr % nlist
            allocate(recvListPtr % rbuffer(d3))
            call MPI_Irecv(recvListPtr % rbuffer, d3, MPI_REALKIND, &
                           recvListPtr % procID, recvListPtr % procID, dminfo % comm, recvListPtr % reqID, mpi_ierr)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            d3 = dim1 * dim2 * sendListPtr % nlist
            allocate(sendListPtr % rbuffer(d3))
            call mpas_pack_send_buf3d_real(1, dim1, 1, dim2, dim3, array, sendListPtr, 1, d3, &
                                   sendListPtr % rbuffer, nPacked, lastPackedIdx)
            call MPI_Isend(sendListPtr % rbuffer, d3, MPI_REALKIND, &
                           sendListPtr % procID, dminfo % my_proc_id, dminfo % comm, sendListPtr % reqID, mpi_ierr)
         end if
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(recvListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            d3 = dim1 * dim2 * recvListPtr % nlist
            call mpas_unpack_recv_buf3d_real(1, dim1, 1, dim2, dim3, array, recvListPtr, 1, d3, &
                                     recvListPtr % rbuffer, nUnpacked, lastUnpackedIdx)
            deallocate(recvListPtr % rbuffer)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(sendListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            deallocate(sendListPtr % rbuffer)
         end if
         sendListPtr => sendListPtr % next
      end do

#endif

   end subroutine mpas_dmpar_exch_halo_field3d_real


end module mpas_dmpar
