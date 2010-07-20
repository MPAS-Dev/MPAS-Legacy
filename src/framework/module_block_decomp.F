module block_decomp

   use dmpar
   use hash

   type graph
      integer :: nVerticesTotal
      integer :: nVertices, maxDegree
      integer :: ghostStart
      integer, dimension(:), pointer :: vertexID
      integer, dimension(:), pointer :: nAdjacent
      integer, dimension(:,:), pointer :: adjacencyList
   end type graph


   contains


   subroutine block_decomp_cells_for_proc(dminfo, partial_global_graph_info, local_cell_list)

      implicit none

      type (dm_info), intent(in) :: dminfo
      type (graph), intent(in) :: partial_global_graph_info
      integer, dimension(:), pointer :: local_cell_list

      integer, dimension(:), pointer :: global_cell_list
      integer, dimension(:), pointer :: global_start

      integer :: i, j, owner, iunit, istatus
      integer, dimension(:), pointer :: local_nvertices
      character (len=256) :: filename

      if (dminfo % nprocs > 1) then

         allocate(local_nvertices(dminfo % nprocs))
         allocate(global_start(dminfo % nprocs))
         allocate(global_cell_list(partial_global_graph_info % nVerticesTotal))

         if (dminfo % my_proc_id == IO_NODE) then

            iunit = 50 + dminfo % my_proc_id
            if (dminfo % nprocs < 10) then
               write(filename,'(a,i1)') 'graph.info.part.', dminfo % nprocs
            else if (dminfo % nprocs < 100) then
               write(filename,'(a,i2)') 'graph.info.part.', dminfo % nprocs
            else if (dminfo % nprocs < 1000) then
               write(filename,'(a,i3)') 'graph.info.part.', dminfo % nprocs
            else if (dminfo % nprocs < 10000) then
               write(filename,'(a,i4)') 'graph.info.part.', dminfo % nprocs
            else if (dminfo % nprocs < 100000) then
               write(filename,'(a,i5)') 'graph.info.part.', dminfo % nprocs
            end if
          
            open(unit=iunit, file=trim(filename), form='formatted', status='old', iostat=istatus)
      
            if (istatus /= 0) then
               write(0,*) 'Could not open block decomposition file for ',dminfo % nprocs,' tasks.'
               write(0,*) 'Filename: ',trim(filename)
               call dmpar_abort(dminfo)
            end if
      
            local_nvertices(:) = 0
            do i=1,partial_global_graph_info % nVerticesTotal
               read(unit=iunit, fmt=*) owner
               local_nvertices(owner+1) = local_nvertices(owner+1) + 1
            end do
      
!            allocate(global_cell_list(partial_global_graph_info % nVerticesTotal))

            global_start(1) = 1
            do i=2,dminfo % nprocs
               global_start(i) = global_start(i-1) + local_nvertices(i-1)
            end do
      
            rewind(unit=iunit)
      
            do i=1,partial_global_graph_info % nVerticesTotal
               read(unit=iunit, fmt=*) owner
               global_cell_list(global_start(owner+1)) = i
               global_start(owner+1) = global_start(owner+1) + 1
            end do

            global_start(1) = 0
            do i=2,dminfo % nprocs
               global_start(i) = global_start(i-1) + local_nvertices(i-1)
            end do

            close(unit=iunit)

            call dmpar_bcast_ints(dminfo, dminfo % nprocs, local_nvertices)
            allocate(local_cell_list(local_nvertices(dminfo % my_proc_id + 1)))

            call dmpar_scatter_ints(dminfo, dminfo % nprocs, local_nvertices(dminfo % my_proc_id + 1), &
                                    global_start, local_nvertices, global_cell_list, local_cell_list)

         else

            call dmpar_bcast_ints(dminfo, dminfo % nprocs, local_nvertices)
            allocate(local_cell_list(local_nvertices(dminfo % my_proc_id + 1)))

            call dmpar_scatter_ints(dminfo, dminfo % nprocs, local_nvertices(dminfo % my_proc_id + 1), &
                                    global_start, local_nvertices, global_cell_list, local_cell_list)

         end if

         deallocate(local_nvertices)
         deallocate(global_start)
         deallocate(global_cell_list)
      else
         allocate(local_cell_list(partial_global_graph_info % nVerticesTotal))
         do i=1,size(local_cell_list)
            local_cell_list(i) = i
         end do
      endif

   end subroutine block_decomp_cells_for_proc


   subroutine block_decomp_partitioned_edge_list(nCells, cellIDList, maxCells, nEdges, cellsOnEdge, edgeIDList, ghostEdgeStart)

      implicit none

      integer, intent(in) :: nCells, maxCells, nEdges
      integer, dimension(nCells), intent(in) :: cellIDList
      integer, dimension(maxCells, nEdges), intent(in) :: cellsOnEdge
      integer, dimension(nEdges), intent(inout) :: edgeIDList
      integer, intent(inout) :: ghostEdgeStart

      integer :: i, j, lastEdge
      integer, dimension(nEdges) :: edgeIDListLocal
      type (hashtable) :: h

      call hash_init(h)

      do i=1,nCells
         ! OPTIMIZATION: Actually, if we can assume that all cellIDs are unique, the if-test is unnecessary
         if (.not. hash_search(h, cellIDList(i))) call hash_insert(h, cellIDList(i))
      end do

      lastEdge = 0
      ghostEdgeStart = nEdges+1

      edgeIDListLocal(:) = edgeIDList(:)

      do i=1,nEdges
         do j=1,maxCells
            if (cellsOnEdge(j,i) /= 0) exit
         end do
if (j > maxCells) write(0,*) 'Error in block_decomp_partitioned_edge_list: edge/vertex is not adjacent to any valid cells'
         if (hash_search(h, cellsOnEdge(j,i))) then
            lastEdge = lastEdge + 1
            edgeIDList(lastEdge) = edgeIDListLocal(i)
         else
            ghostEdgeStart = ghostEdgeStart - 1
            edgeIDList(ghostEdgeStart) = edgeIDListLocal(i)
         end if
if (ghostEdgeStart <= lastEdge) then
   write(0,*) 'block_decomp_partitioned_edge_list: Somehow we have more edges than we thought we should.'
end if
      end do

if (ghostEdgeStart /= lastEdge + 1) then
   write(0,*) 'block_decomp_partitioned_edge_list: Somehow we didn''t have enough edges to fill edgeIDList.'
end if

      call hash_destroy(h)

   end subroutine block_decomp_partitioned_edge_list


   subroutine block_decomp_all_edges_in_block(maxEdges, nCells, nEdgesOnCell, edgesOnCell, nEdges, edgeList)

      implicit none

      integer, intent(in) :: maxEdges, nCells
      integer, dimension(nCells), intent(in) :: nEdgesOnCell
      integer, dimension(maxEdges, nCells), intent(in) :: edgesOnCell
      integer, intent(out) :: nEdges
      integer, dimension(:), pointer :: edgeList

      integer :: i, j, k
      type (hashtable) :: h

      call hash_init(h)

      do i=1,nCells
         do j=1,nEdgesOnCell(i)
            if (.not. hash_search(h, edgesOnCell(j,i))) call hash_insert(h, edgesOnCell(j,i)) 
         end do
      end do

      nEdges = hash_size(h)
      allocate(edgeList(nEdges))

      call hash_destroy(h)

      call hash_init(h)

      k = 0
      do i=1,nCells
         do j=1,nEdgesOnCell(i)
            if (.not. hash_search(h, edgesOnCell(j,i))) then
               k = k + 1
if (k > nEdges) then
   write(0,*) 'block_decomp_all_edges_in_block: Trying to add more edges than expected.'
   return
end if
               edgeList(k) = edgesOnCell(j,i)
               call hash_insert(h, edgesOnCell(j,i)) 
            end if
         end do
      end do

      call hash_destroy(h)

if (k < nEdges) then
   write(0,*) 'block_decomp_all_edges_in_block: Listed fewer edges than expected.'
end if

   end subroutine block_decomp_all_edges_in_block


   subroutine block_decomp_add_halo(dminfo, local_graph_info, local_graph_with_halo)

      implicit none

      type (dm_info), intent(in) :: dminfo
      type (graph), intent(in) :: local_graph_info
      type (graph), intent(out) :: local_graph_with_halo

      integer :: i, j, k
      type (hashtable) :: h


      call hash_init(h)

      do i=1,local_graph_info % nVertices
         do j=1,local_graph_info % nAdjacent(i)
            if (local_graph_info % adjacencyList(j,i) /= 0) then
               if (.not. hash_search(h, local_graph_info % adjacencyList(j,i))) then
                  call hash_insert(h, local_graph_info % adjacencyList(j,i))
               end if
            end if
         end do
      end do 


      local_graph_with_halo % nVertices = local_graph_info % nVertices
      local_graph_with_halo % maxDegree = local_graph_info % maxDegree
      local_graph_with_halo % nVerticesTotal = hash_size(h)
      local_graph_with_halo % ghostStart = local_graph_with_halo % nVertices + 1
      allocate(local_graph_with_halo % vertexID(local_graph_with_halo % nVerticesTotal))
      allocate(local_graph_with_halo % nAdjacent(local_graph_with_halo % nVerticesTotal))
      allocate(local_graph_with_halo % adjacencyList(local_graph_with_halo % maxDegree, local_graph_with_halo % nVerticesTotal))

      call hash_destroy(h)

      call hash_init(h)

      do i=1,local_graph_info % nVertices
if (hash_search(h, local_graph_info % vertexID(i))) &
  write(0,*) 'block_decomp_add_halo: There appear to be duplicates in vertexID list.'
         call hash_insert(h, local_graph_info % vertexID(i)) 
         local_graph_with_halo % vertexID(i) = local_graph_info % vertexID(i) 
         local_graph_with_halo % nAdjacent(i) = local_graph_info % nAdjacent(i) 
         local_graph_with_halo % adjacencyList(:,i) = local_graph_info % adjacencyList(:,i) 
      end do

      k = local_graph_with_halo % ghostStart
if (hash_size(h) /= k-1) &
  write(0,*) 'block_decomp_add_halo: Somehow we don''t have the right number of non-ghost cells.'
      do i=1,local_graph_info % nVertices
         do j=1,local_graph_info % nAdjacent(i)
            if (local_graph_info % adjacencyList(j,i) /= 0) then
               if (.not. hash_search(h, local_graph_info % adjacencyList(j,i))) then
                  call hash_insert(h, local_graph_info % adjacencyList(j,i))
                  local_graph_with_halo % vertexID(k) = local_graph_info % adjacencyList(j,i)
                  k = k + 1
               end if
            end if
         end do
      end do 
if (local_graph_with_halo % nVerticesTotal /= k-1) &
  write(0,*) 'block_decomp_add_halo: Somehow we don''t have the right number of total cells.'

      call hash_destroy(h)

   end subroutine block_decomp_add_halo

end module block_decomp
