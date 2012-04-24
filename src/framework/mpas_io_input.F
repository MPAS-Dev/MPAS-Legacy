module mpas_io_input

   use mpas_grid_types
   use mpas_dmpar
   use mpas_block_decomp
   use mpas_sort
   use mpas_configure
   use mpas_timekeeping
   use mpas_io_streams


#ifdef HAVE_ZOLTAN
   use mpas_zoltan_interface
#endif

   integer, parameter :: STREAM_INPUT=1, STREAM_SFC=2, STREAM_RESTART=3

   type io_input_object
      character (len=1024) :: filename
      integer :: rd_ncid
      integer :: stream

      integer :: time

      type (MPAS_Stream_type) :: io_stream

   end type io_input_object


   type (exchange_list), pointer :: sendCellList, recvCellList
   type (exchange_list), pointer :: sendEdgeList, recvEdgeList
   type (exchange_list), pointer :: sendVertexList, recvVertexList
   type (exchange_list), pointer :: sendVertLevelList, recvVertLevelList
 
   integer :: readCellStart, readCellEnd, nReadCells
   integer :: readEdgeStart, readEdgeEnd, nReadEdges
   integer :: readVertexStart, readVertexEnd, nReadVertices
   integer :: readVertLevelStart, readVertLevelEnd, nReadVertLevels
   

   contains


   subroutine mpas_input_state_for_domain(domain)
   
      implicit none
   
      type (domain_type), pointer :: domain
   
      integer :: i, j, k
      type (io_input_object) :: input_obj
#include "dim_decls.inc"

      character (len=16) :: c_on_a_sphere
      real (kind=RKIND) :: r_sphere_radius

      integer :: ierr
      integer, dimension(:), pointer :: readIndices
      type (MPAS_IO_Handle_type) :: inputHandle
   
      type (field1dInteger) :: indexToCellIDField
      type (field1dInteger) :: indexToEdgeIDField
      type (field1dInteger) :: indexToVertexIDField
      type (field1dInteger) :: nEdgesOnCellField
      type (field2dInteger) :: cellsOnCellField
      type (field2dInteger) :: edgesOnCellField
      type (field2dInteger) :: verticesOnCellField
      type (field2dInteger) :: cellsOnEdgeField
      type (field2dInteger) :: cellsOnVertexField

#ifdef HAVE_ZOLTAN
#ifdef _MPI 
      type (field1dReal) :: xCellField,   yCellField,   zCellField
      type (field1dReal) :: xEdgeField,   yEdgeField,   zEdgeField
      type (field1dReal) :: xVertexField, yVertexField, zVertexField
#endif
#endif

      type (field1DChar) :: xtime
   
      integer, dimension(:),   pointer :: indexToCellID_0Halo
      integer, dimension(:),   pointer :: nEdgesOnCell_0Halo
      integer, dimension(:,:), pointer :: cellsOnCell_0Halo

      integer, dimension(:),   pointer :: nEdgesOnCell_2Halo

      integer, dimension(:,:), pointer :: edgesOnCell_2Halo
      integer, dimension(:,:), pointer :: verticesOnCell_2Halo

      integer, dimension(:,:), pointer :: cellsOnEdge_2Halo
      integer, dimension(:,:), pointer :: cellsOnVertex_2Halo

      integer, dimension(:,:), pointer :: cellIDSorted
      integer, dimension(:,:), pointer :: edgeIDSorted
      integer, dimension(:,:), pointer :: vertexIDSorted

#ifdef HAVE_ZOLTAN
#ifdef _MPI 
      real (kind=RKIND), dimension(:), pointer :: xCell,   yCell,   zCell
      real (kind=RKIND), dimension(:), pointer :: xEdge,   yEdge,   zEdge
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
#endif
#endif
   
      integer, dimension(:), pointer :: local_cell_list, local_edge_list, local_vertex_list
      integer, dimension(:), pointer :: block_id, block_start, block_count
      integer, dimension(:), pointer :: local_vertlevel_list, needed_vertlevel_list
      integer :: nlocal_edges, nlocal_vertices
      type (exchange_list), pointer :: send1Halo, recv1Halo
      type (exchange_list), pointer :: send2Halo, recv2Halo
      type (graph) :: partial_global_graph_info
      type (graph) :: block_graph_0Halo, block_graph_1Halo, block_graph_2Halo
      integer :: ghostEdgeStart, ghostVertexStart

      type (MPAS_Time_type) :: startTime
      type (MPAS_Time_type) :: sliceTime
      type (MPAS_TimeInterval_type) :: timeDiff
      type (MPAS_TimeInterval_type) :: minTimeDiff
      character(len=32) :: timeStamp
      character(len=1024) :: filename

      integer, parameter :: nHalos = 2
      integer, dimension(nHalos+1) :: nCellsCumulative    ! own cells, halo 1 cells, halo 2 cells
      integer, dimension(nHalos+2) :: nEdgesCumulative    ! own edges, own cell's edges, halo 1 edges, halo 2 edges
      integer, dimension(nHalos+2) :: nVerticesCumulative ! own vertices, own cell's vertices, halo 1 vertices, halo 2 vertices

      integer, dimension(nHalos)   :: nCellsHalo          ! halo 1 cells, halo 2 cells
      integer, dimension(nHalos+1) :: nEdgesHalo          ! own cell's edges, halo 1 edges, halo 2 edges
      integer, dimension(nHalos+1) :: nVerticesHalo       ! own cell's vertices, halo 1 vertices, halo 2 vertices

      integer, dimension(:), pointer :: tempIDs
      integer :: ntempIDs, offset

      integer :: nHalo, nOwnCells, nOwnEdges, nOwnVertices, cellCount, edgeCount, vertexCount, iEdge, iVertex
      type (hashtable) :: edgeHash, vertexHash


      if (config_do_restart) then

         ! this get followed by set is to ensure that the time is in standard format
         call mpas_set_time(curr_time=startTime, dateTimeString=config_start_time)
         call mpas_get_time(curr_time=startTime, dateTimeString=timeStamp)

         call mpas_insert_string_suffix(trim(config_restart_name), timeStamp, filename)

         input_obj % filename = trim(filename)
         input_obj % stream = STREAM_RESTART
      else
         input_obj % filename = trim(config_input_name)
         input_obj % stream = STREAM_INPUT
      end if
      inputHandle = MPAS_io_open(trim(input_obj % filename), MPAS_IO_READ, MPAS_IO_PNETCDF, ierr)
      if (ierr /= MPAS_IO_NOERR) then
         write(0,*) ' '
         if (input_obj % stream == STREAM_RESTART) then
            write(0,*) 'Error opening restart file ''', trim(input_obj % filename), ''''
         else if (input_obj % stream == STREAM_INPUT) then
            write(0,*) 'Error opening input file ''', trim(input_obj % filename), ''''
         else if (input_obj % stream == STREAM_SFC) then
            write(0,*) 'Error opening sfc file ''', trim(input_obj % filename), ''''
         end if
         write(0,*) ' '
         call mpas_dmpar_abort(domain % dminfo)
      end if
   

      !
      ! Read global number of cells/edges/vertices
      !
#include "read_dims.inc"
   
      !
      ! Determine the range of cells/edges/vertices that a processor will initially read
      !   from the input file
      !
      call mpas_dmpar_get_index_range(domain % dminfo, 1, nCells, readCellStart, readCellEnd)   
      nReadCells = readCellEnd - readCellStart + 1
   
      call mpas_dmpar_get_index_range(domain % dminfo, 1, nEdges, readEdgeStart, readEdgeEnd)   
      nReadEdges = readEdgeEnd - readEdgeStart + 1
   
      call mpas_dmpar_get_index_range(domain % dminfo, 1, nVertices, readVertexStart, readVertexEnd)   
      nReadVertices = readVertexEnd - readVertexStart + 1

      readVertLevelStart = 1
      readVertLevelEnd = nVertLevels
      nReadVertLevels = nVertLevels
   
   
      !
      ! Allocate and read fields that we will need in order to ultimately work out
      !   which cells/edges/vertices are owned by each block, and which are ghost
      !

      ! Global cell indices
      allocate(indexToCellIDField % ioinfo)
      indexToCellIDField % ioinfo % fieldName = 'indexToCellID'
      indexToCellIDField % ioinfo % start(1) = readCellStart
      indexToCellIDField % ioinfo % count(1) = nReadCells
      allocate(indexToCellIDField % array(nReadCells))
      allocate(readIndices(nReadCells))
      do i=1,nReadCells
         readIndices(i) = i + readCellStart - 1
      end do
      call MPAS_io_inq_var(inputHandle, 'indexToCellID', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'indexToCellID', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'indexToCellID', indexToCellIDField % array, ierr)
   
#ifdef HAVE_ZOLTAN
#ifdef _MPI 
      ! Cell x-coordinates (in 3d Cartesian space)
      allocate(xCellField % ioinfo)
      xCellField % ioinfo % fieldName = 'xCell'
      xCellField % ioinfo % start(1) = readCellStart
      xCellField % ioinfo % count(1) = nReadCells
      allocate(xCellField % array(nReadCells))
      call MPAS_io_inq_var(inputHandle, 'xCell', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'xCell', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'xCell', xCellField % array, ierr)

      ! Cell y-coordinates (in 3d Cartesian space)
      allocate(yCellField % ioinfo)
      yCellField % ioinfo % fieldName = 'yCell'
      yCellField % ioinfo % start(1) = readCellStart
      yCellField % ioinfo % count(1) = nReadCells
      allocate(yCellField % array(nReadCells))
      call MPAS_io_inq_var(inputHandle, 'yCell', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'yCell', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'yCell', yCellField % array, ierr)

      ! Cell z-coordinates (in 3d Cartesian space)
      allocate(zCellField % ioinfo)
      zCellField % ioinfo % fieldName = 'zCell'
      zCellField % ioinfo % start(1) = readCellStart
      zCellField % ioinfo % count(1) = nReadCells
      allocate(zCellField % array(nReadCells))
      call MPAS_io_inq_var(inputHandle, 'zCell', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'zCell', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'zCell', zCellField % array, ierr)
#endif
#endif
      deallocate(readIndices)


      ! Global edge indices
      allocate(indexToEdgeIDField % ioinfo)
      indexToEdgeIDField % ioinfo % fieldName = 'indexToEdgeID'
      indexToEdgeIDField % ioinfo % start(1) = readEdgeStart
      indexToEdgeIDField % ioinfo % count(1) = nReadEdges
      allocate(indexToEdgeIDField % array(nReadEdges))
      allocate(indexToEdgeIDField % array(nReadEdges))
      allocate(readIndices(nReadEdges))
      do i=1,nReadEdges
         readIndices(i) = i + readEdgeStart - 1
      end do
      call MPAS_io_inq_var(inputHandle, 'indexToEdgeID', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'indexToEdgeID', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'indexToEdgeID', indexToEdgeIDField % array, ierr)
   
#ifdef HAVE_ZOLTAN
#ifdef _MPI 
      ! Edge x-coordinates (in 3d Cartesian space)
      allocate(xEdgeField % ioinfo)
      xEdgeField % ioinfo % fieldName = 'xEdge'
      xEdgeField % ioinfo % start(1) = readEdgeStart
      xEdgeField % ioinfo % count(1) = nReadEdges
      allocate(xEdgeField % array(nReadEdges))
      call MPAS_io_inq_var(inputHandle, 'xEdge', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'xEdge', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'xEdge', xEdgeField % array, ierr)

      ! Edge y-coordinates (in 3d Cartesian space)
      allocate(yEdgeField % ioinfo)
      yEdgeField % ioinfo % fieldName = 'yEdge'
      yEdgeField % ioinfo % start(1) = readEdgeStart
      yEdgeField % ioinfo % count(1) = nReadEdges
      allocate(yEdgeField % array(nReadEdges))
      call MPAS_io_inq_var(inputHandle, 'yEdge', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'yEdge', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'yEdge', yEdgeField % array, ierr)

      ! Edge z-coordinates (in 3d Cartesian space)
      allocate(zEdgeField % ioinfo)
      zEdgeField % ioinfo % fieldName = 'zEdge'
      zEdgeField % ioinfo % start(1) = readEdgeStart
      zEdgeField % ioinfo % count(1) = nReadEdges
      allocate(zEdgeField % array(nReadEdges))
      call MPAS_io_inq_var(inputHandle, 'zEdge', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'zEdge', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'zEdge', zEdgeField % array, ierr)
#endif
#endif
      deallocate(readIndices)


      ! Global vertex indices
      allocate(indexToVertexIDField % ioinfo)
      indexToVertexIDField % ioinfo % fieldName = 'indexToVertexID'
      indexToVertexIDField % ioinfo % start(1) = readVertexStart
      indexToVertexIDField % ioinfo % count(1) = nReadVertices
      allocate(indexToVertexIDField % array(nReadVertices))
      allocate(readIndices(nReadVertices))
      do i=1,nReadVertices
         readIndices(i) = i + readVertexStart - 1
      end do
      call MPAS_io_inq_var(inputHandle, 'indexToVertexID', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'indexToVertexID', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'indexToVertexID', indexToVertexIDField % array, ierr)
   
#ifdef HAVE_ZOLTAN
#ifdef _MPI
      ! Vertex x-coordinates (in 3d Cartesian space)
      allocate(xVertexField % ioinfo)
      xVertexField % ioinfo % fieldName = 'xVertex'
      xVertexField % ioinfo % start(1) = readVertexStart
      xVertexField % ioinfo % count(1) = nReadVertices
      allocate(xVertexField % array(nReadVertices))
      call MPAS_io_inq_var(inputHandle, 'xVertex', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'xVertex', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'xVertex', xVertexField % array, ierr)

      ! Vertex y-coordinates (in 3d Cartesian space)
      allocate(yVertexField % ioinfo)
      yVertexField % ioinfo % fieldName = 'yVertex'
      yVertexField % ioinfo % start(1) = readVertexStart
      yVertexField % ioinfo % count(1) = nReadVertices
      allocate(yVertexField % array(nReadVertices))
      call MPAS_io_inq_var(inputHandle, 'yVertex', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'yVertex', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'yVertex', yVertexField % array, ierr)

      ! Vertex z-coordinates (in 3d Cartesian space)
      allocate(zVertexField % ioinfo)
      zVertexField % ioinfo % fieldName = 'zVertex'
      zVertexField % ioinfo % start(1) = readVertexStart
      zVertexField % ioinfo % count(1) = nReadVertices
      allocate(zVertexField % array(nReadVertices))
      call MPAS_io_inq_var(inputHandle, 'zVertex', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'zVertex', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'zVertex', zVertexField % array, ierr)
#endif
#endif
      deallocate(readIndices)

      ! Number of cell/edges/vertices adjacent to each cell
      allocate(nEdgesOnCellField % ioinfo)
      nEdgesOnCellField % ioinfo % fieldName = 'nEdgesOnCell'
      nEdgesOnCellField % ioinfo % start(1) = readCellStart
      nEdgesOnCellField % ioinfo % count(1) = nReadCells
      allocate(nEdgesOnCellField % array(nReadCells))
      allocate(readIndices(nReadCells))
      do i=1,nReadCells
         readIndices(i) = i + readCellStart - 1
      end do
      call MPAS_io_inq_var(inputHandle, 'nEdgesOnCell', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'nEdgesOnCell', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'nEdgesOnCell', nEdgesOnCellField % array, ierr)
   
      ! Global indices of cells adjacent to each cell
      allocate(cellsOnCellField % ioinfo)
      cellsOnCellField % ioinfo % fieldName = 'cellsOnCell'
      cellsOnCellField % ioinfo % start(1) = 1
      cellsOnCellField % ioinfo % start(2) = readCellStart
      cellsOnCellField % ioinfo % count(1) = maxEdges
      cellsOnCellField % ioinfo % count(2) = nReadCells
      allocate(cellsOnCellField % array(maxEdges,nReadCells))
      call MPAS_io_inq_var(inputHandle, 'cellsOnCell', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'cellsOnCell', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'cellsOnCell', cellsOnCellField % array, ierr)
   
      ! Global indices of edges adjacent to each cell
      allocate(edgesOnCellField % ioinfo)
      edgesOnCellField % ioinfo % fieldName = 'edgesOnCell'
      edgesOnCellField % ioinfo % start(1) = 1
      edgesOnCellField % ioinfo % start(2) = readCellStart
      edgesOnCellField % ioinfo % count(1) = maxEdges
      edgesOnCellField % ioinfo % count(2) = nReadCells
      allocate(edgesOnCellField % array(maxEdges,nReadCells))
      call MPAS_io_inq_var(inputHandle, 'edgesOnCell', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'edgesOnCell', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'edgesOnCell', edgesOnCellField % array, ierr)
   
      ! Global indices of vertices adjacent to each cell
      allocate(verticesOnCellField % ioinfo)
      verticesOnCellField % ioinfo % fieldName = 'verticesOnCell'
      verticesOnCellField % ioinfo % start(1) = 1
      verticesOnCellField % ioinfo % start(2) = readCellStart
      verticesOnCellField % ioinfo % count(1) = maxEdges
      verticesOnCellField % ioinfo % count(2) = nReadCells
      allocate(verticesOnCellField % array(maxEdges,nReadCells))
      call MPAS_io_inq_var(inputHandle, 'verticesOnCell', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'verticesOnCell', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'verticesOnCell', verticesOnCellField % array, ierr)
      deallocate(readIndices)
   
      ! Global indices of cells adjacent to each edge
      !    used for determining which edges are owned by a block, where 
      !    iEdge is owned iff cellsOnEdge(1,iEdge) is an owned cell
      allocate(cellsOnEdgeField % ioinfo)
      cellsOnEdgeField % ioinfo % fieldName = 'cellsOnEdge'
      cellsOnEdgeField % ioinfo % start(1) = 1
      cellsOnEdgeField % ioinfo % start(2) = readEdgeStart
      cellsOnEdgeField % ioinfo % count(1) = 2
      cellsOnEdgeField % ioinfo % count(2) = nReadEdges
      allocate(cellsOnEdgeField % array(2,nReadEdges))
      allocate(readIndices(nReadEdges))
      do i=1,nReadEdges
         readIndices(i) = i + readEdgeStart - 1
      end do
      call MPAS_io_inq_var(inputHandle, 'cellsOnEdge', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'cellsOnEdge', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'cellsOnEdge', cellsOnEdgeField % array, ierr)
      deallocate(readIndices)
   
      ! Global indices of cells adjacent to each vertex
      !    used for determining which vertices are owned by a block, where 
      !    iVtx is owned iff cellsOnVertex(1,iVtx) is an owned cell
      allocate(cellsOnVertexField % ioinfo)
      cellsOnVertexField % ioinfo % fieldName = 'cellsOnVertex'
      cellsOnVertexField % ioinfo % start(1) = 1
      cellsOnVertexField % ioinfo % start(2) = readVertexStart
      cellsOnVertexField % ioinfo % count(1) = vertexDegree
      cellsOnVertexField % ioinfo % count(2) = nReadVertices
      allocate(cellsOnVertexField % array(vertexDegree,nReadVertices))
      allocate(readIndices(nReadVertices))
      do i=1,nReadVertices
         readIndices(i) = i + readVertexStart - 1
      end do
      call MPAS_io_inq_var(inputHandle, 'cellsOnVertex', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'cellsOnVertex', readIndices, ierr=ierr)
      call mpas_io_get_var(inputHandle, 'cellsOnVertex', cellsOnVertexField % array, ierr)
      deallocate(readIndices)
   
   
      !
      ! Set up a graph derived data type describing the connectivity for the cells 
      !   that were read by this process
      ! A partial description is passed to the block decomp module by each process,
      !   and the block decomp module returns with a list of global cell indices
      !   that belong to the block on this process
      !
      partial_global_graph_info % nVertices = nReadCells
      partial_global_graph_info % nVerticesTotal = nCells
      partial_global_graph_info % maxDegree = maxEdges
      partial_global_graph_info % ghostStart = nVertices+1
      allocate(partial_global_graph_info % vertexID(nReadCells))
      allocate(partial_global_graph_info % nAdjacent(nReadCells))
      allocate(partial_global_graph_info % adjacencyList(maxEdges, nReadCells))
   
      partial_global_graph_info % vertexID(:) = indexToCellIDField % array(:)
      partial_global_graph_info % nAdjacent(:) = nEdgesOnCellField % array(:)
      partial_global_graph_info % adjacencyList(:,:) = cellsOnCellField % array(:,:)
      
   
      ! TODO: Ensure (by renaming or exchanging) that initial cell range on each proc is contiguous
      !       This situation may occur when reading a restart file with cells/edges/vertices written
      !       in a scrambled order
   

      ! Determine which cells are owned by this process
      call mpas_block_decomp_cells_for_proc(domain % dminfo, partial_global_graph_info, local_cell_list, block_id, block_start, block_count)

      deallocate(partial_global_graph_info % vertexID)
      deallocate(partial_global_graph_info % nAdjacent)
      deallocate(partial_global_graph_info % adjacencyList)
   
   
      allocate(indexToCellID_0Halo(size(local_cell_list)))
      allocate(nEdgesOnCell_0Halo(size(local_cell_list)))
      allocate(cellsOnCell_0Halo(maxEdges, size(local_cell_list)))
   
#ifdef HAVE_ZOLTAN
#ifdef _MPI 
      allocate(xCell(size(local_cell_list)))
      allocate(yCell(size(local_cell_list)))
      allocate(zCell(size(local_cell_list)))
#endif
#endif
   
      !
      ! Now that each process has a list of cells that it owns, exchange cell connectivity 
      !   information between the processes that read info for a cell and those that own that cell
      !
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                size(indexToCellIDField % array), size(local_cell_list), &
                                indexToCellIDField % array, local_cell_list, &
                                sendCellList, recvCellList)
   
      call mpas_dmpar_alltoall_field(domain % dminfo, indexToCellIDField % array, indexToCellID_0Halo, &
                                size(indexToCellIDField % array), size(local_cell_list), &
                                sendCellList, recvCellList)
   
      call mpas_dmpar_alltoall_field(domain % dminfo, nEdgesOnCellField % array, nEdgesOnCell_0Halo, &
                                size(indexToCellIDField % array), size(local_cell_list), &
                                sendCellList, recvCellList)
   
      call mpas_dmpar_alltoall_field(domain % dminfo, cellsOnCellField % array, cellsOnCell_0Halo, &
                                size(cellsOnCellField % array, 1), size(indexToCellIDField % array), size(local_cell_list), &
                                sendCellList, recvCellList)
   
#ifdef HAVE_ZOLTAN
#ifdef _MPI 
      call mpas_dmpar_alltoall_field(domain % dminfo, xCellField % array, xCell, &
                                size(xCellField % array), size(local_cell_list), &
                                sendCellList, recvCellList)
   
      call mpas_dmpar_alltoall_field(domain % dminfo, yCellField % array, yCell, &
                                size(yCellField % array), size(local_cell_list), &
                                sendCellList, recvCellList)
   
      call mpas_dmpar_alltoall_field(domain % dminfo, zCellField % array, zCell, &
                                size(zCellField % array), size(local_cell_list), &
                                sendCellList, recvCellList)
#endif
#endif


      deallocate(sendCellList % list)
      deallocate(sendCellList)
      deallocate(recvCellList % list)
      deallocate(recvCellList)



      !
      ! Build a graph of cell connectivity based on cells owned by this process
      !
      block_graph_0Halo % nVerticesTotal = size(local_cell_list)
      block_graph_0Halo % nVertices = size(local_cell_list)
      block_graph_0Halo % maxDegree = maxEdges
      block_graph_0Halo % ghostStart = size(local_cell_list) + 1
      allocate(block_graph_0Halo % vertexID(size(local_cell_list)))
      allocate(block_graph_0Halo % nAdjacent(size(local_cell_list)))
      allocate(block_graph_0Halo % adjacencyList(maxEdges, size(local_cell_list)))
   
      block_graph_0Halo % vertexID(:) = indexToCellID_0Halo(:)
      block_graph_0Halo % nAdjacent(:) = nEdgesOnCell_0Halo(:)
      block_graph_0Halo % adjacencyList(:,:) = cellsOnCell_0Halo(:,:)
   
      ! Get back a graph describing the owned cells plus the cells in the 1-halo
      call mpas_block_decomp_add_halo(domain % dminfo, block_graph_0Halo, block_graph_1Halo)
   
   
      !
      ! Work out exchange lists for 1-halo and exchange cell information for 1-halo
      !
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                block_graph_0Halo % nVertices, block_graph_1Halo % nVerticesTotal, &
                                block_graph_0Halo % vertexID,  block_graph_1Halo % vertexID, &
                                send1Halo, recv1Halo)
   
      call mpas_dmpar_alltoall_field(domain % dminfo, block_graph_0Halo % vertexID, block_graph_1Halo % vertexID, &
                                block_graph_0Halo % nVertices, block_graph_1Halo % nVerticesTotal, &
                                send1Halo, recv1Halo)
   
      call mpas_dmpar_alltoall_field(domain % dminfo, block_graph_0Halo % nAdjacent, block_graph_1Halo % nAdjacent, &
                                block_graph_0Halo % nVertices, block_graph_1Halo % nVerticesTotal, &
                                send1Halo, recv1Halo)
   
      call mpas_dmpar_alltoall_field(domain % dminfo, block_graph_0Halo % adjacencyList, block_graph_1Halo % adjacencyList, &
                                block_graph_0Halo % maxDegree, block_graph_0Halo % nVertices, block_graph_1Halo % nVerticesTotal, &
                                send1Halo, recv1Halo)
   
   
      !
      ! Work out exchange lists for 2-halo and exchange cell information for 2-halo
      !
      block_graph_1Halo % nVertices = block_graph_1Halo % nVerticesTotal
      block_graph_1Halo % ghostStart = block_graph_1Halo % nVerticesTotal + 1
     
      ! Get back a graph describing the owned and 1-halo cells plus the cells in the 2-halo
      call mpas_block_decomp_add_halo(domain % dminfo, block_graph_1Halo, block_graph_2Halo)
   
      block_graph_2Halo % nVertices = block_graph_0Halo % nVertices
      block_graph_2Halo % ghostStart = block_graph_2Halo % nVertices + 1

      nOwnCells = block_graph_2Halo % nVertices

#ifdef HAVE_ZOLTAN
#ifdef _MPI 
      !! For now, only use Zoltan with MPI
      !! Zoltan initialization
      call mpas_zoltan_start()

      !! Zoltan hook for cells
      call mpas_zoltan_order_loc_hsfc_cells(block_graph_2Halo%nVertices,block_graph_2Halo%VertexID,3,xCell,yCell,zCell)
#endif
#endif

      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                block_graph_0Halo % nVertices, block_graph_2Halo % nVerticesTotal, &
                                block_graph_0Halo % vertexID,  block_graph_2Halo % vertexID, &
                                send2Halo, recv2Halo)
   
      call mpas_dmpar_alltoall_field(domain % dminfo, block_graph_0Halo % vertexID, block_graph_2Halo % vertexID, &
                                block_graph_0Halo % nVertices, block_graph_2Halo % nVerticesTotal, &
                                send2Halo, recv2Halo)
   
      call mpas_dmpar_alltoall_field(domain % dminfo, block_graph_0Halo % nAdjacent, block_graph_2Halo % nAdjacent, &
                                block_graph_0Halo % nVertices, block_graph_2Halo % nVerticesTotal, &
                                send2Halo, recv2Halo)
   
      call mpas_dmpar_alltoall_field(domain % dminfo, block_graph_0Halo % adjacencyList, block_graph_2Halo % adjacencyList, &
                                block_graph_0Halo % maxDegree, block_graph_0Halo % nVertices, block_graph_2Halo % nVerticesTotal, &
                                send2Halo, recv2Halo)


   
      !
      ! Knowing which cells are in block and the 2-halo, we can exchange lists of which edges are
      !   on each cell and which vertices are on each cell from the processes that read these
      !   fields for each cell to the processes that own the cells
      !
      allocate(nEdgesOnCell_2Halo(block_graph_2Halo % nVerticesTotal))
      allocate(edgesOnCell_2Halo(maxEdges, block_graph_2Halo % nVerticesTotal))
      allocate(verticesOnCell_2Halo(maxEdges, block_graph_2Halo % nVerticesTotal))
   
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                size(indexToCellIDField % array), block_graph_2Halo % nVerticesTotal, &
                                indexToCellIDField % array, block_graph_2Halo % vertexID, &
                                sendCellList, recvCellList)
   
      call mpas_dmpar_alltoall_field(domain % dminfo, nEdgesOnCellField % array, nEdgesOnCell_2Halo, &
                                size(indexToCellIDField % array), size(local_cell_list), &
                                sendCellList, recvCellList)

      call mpas_dmpar_alltoall_field(domain % dminfo, edgesOnCellField % array, edgesOnCell_2Halo, &
                                maxEdges, nReadCells, block_graph_2Halo % nVerticesTotal, &
                                sendCellList, recvCellList)
   
      call mpas_dmpar_alltoall_field(domain % dminfo, verticesOnCellField % array, verticesOnCell_2Halo, &
                                maxEdges, nReadCells, block_graph_2Halo % nVerticesTotal, &
                                sendCellList, recvCellList)

   
      ! 
      ! Get a list of which edges and vertices are adjacent to cells (including halo cells) in block
      ! 
      call mpas_block_decomp_all_edges_in_block(maxEdges, block_graph_2Halo % nVerticesTotal, block_graph_2Halo % nAdjacent, &
                                           edgesOnCell_2Halo, nlocal_edges, local_edge_list)
      call mpas_block_decomp_all_edges_in_block(maxEdges, block_graph_2Halo % nVerticesTotal, block_graph_2Halo % nAdjacent, &
                                           verticesOnCell_2Halo, nlocal_vertices, local_vertex_list)
   
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                size(indexToEdgeIDField % array), nlocal_edges, &
                                indexToEdgeIDField % array, local_edge_list, &
                                sendEdgeList, recvEdgeList)
   
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                size(indexToVertexIDField % array), nlocal_vertices, &
                                indexToVertexIDField % array, local_vertex_list, &
                                sendVertexList, recvVertexList)
   
   
   
      ! 
      ! Work out which edges and vertices are owned by this process, and which are ghost
      ! 
      allocate(cellsOnEdge_2Halo(2,nlocal_edges))
      allocate(cellsOnVertex_2Halo(vertexDegree,nlocal_vertices))
   
      call mpas_dmpar_alltoall_field(domain % dminfo, cellsOnEdgeField % array, cellsOnEdge_2Halo, &
                                2, size(cellsOnEdgeField % array, 2), nlocal_edges, &
                                sendEdgeList, recvEdgeList)
   
      call mpas_dmpar_alltoall_field(domain % dminfo, cellsOnVertexField % array, cellsOnVertex_2Halo, &
                                vertexDegree, size(cellsOnVertexField % array, 2), nlocal_vertices, &
                                sendVertexList, recvVertexList)
   
   
      call mpas_block_decomp_partitioned_edge_list(nOwnCells, &
                                              block_graph_2Halo % vertexID(1:nOwnCells), &
                                              2, nlocal_edges, cellsOnEdge_2Halo, local_edge_list, ghostEdgeStart)

      call mpas_block_decomp_partitioned_edge_list(nOwnCells, &
                                              block_graph_2Halo % vertexID(1:nOwnCells), &
                                              vertexDegree, nlocal_vertices, cellsOnVertex_2Halo, local_vertex_list, ghostVertexStart)

      !------- set owned and halo cell indices -------!  
      
      nCellsCumulative(1) = nOwnCells
      nCellsCumulative(2) = block_graph_1Halo % nVerticesTotal
      nCellsCumulative(3) = block_graph_2Halo % nVerticesTotal

      !------- determin the perimeter and owned edges of own cells and halos -------!  

      nOwnEdges = ghostEdgeStart-1
      nOwnVertices = ghostVertexStart-1

      ! skip the own edges found at the beginning of local_edge_list
      call mpas_hash_init(edgeHash)
      do i=1,nOwnEdges
         call mpas_hash_insert(edgeHash, local_edge_list(i))
      end do

      ! skip the own vertices found at the beginning of local_vertex_list
      call mpas_hash_init(vertexHash)
      do i=1,nOwnVertices
         call mpas_hash_insert(vertexHash, local_vertex_list(i))
      end do

      cellCount = 1              !tracks the index of the local cell array
      edgeCount = nOwnEdges      !tracks where to insert the next local edge
      vertexCount = nOwnVertices !tracks where to insert the next local vertex

      nEdgesCumulative(1) = nOwnEdges
      nVerticesCumulative(1) = nOwnVertices

      !Order the local_edge_list and local_vertex_list accordingly and set the bounds of each perimeter ---- 
      do i = 1, nHalos + 1 ! for the own cells and each halo...
         do j = cellCount, nCellsCumulative(i)

            ! the number of edges on a cell is same to the number of vertices, and therefore
            ! nEdgesOnCell_2Halo(j) will be the correct upper bound for both both edges and vertices on cell
            do k = 1, nEdgesOnCell_2Halo(j)
               iEdge = edgesOnCell_2Halo(k,j)
               if (.not. mpas_hash_search(edgeHash, iEdge)) then
                  edgeCount = edgeCount + 1
                  local_edge_list(edgeCount) = iEdge
                  call mpas_hash_insert(edgeHash, iEdge)
               end if

               iVertex = verticesOnCell_2Halo(k,j)
               if (.not. mpas_hash_search(vertexHash, iVertex)) then
                  vertexCount = vertexCount + 1
                  local_vertex_list(vertexCount) = iVertex
                  call mpas_hash_insert(vertexHash, iVertex)
               end if
            end do

         end do

         cellCount = nCellsCumulative(i) + 1
         nEdgesCumulative(i+1) = edgeCount
         nVerticesCumulative(i+1) = vertexCount
      end do

      do i = 1, nHalos
         nCellsHalo(i) = nCellsCumulative(i+1) - nCellsCumulative(i)
      end do

      do i = 1, nHalos + 1
         nEdgesHalo(i) = nEdgesCumulative(i+1) - nEdgesCumulative(i)
      end do

      do i = 1, nHalos + 1
         nVerticesHalo(i) = nVerticesCumulative(i+1) - nVerticesCumulative(i)
      end do

      call mpas_hash_destroy(edgeHash)
      call mpas_hash_destroy(vertexHash)


      ! At this point, local_edge_list(1:nOwnEdges) contains all of the owned edges for this block
      !   and local_edge_list(ghostEdgeStart:nlocal_edges) contains all of the ghost edges

      ! At this point, local_vertex_list(1;ghostVertexStart-1) contains all of the owned vertices for this block
      !   and local_vertex_list(ghostVertexStart:nlocal_vertices) contains all of the ghost vertices

      ! Also, at this point, block_graph_2Halo % vertexID(1:block_graph_2Halo%nVertices) contains all of the owned
      !   cells for this block, and block_graph_2Halo % vertexID(block_graph_2Halo%nVertices+1:block_graph_2Halo%nVerticesTotal)
      !   contains all of the ghost cells


      deallocate(sendEdgeList % list)
      deallocate(sendEdgeList)
      deallocate(recvEdgeList % list)
      deallocate(recvEdgeList)
   
      deallocate(sendVertexList % list)
      deallocate(sendVertexList)
      deallocate(recvVertexList % list)
      deallocate(recvVertexList)
   
#ifdef HAVE_ZOLTAN
#ifdef _MPI 
      allocate(xEdge(nlocal_edges))
      allocate(yEdge(nlocal_edges))
      allocate(zEdge(nlocal_edges))
      allocate(xVertex(nlocal_vertices))
      allocate(yVertex(nlocal_vertices))
      allocate(zVertex(nlocal_vertices))
#endif
#endif
    
      !
      ! Knowing which edges/vertices are owned by this block and which are actually read
      !   from the input or restart file, we can build exchange lists to perform 
      !   all-to-all field exchanges from process that reads a field to the processes that
      !   need them
      !
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                size(indexToEdgeIDField % array), nlocal_edges, &
                                indexToEdgeIDField % array, local_edge_list, &
                                sendEdgeList, recvEdgeList)
   
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                size(indexToVertexIDField % array), nlocal_vertices, &
                                indexToVertexIDField % array, local_vertex_list, &
                                sendVertexList, recvVertexList)

#ifdef HAVE_ZOLTAN
#ifdef _MPI 
      call mpas_dmpar_alltoall_field(domain % dminfo, xEdgeField % array, xEdge, &
                                size(xEdgeField % array), nlocal_edges, &
                                sendEdgeList, recvEdgeList)
      call mpas_dmpar_alltoall_field(domain % dminfo, yEdgeField % array, yEdge, &
                                size(yEdgeField % array), nlocal_edges, &
                                sendEdgeList, recvEdgeList)
      call mpas_dmpar_alltoall_field(domain % dminfo, zEdgeField % array, zEdge, &
                                size(zEdgeField % array), nlocal_edges, &
                                sendEdgeList, recvEdgeList)

      call mpas_dmpar_alltoall_field(domain % dminfo, xVertexField % array, xVertex, &
                                size(xVertexField % array), nlocal_vertices, &
                                sendVertexList, recvVertexList)
      call mpas_dmpar_alltoall_field(domain % dminfo, yVertexField % array, yVertex, &
                                size(yVertexField % array), nlocal_vertices, &
                                sendVertexList, recvVertexList)
      call mpas_dmpar_alltoall_field(domain % dminfo, zVertexField % array, zVertex, &
                                size(zVertexField % array), nlocal_vertices, &
                                sendVertexList, recvVertexList)
      !!!!!!!!!!!!!!!!!!
      !! Reorder edges
      !!!!!!!!!!!!!!!!!!
      call mpas_zoltan_order_loc_hsfc_edges(nOwnEdges,local_edge_list,3,xEdge,yEdge,zEdge)
      !!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!
      !! Reorder vertices
      !!!!!!!!!!!!!!!!!!
      call mpas_zoltan_order_loc_hsfc_verts(nOwnVertices,local_vertex_list,3,xVertex,yVertex,zVertex)
      !!!!!!!!!!!!!!!!!!

      deallocate(sendEdgeList % list)
      deallocate(sendEdgeList)
      deallocate(recvEdgeList % list)
      deallocate(recvEdgeList)
   
      deallocate(sendVertexList % list)
      deallocate(sendVertexList)
      deallocate(recvVertexList % list)
      deallocate(recvVertexList)
    
      !
      ! Knowing which edges/vertices are owned by this block and which are actually read
      !   from the input or restart file, we can build exchange lists to perform 
      !   all-to-all field exchanges from process that reads a field to the processes that
      !   need them
      !
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                size(indexToEdgeIDField % array), nlocal_edges, &
                                indexToEdgeIDField % array, local_edge_list, &
                                sendEdgeList, recvEdgeList)
   
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                size(indexToVertexIDField % array), nlocal_vertices, &
                                indexToVertexIDField % array, local_vertex_list, &
                                sendVertexList, recvVertexList)

#endif
#endif

      ! 
      ! Build ownership and exchange lists for vertical levels
      ! Essentially, process 0 owns all vertical levels when reading and writing,
      ! and it distributes them or gathers them to/from all other processes
      ! 
      if (domain % dminfo % my_proc_id == 0) then
         allocate(local_vertlevel_list(nVertLevels))
         do i=1,nVertLevels
            local_vertlevel_list(i) = i
         end do
      else
         allocate(local_vertlevel_list(0))
      end if
      allocate(needed_vertlevel_list(nVertLevels))
      do i=1,nVertLevels
         needed_vertlevel_list(i) = i
      end do

      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                size(local_vertlevel_list), size(needed_vertlevel_list), &
                                local_vertlevel_list, needed_vertlevel_list, &
                                sendVertLevelList, recvVertLevelList)

      deallocate(local_vertlevel_list)
      deallocate(needed_vertlevel_list)


      !
      ! Read and distribute all fields given ownership lists and exchange lists (maybe already in block?)
      !
      allocate(domain % blocklist)

      nCells = block_graph_2Halo % nVerticesTotal
      nEdges = nlocal_edges
      nVertices = nlocal_vertices

      call mpas_allocate_block(domain % blocklist, domain, domain%dminfo%my_proc_id, &
#include "dim_dummy_args.inc"
                         )

!!!!!!!!!!MGD HERE WE NEED TO READ IN indexTo*ID fields !!!!!!!!!!!!!!!!!
      call MPAS_io_inq_var(inputHandle, 'indexToCellID', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'indexToCellID', local_cell_list(1:nOwnCells), ierr=ierr)
      call mpas_io_get_var(inputHandle, 'indexToCellID', domain % blocklist % mesh % indexToCellID % array, ierr)

      call MPAS_io_inq_var(inputHandle, 'indexToEdgeID', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'indexToEdgeID', local_edge_list(1:nOwnEdges), ierr=ierr)
      call mpas_io_get_var(inputHandle, 'indexToEdgeID', domain % blocklist % mesh % indexToEdgeID % array, ierr)

      call MPAS_io_inq_var(inputHandle, 'indexToVertexID', ierr=ierr)
      call MPAS_io_set_var_indices(inputHandle, 'indexToVertexID', local_vertex_list(1:nOwnVertices), ierr=ierr)
      call mpas_io_get_var(inputHandle, 'indexToVertexID', domain % blocklist % mesh % indexToVertexID % array, ierr)

      domain % blocklist % mesh % nCellsSolve = nOwnCells
      domain % blocklist % mesh % nEdgesSolve = nOwnEdges
      domain % blocklist % mesh % nVerticesSolve = nOwnVertices
      domain % blocklist % mesh % nVertLevelsSolve = domain % blocklist % mesh % nVertLevels   ! No vertical decomp yet...

      call mpas_io_input_init(input_obj, domain % blocklist, domain % dminfo)


      !
      ! Read attributes
      !
      call MPAS_readStreamAtt(input_obj % io_stream, 'sphere_radius', r_sphere_radius, ierr)
      if (ierr /= MPAS_STREAM_NOERR) then
         write(0,*) 'Warning: Attribute sphere_radius not found in '//trim(input_obj % filename)
         write(0,*) '   Setting sphere_radius to 1.0'
         domain % blocklist % mesh % sphere_radius = 1.0
      else
         domain % blocklist % mesh % sphere_radius = r_sphere_radius
      end if

      call MPAS_readStreamAtt(input_obj % io_stream, 'on_a_sphere', c_on_a_sphere, ierr)
      if (ierr /= MPAS_STREAM_NOERR) then
         write(0,*) 'Warning: Attribute on_a_sphere not found in '//trim(input_obj % filename)
         write(0,*) '   Setting on_a_sphere to ''YES'''
         domain % blocklist % mesh % on_a_sphere = .true.
      else
         if (index(c_on_a_sphere, 'YES') /= 0) then
            domain % blocklist % mesh % on_a_sphere = .true.
         else
            domain % blocklist % mesh % on_a_sphere = .false.
         end if
      end if

      if (.not. config_do_restart) then
         input_obj % time = 1
      else
         !
         ! If doing a restart, we need to decide which time slice to read from the 
         !   restart file
         !
         input_obj % time = MPAS_seekStream(input_obj % io_stream, config_start_time, MPAS_STREAM_EXACT_TIME, timeStamp, ierr)
         if (ierr == MPAS_IO_ERR) then
            write(0,*) 'Error: restart file '//trim(filename)//' did not contain time '//trim(config_start_time)
            call mpas_dmpar_abort(domain % dminfo)
         end if
write(0,*) 'MGD DEBUGGING time = ', input_obj % time
         write(0,*) 'Restarting model from time ', timeStamp

      end if


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! Do the actual work of reading all fields in from the input or restart file
      ! For each field:
      !   1) Each process reads a contiguous range of cell/edge/vertex indices, which
      !      may not correspond with the cells/edges/vertices that are owned by the
      !      process
      !   2) All processes then send the global indices that were read to the 
      !      processes that own those indices based on 
      !      {send,recv}{Cell,Edge,Vertex,VertLevel}List
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      call mpas_read_and_distribute_fields(input_obj)

      call mpas_io_input_finalize(input_obj, domain % dminfo)

      call MPAS_io_close(inputHandle, ierr)

   
      !
      ! Work out halo exchange lists for cells, edges, and vertices
      ! NB: The next pointer in each element of, e.g., cellsToSend, acts as the head pointer of
      !     the list, since Fortran does not allow arrays of pointers
      !

      !--------- Create Cell Exchange Lists ---------!

      ! pass in neededList of ownedCells and halo layer 1 cells
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                nOwnCells, nCellsCumulative(2), &
                                block_graph_2Halo % vertexID(1:nOwnCells), block_graph_2Halo % vertexID(1 : nCellsCumulative(2)), &
                                domain % blocklist % parinfo % cellsToSend(1) % next, domain % blocklist % parinfo % cellsToRecv(1) % next)

      ! pass in neededList of ownedCells and halo layer 2 cells; offset of number of halo 1 cells is required
      offset = nCellsHalo(1)
      nTempIDs = nOwnCells + nCellsHalo(2)
      allocate(tempIDs(nTempIDs))
      tempIDs(1:nOwnCells) = block_graph_2Halo % vertexID(1:nOwnCells)
      tempIDs(nOwnCells+1:nTempIDs) = block_graph_2Halo % vertexID(nCellsCumulative(2)+1 : nCellsCumulative(3))
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                nOwnCells, nTempIDs, &
                                block_graph_2Halo % vertexID(1:nOwnCells), tempIDs, &
                                domain % blocklist % parinfo % cellsToSend(2) % next, domain % blocklist % parinfo % cellsToRecv(2) % next, &
                                offset)
      deallocate(tempIDs)


      !--------- Create Edge Exchange Lists ---------!

      ! pass in neededList of ownedEdges and ownedCell perimeter edges
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                nOwnEdges, nEdgesCumulative(2), &
                                local_edge_list(1:nOwnEdges), local_edge_list(1 : nEdgesCumulative(2)), &
                                domain % blocklist % parinfo % edgesToSend(1) % next, domain % blocklist % parinfo % edgesToRecv(1) % next)

      ! pass in neededList of owned edges and yet-to-be-included edges from halo 1 cells; offset of number of ownedCell perimeter edges is required
      offset = nEdgesHalo(1)
      nTempIDs = nOwnEdges + nEdgesHalo(2)
      allocate(tempIDs(nTempIDs))
      tempIDs(1:nOwnEdges) = local_edge_list(1:nOwnEdges)
      tempIDs(nOwnEdges+1:nTempIDs) = local_edge_list(nEdgesCumulative(2)+1 : nEdgesCumulative(3))
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                nOwnEdges, nTempIDs, &
                                local_edge_list(1:nOwnEdges), tempIDs, &  
                                domain % blocklist % parinfo % edgesToSend(2) % next, domain % blocklist % parinfo % edgesToRecv(2) % next, &
                                offset)
      deallocate(tempIDs)

      ! pass in neededList of owned edges and yet-to-be-included edges from halo 2 cells; offset of number of ownedCell perimeter edges and halo 1 edges is required
      offset = nEdgesHalo(1) + nEdgesHalo(2)
      nTempIDs = nOwnEdges + nEdgesHalo(3)
      allocate(tempIDs(nTempIDs))
      tempIDs(1:nOwnEdges) = local_edge_list(1:nOwnEdges)
      tempIDs(nOwnEdges+1:nTempIDs) = local_edge_list(nEdgesCumulative(3)+1 : nEdgesCumulative(4))
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                nOwnEdges, nTempIDs, &
                                local_edge_list(1:nOwnEdges), tempIDs, &  
                                domain % blocklist % parinfo % edgesToSend(3) % next, domain % blocklist % parinfo % edgesToRecv(3) % next, &
                                offset)
      deallocate(tempIDs)


      !--------- Create Vertex Exchange Lists ---------!


      ! pass in neededList of ownedVertices and ownedCell perimeter vertices
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                nOwnVertices, nVerticesCumulative(2), &
                                local_vertex_list(1:nOwnVertices), local_vertex_list(1 : nVerticesCumulative(2)), &
                                domain % blocklist % parinfo % verticesToSend(1) % next, domain % blocklist % parinfo % verticesToRecv(1) % next)

      ! pass in neededList of owned vertices and yet-to-be-included vertices from halo 1 cells; offset of number of ownedCell perimeter vertices is required
      offset = nVerticesHalo(1)
      nTempIDs = nOwnVertices + nVerticesHalo(2)
      allocate(tempIDs(nTempIDs))
      tempIDs(1:nOwnVertices) = local_vertex_list(1:nOwnVertices)
      tempIDs(nOwnVertices+1:nTempIDs) = local_vertex_list(nVerticesCumulative(2)+1 : nVerticesCumulative(3))
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                nOwnVertices, nTempIDs, &
                                local_vertex_list(1:nOwnVertices), tempIDs, &  
                                domain % blocklist % parinfo % verticesToSend(2) % next, domain % blocklist % parinfo % verticesToRecv(2) % next, &
                                offset)
      deallocate(tempIDs)

      ! pass in neededList of owned vertices and yet-to-be-included vertices from halo 2 cells; offset of number of ownedCell perimeter vertices and halo 1 vertices is required
      offset = nVerticesHalo(1) + nVerticesHalo(2)
      nTempIDs = nOwnVertices + nVerticesHalo(3)
      allocate(tempIDs(nTempIDs))
      tempIDs(1:nOwnVertices) = local_vertex_list(1:nOwnVertices)
      tempIDs(nOwnVertices+1:nTempIDs) = local_vertex_list(nVerticesCumulative(3)+1 : nVerticesCumulative(4))
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                nOwnVertices, nTempIDs, &
                                local_vertex_list(1:nOwnVertices), tempIDs, &  
                                domain % blocklist % parinfo % verticesToSend(3) % next, domain % blocklist % parinfo % verticesToRecv(3) % next, &
                                offset)
      deallocate(tempIDs)


      domain % blocklist % mesh % nCellsSolve = nOwnCells
      domain % blocklist % mesh % nEdgesSolve = nOwnEdges
      domain % blocklist % mesh % nVerticesSolve = ghostVertexStart-1
      domain % blocklist % mesh % nVertLevelsSolve = domain % blocklist % mesh % nVertLevels   ! No vertical decomp yet...

      ! Link the sendList and recvList pointers in each field type to the appropriate lists 
      !   in parinfo, e.g., cellsToSend and cellsToRecv; in future, it can also be extended to 
      !   link blocks of fields to eachother
      call mpas_create_field_links(domain % blocklist)


      !
      ! Exchange halos for all of the fields that were read from the input file
      !
      call mpas_exch_input_field_halos(domain, input_obj)

   
      !
      ! Rename vertices in cellsOnCell, edgesOnCell, etc. to local indices
      !
      allocate(cellIDSorted(2,domain % blocklist % mesh % nCells))
      allocate(edgeIDSorted(2,domain % blocklist % mesh % nEdges))
      allocate(vertexIDSorted(2,domain % blocklist % mesh % nVertices))

      do i=1,domain % blocklist % mesh % nCells
         cellIDSorted(1,i) = domain % blocklist % mesh % indexToCellID % array(i)
         cellIDSorted(2,i) = i
      end do
      call quicksort(block_graph_2Halo % nVerticesTotal, cellIDSorted)

      do i=1,domain % blocklist % mesh % nEdges
         edgeIDSorted(1,i) = domain % blocklist % mesh % indexToEdgeID % array(i)
         edgeIDSorted(2,i) = i
      end do
      call quicksort(nlocal_edges, edgeIDSorted)

      do i=1,domain % blocklist % mesh % nVertices
         vertexIDSorted(1,i) = domain % blocklist % mesh % indexToVertexID % array(i)
         vertexIDSorted(2,i) = i
      end do
      call quicksort(nlocal_vertices, vertexIDSorted)


      do i=1,domain % blocklist % mesh % nCells
         do j=1,domain % blocklist % mesh % nEdgesOnCell % array(i)

            k = mpas_binary_search(cellIDSorted, 2, 1, domain % blocklist % mesh % nCells, &
                              domain % blocklist % mesh % cellsOnCell % array(j,i))
            if (k <= domain % blocklist % mesh % nCells) then
               domain % blocklist % mesh % cellsOnCell % array(j,i) = cellIDSorted(2,k)
            else
               domain % blocklist % mesh % cellsOnCell % array(j,i) = domain % blocklist % mesh % nCells + 1
            end if

            k = mpas_binary_search(edgeIDSorted, 2, 1, domain % blocklist % mesh % nEdges, &
                              domain % blocklist % mesh % edgesOnCell % array(j,i))
            if (k <= domain % blocklist % mesh % nEdges) then
               domain % blocklist % mesh % edgesOnCell % array(j,i) = edgeIDSorted(2,k)
            else
               domain % blocklist % mesh % edgesOnCell % array(j,i) = domain % blocklist % mesh % nEdges + 1
            end if

            k = mpas_binary_search(vertexIDSorted, 2, 1, domain % blocklist % mesh % nVertices, &
                              domain % blocklist % mesh % verticesOnCell % array(j,i))
            if (k <= domain % blocklist % mesh % nVertices) then
               domain % blocklist % mesh % verticesOnCell % array(j,i) = vertexIDSorted(2,k)
            else
               domain % blocklist % mesh % verticesOnCell % array(j,i) = domain % blocklist % mesh % nVertices + 1
            end if

         end do
      end do

      do i=1,domain % blocklist % mesh % nEdges
         do j=1,2

            k = mpas_binary_search(cellIDSorted, 2, 1, domain % blocklist % mesh % nCells, &
                              domain % blocklist % mesh % cellsOnEdge % array(j,i))
            if (k <= domain % blocklist % mesh % nCells) then
               domain % blocklist % mesh % cellsOnEdge % array(j,i) = cellIDSorted(2,k)
            else
               domain % blocklist % mesh % cellsOnEdge % array(j,i) = domain % blocklist % mesh % nCells + 1
            end if

            k = mpas_binary_search(vertexIDSorted, 2, 1, domain % blocklist % mesh % nVertices, &
                              domain % blocklist % mesh % verticesOnEdge % array(j,i))
            if (k <= domain % blocklist % mesh % nVertices) then
               domain % blocklist % mesh % verticesOnEdge % array(j,i) = vertexIDSorted(2,k)
            else
               domain % blocklist % mesh % verticesOnEdge % array(j,i) = domain % blocklist % mesh % nVertices + 1
            end if

         end do

         do j=1,domain % blocklist % mesh % nEdgesOnEdge % array(i)

            k = mpas_binary_search(edgeIDSorted, 2, 1, domain % blocklist % mesh % nEdges, &
                              domain % blocklist % mesh % edgesOnEdge % array(j,i))
            if (k <= domain % blocklist % mesh % nEdges) then
               domain % blocklist % mesh % edgesOnEdge % array(j,i) = edgeIDSorted(2,k)
            else
               domain % blocklist % mesh % edgesOnEdge % array(j,i) = domain % blocklist % mesh % nEdges + 1
            end if

         end do
      end do

      do i=1,domain % blocklist % mesh % nVertices
         do j=1,vertexDegree

            k = mpas_binary_search(cellIDSorted, 2, 1, domain % blocklist % mesh % nCells, &
                              domain % blocklist % mesh % cellsOnVertex % array(j,i))
            if (k <= domain % blocklist % mesh % nCells) then
               domain % blocklist % mesh % cellsOnVertex % array(j,i) = cellIDSorted(2,k)
            else
               domain % blocklist % mesh % cellsOnVertex % array(j,i) = domain % blocklist % mesh % nCells + 1
            end if

            k = mpas_binary_search(edgeIDSorted, 2, 1, domain % blocklist % mesh % nEdges, &
                              domain % blocklist % mesh % edgesOnVertex % array(j,i))
            if (k <= domain % blocklist % mesh % nEdges) then
               domain % blocklist % mesh % edgesOnVertex % array(j,i) = edgeIDSorted(2,k)
            else
               domain % blocklist % mesh % edgesOnVertex % array(j,i) = domain % blocklist % mesh % nEdges + 1
            end if

         end do
      end do

      deallocate(cellIDSorted)
      deallocate(edgeIDSorted)
      deallocate(vertexIDSorted)

   
      !
      ! Deallocate fields, graphs, and other memory
      !
      deallocate(indexToCellIDField % ioinfo)
      deallocate(indexToCellIDField % array)
#ifdef HAVE_ZOLTAN
#ifdef _MPI 
      deallocate(xCellField % ioinfo)
      deallocate(xCellField % array)
      deallocate(yCellField % ioinfo)
      deallocate(yCellField % array)
      deallocate(zCellField % ioinfo)
      deallocate(zCellField % array)
#endif
#endif
      deallocate(indexToEdgeIDField % ioinfo)
      deallocate(indexToEdgeIDField % array)
      deallocate(indexToVertexIDField % ioinfo)
      deallocate(indexToVertexIDField % array)
      deallocate(cellsOnCellField % ioinfo)
      deallocate(cellsOnCellField % array)
      deallocate(edgesOnCellField % ioinfo)
      deallocate(edgesOnCellField % array)
      deallocate(verticesOnCellField % ioinfo)
      deallocate(verticesOnCellField % array)
      deallocate(cellsOnEdgeField % ioinfo)
      deallocate(cellsOnEdgeField % array)
      deallocate(cellsOnVertexField % ioinfo)
      deallocate(cellsOnVertexField % array)
      deallocate(cellsOnCell_0Halo)
      deallocate(nEdgesOnCell_0Halo)
      deallocate(indexToCellID_0Halo)
      deallocate(cellsOnEdge_2Halo)
      deallocate(cellsOnVertex_2Halo)
      deallocate(nEdgesOnCell_2Halo)
      deallocate(edgesOnCell_2Halo)
      deallocate(verticesOnCell_2Halo)
      deallocate(block_graph_0Halo % vertexID)
      deallocate(block_graph_0Halo % nAdjacent)
      deallocate(block_graph_0Halo % adjacencyList)
#ifdef HAVE_ZOLTAN
#ifdef _MPI
      deallocate(xCell)
      deallocate(yCell)
      deallocate(zCell)
#endif
#endif
   end subroutine mpas_input_state_for_domain


   !CR:TODO: an identical subroutine is found in module_io_output - merge
   subroutine mpas_insert_string_suffix(stream, suffix, filename)

      implicit none

      character (len=*), intent(in) :: stream
      character (len=*), intent(in) :: suffix
      character (len=*), intent(out) :: filename
      integer :: length, i

      filename = trim(stream) // '.' // trim(suffix)

      length = len_trim(stream)
      do i=length-1,1,-1
         if(stream(i:i) == '.') then
            filename = trim(stream(:i)) // trim(suffix) // trim(stream(i:))
            exit
         end if
      end do

      do i=1,len_trim(filename)
         if (filename(i:i) == ':') filename(i:i) = '.'
      end do

   end subroutine mpas_insert_string_suffix


   subroutine mpas_read_and_distribute_fields(input_obj)
      
      implicit none

      type (io_input_object), intent(inout) :: input_obj

      integer :: ierr


      call MPAS_readStream(input_obj % io_stream, 1, ierr)


   end subroutine mpas_read_and_distribute_fields



   subroutine mpas_io_input_init(input_obj, blocklist, dminfo)
 
      implicit none

      type (io_input_object), intent(inout) :: input_obj
      type (block_type), intent(in) :: blocklist
      type (dm_info), intent(in) :: dminfo
 
      integer :: nferr
 
      call MPAS_createStream(input_obj % io_stream, trim(input_obj % filename), MPAS_IO_PNETCDF, MPAS_IO_READ, 1, nferr)
      if (nferr /= MPAS_STREAM_NOERR) then
         write(0,*) ' '
         if (input_obj % stream == STREAM_RESTART) then
            write(0,*) 'Error opening restart file ''', trim(input_obj % filename), ''''
         else if (input_obj % stream == STREAM_INPUT) then
            write(0,*) 'Error opening input file ''', trim(input_obj % filename), ''''
         else if (input_obj % stream == STREAM_SFC) then
            write(0,*) 'Error opening sfc file ''', trim(input_obj % filename), ''''
         end if
         write(0,*) ' '
         call mpas_dmpar_abort(dminfo)
      end if

#include "add_input_fields.inc"

   end subroutine mpas_io_input_init

  
   subroutine mpas_io_input_get_dimension(input_obj, dimname, dimsize)

      implicit none

      type (io_input_object), intent(in) :: input_obj
      character (len=*), intent(in) :: dimname
      integer, intent(out) :: dimsize

!include "get_dimension_by_name.inc"

   end subroutine mpas_io_input_get_dimension

   
   subroutine mpas_io_input_get_att_real(input_obj, attname, attvalue)
      
      implicit none

      type (io_input_object), intent(in) :: input_obj
      character (len=*), intent(in) :: attname
      real (kind=RKIND), intent(out) :: attvalue

      integer :: nferr

   end subroutine mpas_io_input_get_att_real

   
   subroutine mpas_io_input_get_att_text(input_obj, attname, attvalue)
      
      implicit none

      type (io_input_object), intent(in) :: input_obj
      character (len=*), intent(in) :: attname
      character (len=*), intent(out) :: attvalue

      integer :: nferr

   end subroutine mpas_io_input_get_att_text


   subroutine mpas_exch_input_field_halos(domain, input_obj)

      implicit none

      type (domain_type), intent(inout) :: domain
      type (io_input_object), intent(inout) :: input_obj

#include "exchange_input_field_halos.inc"

   end subroutine mpas_exch_input_field_halos


   subroutine mpas_io_input_finalize(input_obj, dminfo)
 
      implicit none
 
      type (io_input_object), intent(inout) :: input_obj
      type (dm_info), intent(in) :: dminfo

      integer :: nferr
 
      call MPAS_closeStream(input_obj % io_stream, nferr)
 
   end subroutine mpas_io_input_finalize
 
end module mpas_io_input
