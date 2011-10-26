module mpas_io_input

   use mpas_grid_types
   use mpas_dmpar
   use mpas_block_decomp
   use mpas_sort
   use mpas_configure
   use mpas_timekeeping


#ifdef HAVE_ZOLTAN
   use mpas_zoltan_interface
#endif

   integer, parameter :: STREAM_INPUT=1, STREAM_SFC=2, STREAM_RESTART=3

   type io_input_object
      character (len=1024) :: filename
      integer :: rd_ncid
      integer :: stream

      integer :: time

#include "io_input_obj_decls.inc"
   end type io_input_object


   interface mpas_io_input_field
      module procedure mpas_io_input_field0d_real
      module procedure mpas_io_input_field1d_real
      module procedure mpas_io_input_field2d_real
      module procedure mpas_io_input_field3d_real
      module procedure mpas_io_input_field1d_integer
      module procedure mpas_io_input_field2d_integer
      module procedure mpas_io_input_field0d_char
      module procedure mpas_io_input_field1d_char
   end interface mpas_io_input_field

   interface mpas_io_input_field_time
      module procedure mpas_io_input_field0d_real_time
      module procedure mpas_io_input_field1d_real_time
      module procedure mpas_io_input_field2d_real_time
      module procedure mpas_io_input_field3d_real_time
      module procedure mpas_io_input_field1d_integer_time
      module procedure mpas_io_input_field0d_char_time
      module procedure mpas_io_input_field1d_char_time
   end interface mpas_io_input_field_time

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
   
      integer, dimension(:), pointer :: indexToCellID_0Halo
      integer, dimension(:), pointer :: nEdgesOnCell_0Halo
      integer, dimension(:,:), pointer :: cellsOnCell_0Halo
   
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
      call mpas_io_input_init(input_obj, domain % dminfo)
   

      !
      ! Read global number of cells/edges/vertices
      !
#include "read_dims.inc"
   
      !
      ! Determine the range of cells/edges/vertices that a processor will initially read
      !   from the input file
      !
      call mpas_dmpar_get_index_range(domain % dminfo, 1, nCells, readCellStart, readCellEnd)   
      nReadCells    = readCellEnd - readCellStart + 1
   
      call mpas_dmpar_get_index_range(domain % dminfo, 1, nEdges, readEdgeStart, readEdgeEnd)   
      nReadEdges    = readEdgeEnd - readEdgeStart + 1
   
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
      call mpas_io_input_field(input_obj, indexToCellIDField)
   
#ifdef HAVE_ZOLTAN
#ifdef _MPI 
      ! Cell x-coordinates (in 3d Cartesian space)
      allocate(xCellField % ioinfo)
      xCellField % ioinfo % fieldName = 'xCell'
      xCellField % ioinfo % start(1) = readCellStart
      xCellField % ioinfo % count(1) = nReadCells
      allocate(xCellField % array(nReadCells))
      call mpas_io_input_field(input_obj, xCellField)

      ! Cell y-coordinates (in 3d Cartesian space)
      allocate(yCellField % ioinfo)
      yCellField % ioinfo % fieldName = 'yCell'
      yCellField % ioinfo % start(1) = readCellStart
      yCellField % ioinfo % count(1) = nReadCells
      allocate(yCellField % array(nReadCells))
      call mpas_io_input_field(input_obj, yCellField)

      ! Cell z-coordinates (in 3d Cartesian space)
      allocate(zCellField % ioinfo)
      zCellField % ioinfo % fieldName = 'zCell'
      zCellField % ioinfo % start(1) = readCellStart
      zCellField % ioinfo % count(1) = nReadCells
      allocate(zCellField % array(nReadCells))
      call mpas_io_input_field(input_obj, zCellField)
#endif
#endif


      ! Global edge indices
      allocate(indexToEdgeIDField % ioinfo)
      indexToEdgeIDField % ioinfo % fieldName = 'indexToEdgeID'
      indexToEdgeIDField % ioinfo % start(1) = readEdgeStart
      indexToEdgeIDField % ioinfo % count(1) = nReadEdges
      allocate(indexToEdgeIDField % array(nReadEdges))
      call mpas_io_input_field(input_obj, indexToEdgeIDField)
   
#ifdef HAVE_ZOLTAN
#ifdef _MPI 
      ! Edge x-coordinates (in 3d Cartesian space)
      allocate(xEdgeField % ioinfo)
      xEdgeField % ioinfo % fieldName = 'xEdge'
      xEdgeField % ioinfo % start(1) = readEdgeStart
      xEdgeField % ioinfo % count(1) = nReadEdges
      allocate(xEdgeField % array(nReadEdges))
      call mpas_io_input_field(input_obj, xEdgeField)

      ! Edge y-coordinates (in 3d Cartesian space)
      allocate(yEdgeField % ioinfo)
      yEdgeField % ioinfo % fieldName = 'yEdge'
      yEdgeField % ioinfo % start(1) = readEdgeStart
      yEdgeField % ioinfo % count(1) = nReadEdges
      allocate(yEdgeField % array(nReadEdges))
      call mpas_io_input_field(input_obj, yEdgeField)

      ! Edge z-coordinates (in 3d Cartesian space)
      allocate(zEdgeField % ioinfo)
      zEdgeField % ioinfo % fieldName = 'zEdge'
      zEdgeField % ioinfo % start(1) = readEdgeStart
      zEdgeField % ioinfo % count(1) = nReadEdges
      allocate(zEdgeField % array(nReadEdges))
      call mpas_io_input_field(input_obj, zEdgeField)
#endif
#endif

      ! Global vertex indices
      allocate(indexToVertexIDField % ioinfo)
      indexToVertexIDField % ioinfo % fieldName = 'indexToVertexID'
      indexToVertexIDField % ioinfo % start(1) = readVertexStart
      indexToVertexIDField % ioinfo % count(1) = nReadVertices
      allocate(indexToVertexIDField % array(nReadVertices))
      call mpas_io_input_field(input_obj, indexToVertexIDField)
   
#ifdef HAVE_ZOLTAN
#ifdef _MPI
      ! Vertex x-coordinates (in 3d Cartesian space)
      allocate(xVertexField % ioinfo)
      xVertexField % ioinfo % fieldName = 'xVertex'
      xVertexField % ioinfo % start(1) = readVertexStart
      xVertexField % ioinfo % count(1) = nReadVertices
      allocate(xVertexField % array(nReadVertices))
      call mpas_io_input_field(input_obj, xVertexField)

      ! Vertex y-coordinates (in 3d Cartesian space)
      allocate(yVertexField % ioinfo)
      yVertexField % ioinfo % fieldName = 'yVertex'
      yVertexField % ioinfo % start(1) = readVertexStart
      yVertexField % ioinfo % count(1) = nReadVertices
      allocate(yVertexField % array(nReadVertices))
      call mpas_io_input_field(input_obj, yVertexField)

      ! Vertex z-coordinates (in 3d Cartesian space)
      allocate(zVertexField % ioinfo)
      zVertexField % ioinfo % fieldName = 'zVertex'
      zVertexField % ioinfo % start(1) = readVertexStart
      zVertexField % ioinfo % count(1) = nReadVertices
      allocate(zVertexField % array(nReadVertices))
      call mpas_io_input_field(input_obj, zVertexField)
#endif
#endif

      ! Number of cell/edges/vertices adjacent to each cell
      allocate(nEdgesOnCellField % ioinfo)
      nEdgesOnCellField % ioinfo % fieldName = 'nEdgesOnCell'
      nEdgesOnCellField % ioinfo % start(1) = readCellStart
      nEdgesOnCellField % ioinfo % count(1) = nReadCells
      allocate(nEdgesOnCellField % array(nReadCells))
      call mpas_io_input_field(input_obj, nEdgesOnCellField)
   
      ! Global indices of cells adjacent to each cell
      allocate(cellsOnCellField % ioinfo)
      cellsOnCellField % ioinfo % fieldName = 'cellsOnCell'
      cellsOnCellField % ioinfo % start(1) = 1
      cellsOnCellField % ioinfo % start(2) = readCellStart
      cellsOnCellField % ioinfo % count(1) = maxEdges
      cellsOnCellField % ioinfo % count(2) = nReadCells
      allocate(cellsOnCellField % array(maxEdges,nReadCells))
      call mpas_io_input_field(input_obj, cellsOnCellField)
   
      ! Global indices of edges adjacent to each cell
      allocate(edgesOnCellField % ioinfo)
      edgesOnCellField % ioinfo % fieldName = 'edgesOnCell'
      edgesOnCellField % ioinfo % start(1) = 1
      edgesOnCellField % ioinfo % start(2) = readCellStart
      edgesOnCellField % ioinfo % count(1) = maxEdges
      edgesOnCellField % ioinfo % count(2) = nReadCells
      allocate(edgesOnCellField % array(maxEdges,nReadCells))
      call mpas_io_input_field(input_obj, edgesOnCellField)
   
      ! Global indices of vertices adjacent to each cell
      allocate(verticesOnCellField % ioinfo)
      verticesOnCellField % ioinfo % fieldName = 'verticesOnCell'
      verticesOnCellField % ioinfo % start(1) = 1
      verticesOnCellField % ioinfo % start(2) = readCellStart
      verticesOnCellField % ioinfo % count(1) = maxEdges
      verticesOnCellField % ioinfo % count(2) = nReadCells
      allocate(verticesOnCellField % array(maxEdges,nReadCells))
      call mpas_io_input_field(input_obj, verticesOnCellField)
   
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
      call mpas_io_input_field(input_obj, cellsOnEdgeField)
   
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
      call mpas_io_input_field(input_obj, cellsOnVertexField)
   
   
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
      call mpas_block_decomp_cells_for_proc(domain % dminfo, partial_global_graph_info, local_cell_list)

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
      allocate(edgesOnCell_2Halo(maxEdges, block_graph_2Halo % nVerticesTotal))
      allocate(verticesOnCell_2Halo(maxEdges, block_graph_2Halo % nVerticesTotal))
   
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                size(indexToCellIDField % array), block_graph_2Halo % nVerticesTotal, &
                                indexToCellIDField % array, block_graph_2Halo % vertexID, &
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
   
   
      call mpas_block_decomp_partitioned_edge_list(block_graph_2Halo % nVertices, &
                                              block_graph_2Halo % vertexID(1:block_graph_2Halo % nVertices), &
                                              2, nlocal_edges, cellsOnEdge_2Halo, local_edge_list, ghostEdgeStart)
      call mpas_block_decomp_partitioned_edge_list(block_graph_2Halo % nVertices, &
                                              block_graph_2Halo % vertexID(1:block_graph_2Halo % nVertices), &
                                              vertexDegree, nlocal_vertices, cellsOnVertex_2Halo, local_vertex_list, ghostVertexStart)


      ! At this point, local_edge_list(1;ghostEdgeStart-1) contains all of the owned edges for this block
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
      call mpas_zoltan_order_loc_hsfc_edges(ghostEdgeStart-1,local_edge_list,3,xEdge,yEdge,zEdge)
      !!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!
      !! Reorder vertices
      !!!!!!!!!!!!!!!!!!
      call mpas_zoltan_order_loc_hsfc_verts(ghostVertexStart-1,local_vertex_list,3,xVertex,yVertex,zVertex)
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

      call mpas_allocate_block(domain % blocklist, domain, &
#include "dim_dummy_args.inc"
                         )

      !
      ! Read attributes
      !
      call mpas_io_input_get_att_text(input_obj, 'on_a_sphere', c_on_a_sphere)
      call mpas_io_input_get_att_real(input_obj, 'sphere_radius', r_sphere_radius)
      if (index(c_on_a_sphere, 'YES') /= 0) then
         domain % blocklist % mesh % on_a_sphere = .true.
      else
         domain % blocklist % mesh % on_a_sphere = .false.
      end if
      domain % blocklist % mesh % sphere_radius = r_sphere_radius

      if (.not. config_do_restart) then
         input_obj % time = 1
      else
         input_obj % time = 1

         !
         ! If doing a restart, we need to decide which time slice to read from the 
         !   restart file
         !
         if (input_obj % rdLocalTime <= 0) then
            write(0,*) 'Error: Couldn''t find any times in restart file.'
            call mpas_dmpar_abort(domain % dminfo)
         end if
         if (domain % dminfo % my_proc_id == IO_NODE) then
            allocate(xtime % ioinfo)
            xtime % ioinfo % start(1) = 1
            xtime % ioinfo % count(1) = input_obj % rdLocalTime
            allocate(xtime % array(input_obj % rdLocalTime))

            xtime % ioinfo % fieldName = 'xtime'
            call mpas_io_input_field(input_obj, xtime)

            call mpas_set_timeInterval(interval=minTimeDiff, DD=10000)
            call mpas_set_time(curr_time=startTime, dateTimeString=config_start_time)

            do i=1,input_obj % rdLocalTime
               call mpas_set_time(curr_time=sliceTime, dateTimeString=xtime % array(i))
               timeDiff = abs(sliceTime - startTime)
               if (timeDiff < minTimeDiff) then
                  minTimeDiff = timeDiff
                  input_obj % time = i
               end if
            end do

            ! require restart time to exactly match start time (this error should never be reached as we have by this point opened the restart file with a name containing the startTime)
            if(sliceTime /= startTime) then
               write(0,*) "Error: restart file ", filename, " did not contain time ", config_start_time
               call mpas_dmpar_abort(domain % dminfo)
            end if

            timeStamp = xtime % array(input_obj % time)

            deallocate(xtime % ioinfo)
            deallocate(xtime % array)
         end if

         call mpas_dmpar_bcast_int(domain % dminfo, input_obj % time)
         call mpas_dmpar_bcast_char(domain % dminfo, timeStamp)

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
      call mpas_read_and_distribute_fields(domain % dminfo, input_obj, domain % blocklist, &
                                      readCellStart, nReadCells, readEdgeStart, nReadEdges, readVertexStart, nReadVertices, &
                                      readVertLevelStart, nReadVertLevels, &
                                      sendCellList, recvCellList, sendEdgeList, recvEdgeList, sendVertexList, recvVertexList, &
                                      sendVertLevelList, recvVertLevelList) 


      call mpas_io_input_finalize(input_obj, domain % dminfo)

   
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
!               domain % blocklist % mesh % cellsOnCell % array(j,i) = 0
            end if

            k = mpas_binary_search(edgeIDSorted, 2, 1, domain % blocklist % mesh % nEdges, &
                              domain % blocklist % mesh % edgesOnCell % array(j,i))
            if (k <= domain % blocklist % mesh % nEdges) then
               domain % blocklist % mesh % edgesOnCell % array(j,i) = edgeIDSorted(2,k)
            else
               domain % blocklist % mesh % edgesOnCell % array(j,i) = domain % blocklist % mesh % nEdges + 1
!               domain % blocklist % mesh % edgesOnCell % array(j,i) = 0
            end if

            k = mpas_binary_search(vertexIDSorted, 2, 1, domain % blocklist % mesh % nVertices, &
                              domain % blocklist % mesh % verticesOnCell % array(j,i))
            if (k <= domain % blocklist % mesh % nVertices) then
               domain % blocklist % mesh % verticesOnCell % array(j,i) = vertexIDSorted(2,k)
            else
               domain % blocklist % mesh % verticesOnCell % array(j,i) = domain % blocklist % mesh % nVertices + 1
!               domain % blocklist % mesh % verticesOnCell % array(j,i) = 0
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
!               domain % blocklist % mesh % cellsOnEdge % array(j,i) = 0
            end if

            k = mpas_binary_search(vertexIDSorted, 2, 1, domain % blocklist % mesh % nVertices, &
                              domain % blocklist % mesh % verticesOnEdge % array(j,i))
            if (k <= domain % blocklist % mesh % nVertices) then
               domain % blocklist % mesh % verticesOnEdge % array(j,i) = vertexIDSorted(2,k)
            else
               domain % blocklist % mesh % verticesOnEdge % array(j,i) = domain % blocklist % mesh % nVertices + 1
!               domain % blocklist % mesh % verticesOnEdge % array(j,i) = 0
            end if

         end do

         do j=1,domain % blocklist % mesh % nEdgesOnEdge % array(i)

            k = mpas_binary_search(edgeIDSorted, 2, 1, domain % blocklist % mesh % nEdges, &
                              domain % blocklist % mesh % edgesOnEdge % array(j,i))
            if (k <= domain % blocklist % mesh % nEdges) then
               domain % blocklist % mesh % edgesOnEdge % array(j,i) = edgeIDSorted(2,k)
            else
               domain % blocklist % mesh % edgesOnEdge % array(j,i) = domain % blocklist % mesh % nEdges + 1
!               domain % blocklist % mesh % edgesOnEdge % array(j,i) = 0
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
!               domain % blocklist % mesh % cellsOnVertex % array(j,i) = 0
            end if

            k = mpas_binary_search(edgeIDSorted, 2, 1, domain % blocklist % mesh % nEdges, &
                              domain % blocklist % mesh % edgesOnVertex % array(j,i))
            if (k <= domain % blocklist % mesh % nEdges) then
               domain % blocklist % mesh % edgesOnVertex % array(j,i) = edgeIDSorted(2,k)
            else
               domain % blocklist % mesh % edgesOnVertex % array(j,i) = domain % blocklist % mesh % nEdges + 1
!               domain % blocklist % mesh % edgesOnVertex % array(j,i) = 0
            end if

         end do
      end do

      deallocate(cellIDSorted)
      deallocate(edgeIDSorted)
      deallocate(vertexIDSorted)


      !
      ! Work out halo exchange lists for cells, edges, and vertices
      !
      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                block_graph_2Halo % nVertices, block_graph_2Halo % nVerticesTotal, &
                                block_graph_2Halo % vertexID(1:block_graph_2Halo % nVertices), block_graph_2Halo % vertexID, &
                                domain % blocklist % parinfo % cellsToSend, domain % blocklist % parinfo % cellsToRecv)

      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                ghostEdgeStart-1, nlocal_edges, &
                                local_edge_list(1:ghostEdgeStart-1), local_edge_list, &
                                domain % blocklist % parinfo % edgesToSend, domain % blocklist % parinfo % edgesToRecv)

      call mpas_dmpar_get_owner_list(domain % dminfo, &
                                ghostVertexStart-1, nlocal_vertices, &
                                local_vertex_list(1:ghostVertexStart-1), local_vertex_list, &
                                domain % blocklist % parinfo % verticesToSend, domain % blocklist % parinfo % verticesToRecv)

      domain % blocklist % mesh % nCellsSolve = block_graph_2Halo % nVertices
      domain % blocklist % mesh % nEdgesSolve = ghostEdgeStart-1
      domain % blocklist % mesh % nVerticesSolve = ghostVertexStart-1
      domain % blocklist % mesh % nVertLevelsSolve = domain % blocklist % mesh % nVertLevels   ! No vertical decomp yet...

   
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

   end subroutine mpas_insert_string_suffix


   subroutine mpas_read_and_distribute_fields(dminfo, input_obj, block, &
                                     readCellsStart, readCellsCount, &
                                     readEdgesStart, readEdgesCount, &
                                     readVerticesStart, readVerticesCount, &
                                     readVertLevelsStart, readVertLevelsCount, &
                                     sendCellsList, recvCellsList, &
                                     sendEdgesList, recvEdgesList, &
                                     sendVerticesList, recvVerticesList, &
                                     sendVertLevelsList, recvVertLevelsList)
      
      implicit none

      type (dm_info), intent(in) :: dminfo
      type (io_input_object), intent(in) :: input_obj
      type (block_type), intent(inout) :: block
      integer, intent(in) :: readCellsStart, readCellsCount, readEdgesStart, readEdgesCount, readVerticesStart, readVerticesCount
      integer, intent(in) :: readVertLevelsStart, readVertLevelsCount
      type (exchange_list), pointer :: sendCellsList, recvCellsList
      type (exchange_list), pointer :: sendEdgesList, recvEdgesList
      type (exchange_list), pointer :: sendVerticesList, recvVerticesList
      type (exchange_list), pointer :: sendVertLevelsList, recvVertLevelsList

      type (field1dInteger) :: int1d
      type (field2dInteger) :: int2d
      type (field0dReal) :: real0d
      type (field1dReal) :: real1d
      type (field2dReal) :: real2d
      type (field3dReal) :: real3d
      type (field0dChar) :: char0d
      type (field1dChar) :: char1d

      integer :: i1, i2, i3, i4

      integer, dimension(:), pointer :: super_int1d
      integer, dimension(:,:), pointer :: super_int2d
      real (kind=RKIND) :: super_real0d
      real (kind=RKIND), dimension(:), pointer :: super_real1d
      real (kind=RKIND), dimension(:,:), pointer :: super_real2d
      real (kind=RKIND), dimension(:,:,:), pointer :: super_real3d
      character (len=64) :: super_char0d
      character (len=64), dimension(:), pointer :: super_char1d

      integer :: i, k

#include "nondecomp_dims.inc"

      allocate(int1d % ioinfo)
      allocate(int2d % ioinfo)
      allocate(real0d % ioinfo)
      allocate(real1d % ioinfo)
      allocate(real2d % ioinfo)
      allocate(real3d % ioinfo)
      allocate(char0d % ioinfo)
      allocate(char1d % ioinfo)


#include "io_input_fields.inc"

#include "nondecomp_dims_dealloc.inc"

   end subroutine mpas_read_and_distribute_fields



   subroutine mpas_io_input_init(input_obj, dminfo)
 
      implicit none

      type (io_input_object), intent(inout) :: input_obj
      type (dm_info), intent(in) :: dminfo
 
      include 'netcdf.inc'
 
      integer :: nferr
 
 
#ifdef OFFSET64BIT
      nferr = nf_open(trim(input_obj % filename), ior(NF_SHARE,NF_64BIT_OFFSET), input_obj % rd_ncid)
#else
      nferr = nf_open(trim(input_obj % filename), NF_SHARE, input_obj % rd_ncid)
#endif

      if (nferr /= NF_NOERR) then
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
 
#include "netcdf_read_ids.inc"

   end subroutine mpas_io_input_init

  
   subroutine mpas_io_input_get_dimension(input_obj, dimname, dimsize)

      implicit none

      type (io_input_object), intent(in) :: input_obj
      character (len=*), intent(in) :: dimname
      integer, intent(out) :: dimsize

#include "get_dimension_by_name.inc"

   end subroutine mpas_io_input_get_dimension

   
   subroutine mpas_io_input_get_att_real(input_obj, attname, attvalue)
      
      implicit none

      type (io_input_object), intent(in) :: input_obj
      character (len=*), intent(in) :: attname
      real (kind=RKIND), intent(out) :: attvalue

      include 'netcdf.inc'

      integer :: nferr

      if (RKIND == 8) then
         nferr = nf_get_att_double(input_obj % rd_ncid, NF_GLOBAL, attname, attvalue)
      else
         nferr = nf_get_att_real(input_obj % rd_ncid, NF_GLOBAL, attname, attvalue)
      end if
      if (nferr /= NF_NOERR) then
         write(0,*) 'Warning: Attribute '//trim(attname)//&
           ' not found in '//trim(input_obj % filename)
         if (index(attname, 'sphere_radius') /= 0) then
            write(0,*) '   Setting '//trim(attname)//' to 1.0'
            attvalue = 1.0
         end if
      end if

   end subroutine mpas_io_input_get_att_real

   
   subroutine mpas_io_input_get_att_text(input_obj, attname, attvalue)
      
      implicit none

      type (io_input_object), intent(in) :: input_obj
      character (len=*), intent(in) :: attname
      character (len=*), intent(out) :: attvalue

      include 'netcdf.inc'

      integer :: nferr

      nferr = nf_get_att_text(input_obj % rd_ncid, NF_GLOBAL, attname, attvalue)
      if (nferr /= NF_NOERR) then
         write(0,*) 'Warning: Attribute '//trim(attname)//&
            ' not found in '//trim(input_obj % filename)
         if (index(attname, 'on_a_sphere') /= 0) then
            write(0,*) '   Setting '//trim(attname)//' to ''YES'''
            attvalue = 'YES'
         end if
      end if

   end subroutine mpas_io_input_get_att_text


   subroutine mpas_io_input_field0d_real(input_obj, field)
 
      implicit none

      type (io_input_object), intent(in) :: input_obj      
      type (field0dReal), intent(inout) :: field
 
      include 'netcdf.inc'
 
      integer :: nferr
      integer :: varID
      integer, dimension(1) :: start1, count1
 
      start1(1) = 1
      count1(1) = 1

#include "input_field0dreal.inc"

#if (RKIND == 8)
      nferr = nf_get_vara_double(input_obj % rd_ncid, varID, start1, count1, field % scalar)
#else
      nferr = nf_get_vara_real(input_obj % rd_ncid, varID, start1, count1, field % scalar)
#endif
 
   end subroutine mpas_io_input_field0d_real


   subroutine mpas_io_input_field1d_real(input_obj, field)
 
      implicit none

      type (io_input_object), intent(in) :: input_obj      
      type (field1dReal), intent(inout) :: field
 
      include 'netcdf.inc'
 
      integer :: nferr
      integer :: varID
      integer, dimension(1) :: start1, count1
 
      start1(1) = field % ioinfo % start(1)
      count1(1) = field % ioinfo % count(1)

      !
      ! Special case: we may want to read the xtime variable across the
      !   time dimension as a 1d array.
      !
      if (trim(field % ioinfo % fieldName) == 'xtime') then
         varID = input_obj % rdVarIDxtime
      end if
 
#include "input_field1dreal.inc"

#if (RKIND == 8)
      nferr = nf_get_vara_double(input_obj % rd_ncid, varID, start1, count1, field % array)
#else
      nferr = nf_get_vara_real(input_obj % rd_ncid, varID, start1, count1, field % array)
#endif
 
   end subroutine mpas_io_input_field1d_real


   subroutine mpas_io_input_field2d_real(input_obj, field)
 
      implicit none

      type (io_input_object), intent(in) :: input_obj      
      type (field2dReal), intent(inout) :: field
 
      include 'netcdf.inc'
 
      integer :: nferr
      integer :: varID
      integer, dimension(2) :: start2, count2
 
      start2(1) = field % ioinfo % start(1)
      start2(2) = field % ioinfo % start(2)
      count2(1) = field % ioinfo % count(1)
      count2(2) = field % ioinfo % count(2)
 
#include "input_field2dreal.inc"

#if (RKIND == 8)
      nferr = nf_get_vara_double(input_obj % rd_ncid, varID, start2, count2, field % array)
#else
      nferr = nf_get_vara_real(input_obj % rd_ncid, varID, start2, count2, field % array)
#endif

   end subroutine mpas_io_input_field2d_real


   subroutine mpas_io_input_field3d_real(input_obj, field)
 
      implicit none

      type (io_input_object), intent(in) :: input_obj      
      type (field3dReal), intent(inout) :: field
 
      include 'netcdf.inc'
 
      integer :: nferr
      integer :: varID
      integer, dimension(3) :: start3, count3
 
      start3(1) = field % ioinfo % start(1)
      start3(2) = field % ioinfo % start(2)
      start3(3) = field % ioinfo % start(3)
      count3(1) = field % ioinfo % count(1)
      count3(2) = field % ioinfo % count(2)
      count3(3) = field % ioinfo % count(3)
 
#include "input_field3dreal.inc"

#if (RKIND == 8)
      nferr = nf_get_vara_double(input_obj % rd_ncid, varID, start3, count3, field % array)
#else
      nferr = nf_get_vara_real(input_obj % rd_ncid, varID, start3, count3, field % array)
#endif

   end subroutine mpas_io_input_field3d_real


   subroutine mpas_io_input_field0d_real_time(input_obj, field)
 
      implicit none

      type (io_input_object), intent(in) :: input_obj      
      type (field0dReal), intent(inout) :: field
 
      include 'netcdf.inc'
 
      integer :: nferr
      integer :: varID
      integer, dimension(1) :: start1, count1
 
      start1(1) = input_obj % time
      count1(1) = 1
 
#include "input_field0dreal_time.inc"

#if (RKIND == 8)
      nferr = nf_get_vara_double(input_obj % rd_ncid, varID, start1, count1, field % scalar)
#else
      nferr = nf_get_vara_real(input_obj % rd_ncid, varID, start1, count1, field % scalar)
#endif

   end subroutine mpas_io_input_field0d_real_time


   subroutine mpas_io_input_field1d_real_time(input_obj, field)
 
      implicit none

      type (io_input_object), intent(in) :: input_obj      
      type (field1dReal), intent(inout) :: field
 
      include 'netcdf.inc'
 
      integer :: nferr
      integer :: varID
      integer, dimension(2) :: start2, count2
 
      start2(1) = field % ioinfo % start(1)
      start2(2) = input_obj % time
      count2(1) = field % ioinfo % count(1)
      count2(2) = 1
 
#include "input_field1dreal_time.inc"

#if (RKIND == 8)
      nferr = nf_get_vara_double(input_obj % rd_ncid, varID, start2, count2, field % array)
#else
      nferr = nf_get_vara_real(input_obj % rd_ncid, varID, start2, count2, field % array)
#endif

   end subroutine mpas_io_input_field1d_real_time


   subroutine mpas_io_input_field2d_real_time(input_obj, field)
 
      implicit none

      type (io_input_object), intent(in) :: input_obj      
      type (field2dReal), intent(inout) :: field
 
      include 'netcdf.inc'
 
      integer :: nferr
      integer :: varID
      integer, dimension(3) :: start3, count3
 
      start3(1) = field % ioinfo % start(1)
      start3(2) = field % ioinfo % start(2)
      start3(3) = input_obj % time
      count3(1) = field % ioinfo % count(1)
      count3(2) = field % ioinfo % count(2)
      count3(3) = 1
 
#include "input_field2dreal_time.inc"

#if (RKIND == 8)
      nferr = nf_get_vara_double(input_obj % rd_ncid, varID, start3, count3, field % array)
#else
      nferr = nf_get_vara_real(input_obj % rd_ncid, varID, start3, count3, field % array)
#endif

   end subroutine mpas_io_input_field2d_real_time


   subroutine mpas_io_input_field3d_real_time(input_obj, field)
 
      implicit none

      type (io_input_object), intent(in) :: input_obj      
      type (field3dReal), intent(inout) :: field
 
      include 'netcdf.inc'
 
      integer :: nferr
      integer :: varID
      integer, dimension(4) :: start4, count4
 
      start4(1) = field % ioinfo % start(1)
      start4(2) = field % ioinfo % start(2)
      start4(3) = field % ioinfo % start(3)
      start4(4) = input_obj % time
      count4(1) = field % ioinfo % count(1)
      count4(2) = field % ioinfo % count(2)
      count4(3) = field % ioinfo % count(3)
      count4(4) = 1
 
#include "input_field3dreal_time.inc"

#if (RKIND == 8)
      nferr = nf_get_vara_double(input_obj % rd_ncid, varID, start4, count4, field % array)
#else
      nferr = nf_get_vara_real(input_obj % rd_ncid, varID, start4, count4, field % array)
#endif

   end subroutine mpas_io_input_field3d_real_time


   subroutine mpas_io_input_field1d_integer(input_obj, field)
 
      implicit none

      type (io_input_object), intent(in) :: input_obj      
      type (field1dInteger), intent(inout) :: field
 
      include 'netcdf.inc'
 
      integer :: nferr
      integer :: varID
      integer, dimension(1) :: start1, count1
 
      start1(1) = field % ioinfo % start(1)
      count1(1) = field % ioinfo % count(1)
      
#include "input_field1dinteger.inc"

      nferr = nf_get_vara_int(input_obj % rd_ncid, varID, start1, count1, field % array)
 
   end subroutine mpas_io_input_field1d_integer


   subroutine mpas_io_input_field2d_integer(input_obj, field)
 
      implicit none

      type (io_input_object), intent(in) :: input_obj      
      type (field2dInteger), intent(inout) :: field
 
      include 'netcdf.inc'
 
      integer :: nferr
      integer :: varID
      integer, dimension(2) :: start2, count2
 
      start2(1) = field % ioinfo % start(1)
      start2(2) = field % ioinfo % start(2)
      count2(1) = field % ioinfo % count(1)
      count2(2) = field % ioinfo % count(2)

#include "input_field2dinteger.inc"

      nferr = nf_get_vara_int(input_obj % rd_ncid, varID, start2, count2, field % array)

   end subroutine mpas_io_input_field2d_integer


   subroutine mpas_io_input_field1d_integer_time(input_obj, field)
 
      implicit none

      type (io_input_object), intent(in) :: input_obj      
      type (field1dInteger), intent(inout) :: field
 
      include 'netcdf.inc'
 
      integer :: nferr
      integer :: varID
      integer, dimension(2) :: start2, count2
 
      start2(1) = field % ioinfo % start(1)
      start2(2) = input_obj % time
      count2(1) = field % ioinfo % count(1)
      count2(2) = 1
 
#include "input_field1dinteger_time.inc"

      nferr = nf_get_vara_int(input_obj % rd_ncid, varID, start2, count2, field % array)

   end subroutine mpas_io_input_field1d_integer_time


   subroutine mpas_io_input_field0d_char_time(input_obj, field)
 
      implicit none

      type (io_input_object), intent(in) :: input_obj      
      type (field0dChar), intent(inout) :: field
 
      include 'netcdf.inc'
 
      integer :: nferr
      integer :: varID
      integer, dimension(2) :: start1, count1
 
      start1(1) = 1
      count1(1) = 64
      start1(2) = input_obj % time
      count1(2) = 1
 
#include "input_field0dchar_time.inc"

      nferr = nf_get_vara_text(input_obj % rd_ncid, varID, start1, count1, field % scalar)

   end subroutine mpas_io_input_field0d_char_time


   subroutine mpas_io_input_field1d_char_time(input_obj, field)
 
      implicit none

      type (io_input_object), intent(in) :: input_obj      
      type (field1dChar), intent(inout) :: field
 
      include 'netcdf.inc'
 
      integer :: nferr
      integer :: varID
      integer, dimension(3) :: start2, count2
 
      start2(1) = 1
      start2(2) = field % ioinfo % start(1)
      start2(3) = input_obj % time
      count2(1) = 64
      count2(2) = field % ioinfo % count(1)
      count2(3) = 1
 
#include "input_field1dchar_time.inc"

      nferr = nf_get_vara_text(input_obj % rd_ncid, varID, start2, count2, field % array)

   end subroutine mpas_io_input_field1d_char_time


   subroutine mpas_io_input_field0d_char(input_obj, field)
 
      implicit none

      type (io_input_object), intent(in) :: input_obj      
      type (field0dChar), intent(inout) :: field
 
      include 'netcdf.inc'
 
      integer :: nferr
      integer :: varID
      integer, dimension(2) :: start1, count1
 
      start1(1) = 1
      count1(1) = 64
      start1(2) = 1
      count1(2) = 1

#include "input_field0dchar.inc"

      nferr = nf_get_vara_text(input_obj % rd_ncid, varID, start1, count1, field % scalar)
 
   end subroutine mpas_io_input_field0d_char


   subroutine mpas_io_input_field1d_char(input_obj, field)
 
      implicit none

      type (io_input_object), intent(in) :: input_obj      
      type (field1dChar), intent(inout) :: field
 
      include 'netcdf.inc'
 
      integer :: nferr
      integer :: varID
      integer, dimension(2) :: start1, count1
 
      start1(1) = 1
      count1(1) = 64
      start1(2) = field % ioinfo % start(1)
      count1(2) = field % ioinfo % count(1)

      !
      ! Special case: we may want to read the xtime variable across the
      !   time dimension as a 1d array.
      !
      if (trim(field % ioinfo % fieldName) == 'xtime') then
         varID = input_obj % rdVarIDxtime
      end if
 
#include "input_field1dchar.inc"

      nferr = nf_get_vara_text(input_obj % rd_ncid, varID, start1, count1, field % array)
 
   end subroutine mpas_io_input_field1d_char


   subroutine mpas_io_input_finalize(input_obj, dminfo)
 
      implicit none
 
      type (io_input_object), intent(inout) :: input_obj
      type (dm_info), intent(in) :: dminfo

      include 'netcdf.inc'
 
      integer :: nferr
 
      nferr = nf_close(input_obj % rd_ncid)
 
   end subroutine mpas_io_input_finalize
 
end module mpas_io_input
