module zoltan_interface
   use zoltan

   implicit none

   include 'mpif.h'

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Data for reordering cells
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer :: numCells
   integer, dimension(:), pointer :: cellIDs
   integer :: geomDim
   real (kind=RKIND), dimension(:), pointer :: cellCoordX, cellCoordY, cellCoordZ
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Data for reordering edges
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer :: numEdges
   integer, dimension(:), pointer :: edgeIDs
   real (kind=RKIND), dimension(:), pointer :: edgeCoordX, edgeCoordY, edgeCoordZ  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Data for reordering vertices
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer :: numVerts
   integer, dimension(:), pointer :: vertIDs
   real (kind=RKIND), dimension(:), pointer :: vertCoordX, vertCoordY, vertCoordZ  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Perhaps not necessary, but implemented in case it helps
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine zoltanStart()  

      integer(Zoltan_INT) :: error
      real(Zoltan_FLOAT) :: version

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Body of subroutine
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      error = Zoltan_Initialize(version)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
   end subroutine
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine zoltanOrderLocHSFC_Cells(in_numcells,in_cellIDs,in_geomDim,in_cellX, &
                                       in_cellY, in_cellZ)
      implicit none

      integer :: in_numcells
      integer, dimension(:), pointer :: in_cellIDs
      integer :: in_geomDim
      real (kind=RKIND), dimension(:), pointer :: in_cellX, in_cellY, in_cellZ

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! local variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      type(Zoltan_Struct), pointer :: zz_obj
      integer(ZOLTAN_INT) :: ierr

      integer :: numGidEntries, i
      integer(ZOLTAN_INT), allocatable :: global_ids(:), permIndices(:),permGIDs(:)
      real(kind=RKIND), allocatable :: permXs(:),permYs(:),permZs(:)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Body of subroutine
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      numCells = in_numcells
      cellIDs => in_cellIDs
      geomDim = in_geomDim
      cellCoordX => in_cellX
      cellCoordY => in_cellY
      cellCoordZ => in_cellZ

      nullify(zz_obj)
      zz_obj => Zoltan_Create(MPI_COMM_SELF)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! General Zoltan Parameters
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ierr = Zoltan_Set_Param(zz_obj, "ORDER_METHOD", "LOCAL_HSFC")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! register query functions
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_OBJ_FN_TYPE,zqfNumCells)
      ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_OBJ_LIST_FN_TYPE,zqfGetCells)
      ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_GEOM_FN_TYPE,zqfGeomDim)
      ierr =  Zoltan_Set_Fn(zz_obj, ZOLTAN_GEOM_FN_TYPE, zqfGetCellGeom)

      numGidEntries=1

      allocate(global_ids(numCells))
      allocate(permIndices(numCells))
      allocate(permGIDs(numCells))
      allocate(permXs(numCells))
      allocate(permYs(numCells))
      allocate(permZs(numCells))

      !! MMW: There might be a way to use cellIDs directly
      do i=1,numCells
        global_ids(i) = cellIDs(i)
      end do

      ierr = Zoltan_Order(zz_obj, numGidEntries, numCells, global_ids, permIndices);

      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !! This is necessary for now until we fix a small bug in Zoltan_Order
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=1,numCells
        permGIDs(i) = global_ids(permIndices(i)+1)
        permXs(i) = cellCoordX(permIndices(i)+1)
        permYs(i) = cellCoordY(permIndices(i)+1)
        permZs(i) = cellCoordZ(permIndices(i)+1)
      end do

      !!do i=1,numCells
      !!   write(*,*) global_ids(i), permGIDs(i)
      !!end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Actually change the ordering of the cells
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=1,numCells
        cellIDs(i) = permGIDs(i)
        cellCoordX(i) = permXs(i)
        cellCoordY(i) = permYs(i)
        cellCoordZ(i) = permZs(i)
      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!

      deallocate(global_ids)
      deallocate(permIndices)
      deallocate(permGIDs)
      deallocate(permXs)
      deallocate(permYs)
      deallocate(permZs)

      call Zoltan_Destroy(zz_obj)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end subroutine zoltanOrderLocHSFC_Cells
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! zoltan query function:
   !!    Returns number of cells
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function zqfNumCells(data, ierr)

      ! Local declarations
      integer(ZOLTAN_INT), intent(in) :: data(*)
      integer(ZOLTAN_INT), intent(out) :: ierr

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      zqfNumCells = numCells
      ierr = ZOLTAN_OK
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end function zqfNumCells
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! zoltan query function: 
   !!    Returns lists of Cell IDs
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine zqfGetCells (data, num_gid_entries, num_lid_entries, global_ids, &
                           local_ids, wgt_dim, obj_wgts, ierr)
     integer(ZOLTAN_INT), intent(in) :: data(*)
     integer(ZOLTAN_INT), intent(in) :: num_gid_entries, num_lid_entries
     integer(ZOLTAN_INT), intent(out) :: global_ids(*), local_ids(*)
     integer(ZOLTAN_INT), intent(in) :: wgt_dim 
     real(ZOLTAN_FLOAT), intent(out) :: obj_wgts(*)
     integer(ZOLTAN_INT), intent(out) :: ierr

     ! local declarations
     integer :: i

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do i= 1, numCells
       global_ids(i) = cellIDs(i)
       local_ids(i) = i
     end do

     ierr = ZOLTAN_OK
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end subroutine zqfGetCells
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Zoltan Query Function:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function zqfGeomDim(data, ierr)
      !use zoltan
      implicit none
      integer(ZOLTAN_INT), intent(in) :: data(*)
      integer(ZOLTAN_INT) :: ierr

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      zqfGeomDim = geomDim
      ierr = ZOLTAN_OK
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end function zqfGeomDim
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Zoltan Query Function:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine zqfGetCellGeom(data, num_gid_entries, num_lid_entries, global_id, &
                             local_id, geom_vec, ierr)
      !use zoltan
      implicit none

      integer(ZOLTAN_INT), intent(in) :: data(*)
      integer(ZOLTAN_INT), intent(in) :: num_gid_entries, num_lid_entries
      integer(ZOLTAN_INT), intent(in) :: global_id, local_id
      real(ZOLTAN_DOUBLE), intent(out) :: geom_vec(*)
      integer(ZOLTAN_INT), intent(out) :: ierr

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Assuming geom_dim is 3
      geom_vec(1) = cellCoordX(local_id)
      geom_vec(2) = cellCoordY(local_id)
      geom_vec(3) = cellCoordZ(local_id)

      ierr = ZOLTAN_OK
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end subroutine zqfGetCellGeom
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! The ordering functions should perhaps be refactored so that there
   !! are not separate functions for cells, edges, and vertices
   !! Not sure if this is worth it with the additional conditionals that would 
   !! be required. 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine zoltanOrderLocHSFC_Edges(in_numedges,in_edgeIDs,in_geomDim,in_edgeX, &
                                       in_edgeY, in_edgeZ)
      implicit none

      integer :: in_numedges
      integer, dimension(:), pointer :: in_edgeIDs
      integer :: in_geomDim
      real (kind=RKIND), dimension(:), pointer :: in_edgeX, in_edgeY, in_edgeZ

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! local variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      type(Zoltan_Struct), pointer :: zz_obj
      integer(ZOLTAN_INT) :: ierr

      integer :: numGidEntries, i
      integer(ZOLTAN_INT), allocatable :: global_ids(:), permIndices(:),permGIDs(:)
      real(kind=RKIND), allocatable :: permXs(:),permYs(:),permZs(:)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Body of subroutine
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      numEdges = in_numedges
      edgeIDs => in_edgeIDs
      geomDim = in_geomDim
      edgeCoordX => in_edgeX
      edgeCoordY => in_edgeY
      edgeCoordZ => in_edgeZ

      nullify(zz_obj)
      zz_obj => Zoltan_Create(MPI_COMM_SELF)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! General Zoltan Parameters
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ierr = Zoltan_Set_Param(zz_obj, "ORDER_METHOD", "LOCAL_HSFC")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! register query functions
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_OBJ_FN_TYPE,zqfNumEdges)
      ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_OBJ_LIST_FN_TYPE,zqfGetEdges)
      ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_GEOM_FN_TYPE,zqfGeomDim)
      ierr =  Zoltan_Set_Fn(zz_obj, ZOLTAN_GEOM_FN_TYPE, zqfGetEdgeGeom)

      numGidEntries=1

      allocate(global_ids(numEdges))
      allocate(permIndices(numEdges))
      allocate(permGIDs(numEdges))
      allocate(permXs(numEdges))
      allocate(permYs(numEdges))
      allocate(permZs(numEdges))

      !! MMW: There might be a way to use edgeIDs directly
      do i=1,numEdges
        global_ids(i) = edgeIDs(i)
      end do

      ierr = Zoltan_Order(zz_obj, numGidEntries, numEdges, global_ids, permIndices);

      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !! This is necessary for now until we fix a small bug in Zoltan_Order
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=1,numEdges
        permGIDs(i) = global_ids(permIndices(i)+1)
        permXs(i) = edgeCoordX(permIndices(i)+1)
        permYs(i) = edgeCoordY(permIndices(i)+1)
        permZs(i) = edgeCoordZ(permIndices(i)+1)
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Actually change the ordering of the edges
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=1,numEdges
        edgeIDs(i) = permGIDs(i)
        edgeCoordX(i) = permXs(i)
        edgeCoordY(i) = permYs(i)
        edgeCoordZ(i) = permZs(i)
      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!

      deallocate(global_ids)
      deallocate(permIndices)
      deallocate(permGIDs)
      deallocate(permXs)
      deallocate(permYs)
      deallocate(permZs)

      call Zoltan_Destroy(zz_obj)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end subroutine zoltanOrderLocHSFC_Edges
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! zoltan query function:
   !!    Returns number of edges
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function zqfNumEdges(data, ierr)
      ! Local declarations
      integer(ZOLTAN_INT), intent(in) :: data(*)
      integer(ZOLTAN_INT), intent(out) :: ierr

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      zqfNumEdges = numEdges
      ierr = ZOLTAN_OK
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end function zqfNumEdges
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! zoltan query function: 
   !!    Returns lists of Edge IDs
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine zqfGetEdges (data, num_gid_entries, num_lid_entries, global_ids, &
                           local_ids, wgt_dim, obj_wgts, ierr)
     integer(ZOLTAN_INT), intent(in) :: data(*)
     integer(ZOLTAN_INT), intent(in) :: num_gid_entries, num_lid_entries
     integer(ZOLTAN_INT), intent(out) :: global_ids(*), local_ids(*)
     integer(ZOLTAN_INT), intent(in) :: wgt_dim 
     real(ZOLTAN_FLOAT), intent(out) :: obj_wgts(*)
     integer(ZOLTAN_INT), intent(out) :: ierr

     ! local declarations
     integer :: i

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do i= 1, numEdges
       global_ids(i) = edgeIDs(i)
       local_ids(i) = i
     end do

     ierr = ZOLTAN_OK
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end subroutine zqfGetEdges
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Zoltan Query Function:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine zqfGetEdgeGeom(data, num_gid_entries, num_lid_entries, global_id, &
                             local_id, geom_vec, ierr)
      !use zoltan
      implicit none

      integer(ZOLTAN_INT), intent(in) :: data(*)
      integer(ZOLTAN_INT), intent(in) :: num_gid_entries, num_lid_entries
      integer(ZOLTAN_INT), intent(in) :: global_id, local_id
      real(ZOLTAN_DOUBLE), intent(out) :: geom_vec(*)
      integer(ZOLTAN_INT), intent(out) :: ierr

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Assuming geom_dim is 3
      geom_vec(1) = edgeCoordX(local_id)
      geom_vec(2) = edgeCoordY(local_id)
      geom_vec(3) = edgeCoordZ(local_id)

      ierr = ZOLTAN_OK
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end subroutine zqfGetEdgeGeom
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine zoltanOrderLocHSFC_Verts(in_numverts,in_vertIDs,in_geomDim,in_vertX, &
                                       in_vertY, in_vertZ)
      implicit none

      integer :: in_numverts
      integer, dimension(:), pointer :: in_vertIDs
      integer :: in_geomDim
      real (kind=RKIND), dimension(:), pointer :: in_vertX, in_vertY, in_vertZ

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! local variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      type(Zoltan_Struct), pointer :: zz_obj
      integer(ZOLTAN_INT) :: ierr

      integer :: numGidEntries, i
      integer(ZOLTAN_INT), allocatable :: global_ids(:), permIndices(:),permGIDs(:)
      real(kind=RKIND), allocatable :: permXs(:),permYs(:),permZs(:)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Body of subroutine
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      numVerts = in_numverts
      vertIDs => in_vertIDs
      geomDim = in_geomDim
      vertCoordX => in_vertX
      vertCoordY => in_vertY
      vertCoordZ => in_vertZ

      nullify(zz_obj)
      zz_obj => Zoltan_Create(MPI_COMM_SELF)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! General Zoltan Parameters
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ierr = Zoltan_Set_Param(zz_obj, "ORDER_METHOD", "LOCAL_HSFC")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! register query functions
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_OBJ_FN_TYPE,zqfNumVerts)
      ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_OBJ_LIST_FN_TYPE,zqfGetVerts)
      ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_GEOM_FN_TYPE,zqfGeomDim)
      ierr =  Zoltan_Set_Fn(zz_obj, ZOLTAN_GEOM_FN_TYPE, zqfGetVertGeom)

      numGidEntries=1

      allocate(global_ids(numVerts))
      allocate(permIndices(numVerts))
      allocate(permGIDs(numVerts))
      allocate(permXs(numVerts))
      allocate(permYs(numVerts))
      allocate(permZs(numVerts))

      !! MMW: There might be a way to use vertIDs directly
      do i=1,numVerts
        global_ids(i) = vertIDs(i)
      end do

      ierr = Zoltan_Order(zz_obj, numGidEntries, numVerts, global_ids, permIndices);

      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !! This is necessary for now until we fix a small bug in Zoltan_Order
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=1,numVerts
        permGIDs(i) = global_ids(permIndices(i)+1)
        permXs(i) = vertCoordX(permIndices(i)+1)
        permYs(i) = vertCoordY(permIndices(i)+1)
        permZs(i) = vertCoordZ(permIndices(i)+1)
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Actually change the ordering of the verts
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=1,numVerts
        vertIDs(i) = permGIDs(i)
        vertCoordX(i) = permXs(i)
        vertCoordY(i) = permYs(i)
        vertCoordZ(i) = permZs(i)
      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!

      deallocate(global_ids)
      deallocate(permIndices)
      deallocate(permGIDs)
      deallocate(permXs)
      deallocate(permYs)
      deallocate(permZs)

      call Zoltan_Destroy(zz_obj)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end subroutine zoltanOrderLocHSFC_Verts
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! zoltan query function:
   !!    Returns number of verts
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function zqfNumVerts(data, ierr)

      ! Local declarations
      integer(ZOLTAN_INT), intent(in) :: data(*)
      integer(ZOLTAN_INT), intent(out) :: ierr

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      zqfNumVerts = numVerts
      ierr = ZOLTAN_OK
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end function zqfNumVerts
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! zoltan query function: 
   !!    Returns lists of Vert IDs
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine zqfGetVerts (data, num_gid_entries, num_lid_entries, global_ids, &
                           local_ids, wgt_dim, obj_wgts, ierr)

     integer(ZOLTAN_INT), intent(in) :: data(*)
     integer(ZOLTAN_INT), intent(in) :: num_gid_entries, num_lid_entries
     integer(ZOLTAN_INT), intent(out) :: global_ids(*), local_ids(*)
     integer(ZOLTAN_INT), intent(in) :: wgt_dim 
     real(ZOLTAN_FLOAT), intent(out) :: obj_wgts(*)
     integer(ZOLTAN_INT), intent(out) :: ierr

     ! local declarations
     integer :: i

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do i= 1, numVerts
       global_ids(i) = vertIDs(i)
       local_ids(i) = i
     end do

     ierr = ZOLTAN_OK
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end subroutine zqfGetVerts
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Zoltan Query Function:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine zqfGetVertGeom(data, num_gid_entries, num_lid_entries, global_id, &
                             local_id, geom_vec, ierr)
      !use zoltan
      implicit none

      integer(ZOLTAN_INT), intent(in) :: data(*)
      integer(ZOLTAN_INT), intent(in) :: num_gid_entries, num_lid_entries
      integer(ZOLTAN_INT), intent(in) :: global_id, local_id
      real(ZOLTAN_DOUBLE), intent(out) :: geom_vec(*)
      integer(ZOLTAN_INT), intent(out) :: ierr

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Assuming geom_dim is 3
      geom_vec(1) = vertCoordX(local_id)
      geom_vec(2) = vertCoordY(local_id)
      geom_vec(3) = vertCoordZ(local_id)

      ierr = ZOLTAN_OK
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end subroutine zqfGetVertGeom
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




end module zoltan_interface
