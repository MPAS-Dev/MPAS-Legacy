module io_output

   use grid_types
   use dmpar
   use sort
   use configure

   integer, parameter :: OUTPUT = 1
   integer, parameter :: RESTART = 2
 
   type io_output_object
      integer :: wr_ncid
      character (len=1024) :: filename

      integer :: time

      integer :: stream

#include "io_output_obj_decls.inc"

      logical :: validExchangeLists
      type (exchange_list), pointer :: sendCellsList, recvCellsList
      type (exchange_list), pointer :: sendEdgesList, recvEdgesList
      type (exchange_list), pointer :: sendVerticesList, recvVerticesList
      type (exchange_list), pointer :: sendVertLevelsList, recvVertLevelsList
      type (exchange_list), pointer :: sendVertLevelsP1List, recvVertLevelsP1List
   end type io_output_object


   interface io_output_field
      module procedure io_output_field0dReal
      module procedure io_output_field1dReal
      module procedure io_output_field2dReal
      module procedure io_output_field3dReal
      module procedure io_output_field1dInteger
      module procedure io_output_field2dInteger
   end interface io_output_field

   interface io_output_field_time
      module procedure io_output_field0dReal_time
      module procedure io_output_field1dReal_time
      module procedure io_output_field2dReal_time
      module procedure io_output_field3dReal_time
   end interface io_output_field_time
 

   contains

 
   subroutine output_state_init(output_obj, domain, stream)

      implicit none

      type (io_output_object), intent(inout) :: output_obj
      type (domain_type), intent(in) :: domain
      character (len=*) :: stream

      type (block_type), pointer :: block_ptr
#include "output_dim_actual_decls.inc"

      block_ptr => domain % blocklist
      nullify(output_obj % sendCellsList)
      nullify(output_obj % recvCellsList)
      nullify(output_obj % sendEdgesList)
      nullify(output_obj % recvEdgesList)
      nullify(output_obj % sendVerticesList)
      nullify(output_obj % recvVerticesList)
      nullify(output_obj % sendVertLevelsList)
      nullify(output_obj % recvVertLevelsList)
      nullify(output_obj % sendVertLevelsP1List)
      nullify(output_obj % recvVertLevelsP1List)
      output_obj % validExchangeLists = .false.

#include "output_dim_inits.inc"

      call dmpar_sum_int(domain % dminfo, block_ptr % mesh % nCellsSolve, nCellsGlobal) 
      call dmpar_sum_int(domain % dminfo, block_ptr % mesh % nEdgesSolve, nEdgesGlobal) 
      call dmpar_sum_int(domain % dminfo, block_ptr % mesh % nVerticesSolve, nVerticesGlobal) 
      nVertLevelsGlobal = block_ptr % mesh % nVertLevels

      if (trim(stream) == 'OUTPUT') then
         output_obj % filename = trim(config_output_name)
         output_obj % stream = OUTPUT
      else if (trim(stream) == 'RESTART') then
         output_obj % filename = trim(config_restart_name)
         output_obj % stream = RESTART
      end if

      ! For now, we assume that a domain consists only of one block,
      !   although in future, work needs to be done to write model state
      !   from many distributed blocks
      call io_output_init(output_obj, domain % dminfo, &
                          block_ptr % mesh, &
#include "output_dim_actual_args.inc"
                         )
   
   end subroutine output_state_init


   subroutine output_state_for_domain(output_obj, domain, itime)
   
      implicit none
   
      type (io_output_object), intent(inout) :: output_obj
      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itime

      integer :: i, j
      integer :: nCellsGlobal
      integer :: nEdgesGlobal
      integer :: nVerticesGlobal
      integer :: nVertLevelsGlobal
      integer, dimension(:), pointer :: neededCellList
      integer, dimension(:), pointer :: neededEdgeList
      integer, dimension(:), pointer :: neededVertexList
      integer, dimension(:), pointer :: neededVertLevelList
      integer, dimension(:), pointer :: neededVertLevelP1List
      integer, dimension(:,:), pointer :: cellsOnCell, edgesOnCell, verticesOnCell, &
                                          cellsOnEdge, verticesOnEdge, edgesOnEdge, cellsOnVertex, edgesOnVertex
      integer, dimension(:,:), pointer :: cellsOnCell_save, edgesOnCell_save, verticesOnCell_save, &
                                          cellsOnEdge_save, verticesOnEdge_save, edgesOnEdge_save, &
                                          cellsOnVertex_save, edgesOnVertex_save
      type (field1dInteger) :: int1d
      type (field2dInteger) :: int2d
      type (field0dReal) :: real0d
      type (field1dReal) :: real1d
      type (field2dReal) :: real2d
      type (field3dReal) :: real3d

      integer :: i1, i2, i3, i4

      integer, dimension(:), pointer :: super_int1d
      integer, dimension(:,:), pointer :: super_int2d
      real :: super_real0d
      real (kind=RKIND), dimension(:), pointer :: super_real1d
      real (kind=RKIND), dimension(:,:), pointer :: super_real2d
      real (kind=RKIND), dimension(:,:,:), pointer :: super_real3d

      output_obj % time = itime

      allocate(int1d % ioinfo)
      allocate(int2d % ioinfo)
      allocate(real0d % ioinfo)
      allocate(real1d % ioinfo)
      allocate(real2d % ioinfo)
      allocate(real3d % ioinfo)

      call dmpar_sum_int(domain % dminfo, domain % blocklist % mesh % nCellsSolve, nCellsGlobal)
      call dmpar_sum_int(domain % dminfo, domain % blocklist % mesh % nEdgesSolve, nEdgesGlobal)
      call dmpar_sum_int(domain % dminfo, domain % blocklist % mesh % nVerticesSolve, nVerticesGlobal)
      nVertLevelsGlobal = domain % blocklist % mesh % nVertLevels

      allocate(cellsOnCell(domain % blocklist % mesh % maxEdges, domain % blocklist % mesh % nCellsSolve))
      allocate(edgesOnCell(domain % blocklist % mesh % maxEdges, domain % blocklist % mesh % nCellsSolve))
      allocate(verticesOnCell(domain % blocklist % mesh % maxEdges, domain % blocklist % mesh % nCellsSolve))
      allocate(cellsOnEdge(2, domain % blocklist % mesh % nEdgesSolve))
      allocate(verticesOnEdge(2, domain % blocklist % mesh % nEdgesSolve))
      allocate(edgesOnEdge(2 * domain % blocklist % mesh % maxEdges, domain % blocklist % mesh % nEdgesSolve))
      allocate(cellsOnVertex(domain % blocklist % mesh % vertexDegree, domain % blocklist % mesh % nVerticesSolve))
      allocate(edgesOnVertex(domain % blocklist % mesh % vertexDegree, domain % blocklist % mesh % nVerticesSolve))


      !
      ! Convert connectivity information from local to global indices
      !
      do i=1,domain % blocklist % mesh % nCellsSolve
         do j=1,domain % blocklist % mesh % nEdgesOnCell % array(i)
            cellsOnCell(j,i) = domain % blocklist % mesh % indexToCellID % array( &
                                                                           domain % blocklist % mesh % cellsOnCell % array(j,i))
            edgesOnCell(j,i) = domain % blocklist % mesh % indexToEdgeID % array( &
                                                                           domain % blocklist % mesh % edgesOnCell % array(j,i))
            verticesOnCell(j,i) = domain % blocklist % mesh % indexToVertexID % array( &
                                                                           domain % blocklist % mesh % verticesOnCell % array(j,i))
         end do
         do j=domain % blocklist % mesh % nEdgesOnCell % array(i)+1,domain % blocklist % mesh % maxEdges
            cellsOnCell(j,i) = domain % blocklist % mesh % indexToCellID % array( &
                                                                           domain % blocklist % mesh % nEdgesOnCell % array(i))
            edgesOnCell(j,i) = domain % blocklist % mesh % indexToEdgeID % array( &
                                                                           domain % blocklist % mesh % nEdgesOnCell % array(i))
            verticesOnCell(j,i) = domain % blocklist % mesh % indexToVertexID % array( &
                                                                           domain % blocklist % mesh % nEdgesOnCell % array(i))
         end do
      end do
      do i=1,domain % blocklist % mesh % nEdgesSolve
         cellsOnEdge(1,i) = domain % blocklist % mesh % indexToCellID % array(domain % blocklist % mesh % cellsOnEdge % array(1,i))
         cellsOnEdge(2,i) = domain % blocklist % mesh % indexToCellID % array(domain % blocklist % mesh % cellsOnEdge % array(2,i))
         verticesOnEdge(1,i) = domain % blocklist % mesh % indexToVertexID % array( &
                                                                           domain % blocklist % mesh % verticesOnEdge % array(1,i))
         verticesOnEdge(2,i) = domain % blocklist % mesh % indexToVertexID % array( &
                                                                           domain % blocklist % mesh % verticesOnEdge % array(2,i))
         do j=1,domain % blocklist % mesh % nEdgesOnEdge % array(i)
            edgesOnEdge(j,i) = domain % blocklist % mesh % indexToEdgeID % array( &
                                                                           domain % blocklist % mesh % edgesOnEdge % array(j,i))
         end do
         do j=domain % blocklist % mesh % nEdgesOnEdge % array(i)+1,2*domain % blocklist % mesh % maxEdges
            if(domain % blocklist % mesh % nEdgesOnEdge % array(i) .eq. 0) then
               edgesOnEdge(j,i) = domain % blocklist % mesh % nEdgesSolve + 1
            else
               edgesOnEdge(j,i) = domain % blocklist % mesh % indexToEdgeID % array( &
                                                                           domain % blocklist % mesh % nEdgesOnEdge % array(i))
            endif
         end do
      end do
      do i=1,domain % blocklist % mesh % nVerticesSolve
         do j=1,domain % blocklist % mesh % vertexDegree
            cellsOnVertex(j,i) = domain % blocklist % mesh % indexToCellID % array( &
                                                                           domain % blocklist % mesh % cellsOnVertex % array(j,i))
            edgesOnVertex(j,i) = domain % blocklist % mesh % indexToEdgeID % array( &
                                                                           domain % blocklist % mesh % edgesOnVertex % array(j,i))
         end do
      end do

      if (domain % dminfo % my_proc_id == 0) then
         allocate(neededCellList(nCellsGlobal))
         allocate(neededEdgeList(nEdgesGlobal))
         allocate(neededVertexList(nVerticesGlobal))
         allocate(neededVertLevelList(nVertLevelsGlobal))
         allocate(neededVertLevelP1List(nVertLevelsGlobal+1))
         do i=1,nCellsGlobal
            neededCellList(i) = i
         end do
         do i=1,nEdgesGlobal
            neededEdgeList(i) = i
         end do
         do i=1,nVerticesGlobal
            neededVertexList(i) = i
         end do
         do i=1,nVertLevelsGlobal
            neededVertLevelList(i) = i
         end do
         do i=1,nVertLevelsGlobal+1
            neededVertLevelP1List(i) = i
         end do
      else
         allocate(neededCellList(0))
         allocate(neededEdgeList(0))
         allocate(neededVertexList(0))
         allocate(neededVertLevelList(0))
         allocate(neededVertLevelP1List(0))
      end if

      if (.not. output_obj % validExchangeLists) then
         call dmpar_get_owner_list(domain % dminfo, &
                                   domain % blocklist % mesh % nCellsSolve, size(neededCellList), &
                                   domain % blocklist % mesh % indexToCellID % array, neededCellList, &
                                   output_obj % sendCellsList, output_obj % recvCellsList)

         call dmpar_get_owner_list(domain % dminfo, &
                                   domain % blocklist % mesh % nEdgesSolve, size(neededEdgeList), &
                                   domain % blocklist % mesh % indexToEdgeID % array, neededEdgeList, &
                                   output_obj % sendEdgesList, output_obj % recvEdgesList)

         call dmpar_get_owner_list(domain % dminfo, &
                                   domain % blocklist % mesh % nVerticesSolve, size(neededVertexList), &
                                   domain % blocklist % mesh % indexToVertexID % array, neededVertexList, &
                                   output_obj % sendVerticesList, output_obj % recvVerticesList)

         call dmpar_get_owner_list(domain % dminfo, &
                                   size(neededVertLevelList), size(neededVertLevelList), &
                                   neededVertLevelList, neededVertLevelList, &
                                   output_obj % sendVertLevelsList, output_obj % recvVertLevelsList)

         call dmpar_get_owner_list(domain % dminfo, &
                                   size(neededVertLevelP1List), size(neededVertLevelP1List), &
                                   neededVertLevelP1List, neededVertLevelP1List, &
                                   output_obj % sendVertLevelsP1List, output_obj % recvVertLevelsP1List)

         output_obj % validExchangeLists = .true.
      end if

      deallocate(neededCellList)
      deallocate(neededEdgeList)
      deallocate(neededVertexList)

      cellsOnCell_save => domain % blocklist % mesh % cellsOnCell % array
      edgesOnCell_save => domain % blocklist % mesh % edgesOnCell % array
      verticesOnCell_save => domain % blocklist % mesh % verticesOnCell % array
      cellsOnEdge_save => domain % blocklist % mesh % cellsOnEdge % array
      verticesOnEdge_save => domain % blocklist % mesh % verticesOnEdge % array
      edgesOnEdge_save => domain % blocklist % mesh % edgesOnEdge % array
      cellsOnVertex_save => domain % blocklist % mesh % cellsOnVertex % array
      edgesOnVertex_save => domain % blocklist % mesh % edgesOnVertex % array

      domain % blocklist % mesh % cellsOnCell % array => cellsOnCell
      domain % blocklist % mesh % edgesOnCell % array => edgesOnCell
      domain % blocklist % mesh % verticesOnCell % array => verticesOnCell
      domain % blocklist % mesh % cellsOnEdge % array => cellsOnEdge
      domain % blocklist % mesh % verticesOnEdge % array => verticesOnEdge
      domain % blocklist % mesh % edgesOnEdge % array => edgesOnEdge
      domain % blocklist % mesh % cellsOnVertex % array => cellsOnVertex
      domain % blocklist % mesh % edgesOnVertex % array => edgesOnVertex

#include "io_output_fields.inc"

      domain % blocklist % mesh % cellsOnCell % array => cellsOnCell_save
      domain % blocklist % mesh % edgesOnCell % array => edgesOnCell_save
      domain % blocklist % mesh % verticesOnCell % array => verticesOnCell_save
      domain % blocklist % mesh % cellsOnEdge % array => cellsOnEdge_save
      domain % blocklist % mesh % verticesOnEdge % array => verticesOnEdge_save
      domain % blocklist % mesh % edgesOnEdge % array => edgesOnEdge_save
      domain % blocklist % mesh % cellsOnVertex % array => cellsOnVertex_save
      domain % blocklist % mesh % edgesOnVertex % array => edgesOnVertex_save

      deallocate(cellsOnCell)
      deallocate(edgesOnCell)
      deallocate(verticesOnCell)
      deallocate(cellsOnEdge)
      deallocate(verticesOnEdge)
      deallocate(edgesOnEdge)
      deallocate(cellsOnVertex)
      deallocate(edgesOnVertex)

   end subroutine output_state_for_domain


   subroutine output_state_finalize(output_obj, dminfo)

      implicit none

      type (io_output_object), intent(inout) :: output_obj
      type (dm_info), intent(in) :: dminfo

      call io_output_finalize(output_obj, dminfo)

   end subroutine output_state_finalize


   subroutine io_output_init( output_obj, &
                              dminfo, &
                              mesh, &
#include "dim_dummy_args.inc"
                            )
 
      implicit none
 
      include 'netcdf.inc'
 
      type (io_output_object), intent(inout) :: output_obj
      type (dm_info), intent(in) :: dminfo
      type (mesh_type), intent(in) :: mesh
#include "dim_dummy_decls.inc"
 
      integer :: nferr
      integer, dimension(10) :: dimlist
 
      if (dminfo % my_proc_id == 0) then
#ifdef OFFSET64BIT
      nferr = nf_create(trim(output_obj % filename), ior(NF_CLOBBER,NF_64BIT_OFFSET), output_obj % wr_ncid)
#else
      nferr = nf_create(trim(output_obj % filename), NF_CLOBBER, output_obj % wr_ncid)
#endif
 
#include "netcdf_def_dims_vars.inc"

      if (mesh % on_a_sphere) then
         nferr = nf_put_att_text(output_obj % wr_ncid, NF_GLOBAL, 'on_a_sphere', 16, 'YES             ')
      else
         nferr = nf_put_att_text(output_obj % wr_ncid, NF_GLOBAL, 'on_a_sphere', 16, 'NO              ')
      end if
      if (RKIND == 8) then
         nferr = nf_put_att_double(output_obj % wr_ncid, NF_GLOBAL, 'sphere_radius', NF_DOUBLE, 1, mesh % sphere_radius)
      else
         nferr = nf_put_att_real(output_obj % wr_ncid, NF_GLOBAL, 'sphere_radius', NF_FLOAT, 1, mesh % sphere_radius)
      end if
 
      nferr = nf_enddef(output_obj % wr_ncid)
      end if
 
   end subroutine io_output_init


   subroutine io_output_field0dReal(output_obj, field)

      implicit none

      type (io_output_object), intent(in) :: output_obj
      type (field0dReal), intent(inout) :: field

      include 'netcdf.inc'

      integer :: nferr
      integer :: varID
      integer, dimension(1) :: start1, count1

      start1(1) = 1
      count1(1) = 1

#include "output_field0dreal.inc"

#if (RKIND == 8)
      nferr = nf_put_vara_double(output_obj % wr_ncid, varID, start1, count1, field % scalar)
#else
      nferr = nf_put_vara_real(output_obj % wr_ncid, varID, start1, count1, field % scalar)
#endif
 
      nferr = nf_sync(output_obj % wr_ncid)

   end subroutine io_output_field0dReal


   subroutine io_output_field1dReal(output_obj, field)

      implicit none

      type (io_output_object), intent(in) :: output_obj
      type (field1dReal), intent(inout) :: field

      include 'netcdf.inc'

      integer :: nferr
      integer :: varID
      integer, dimension(1) :: start1, count1

      start1(1) = field % ioinfo % start(1)
      count1(1) = field % ioinfo % count(1)

#include "output_field1dreal.inc"

#if (RKIND == 8)
      nferr = nf_put_vara_double(output_obj % wr_ncid, VarID, start1, count1, field % array)
#else
      nferr = nf_put_vara_real(output_obj % wr_ncid, VarID, start1, count1, field % array)
#endif
 
      nferr = nf_sync(output_obj % wr_ncid)

   end subroutine io_output_field1dReal
 
 
   subroutine io_output_field2dReal(output_obj, field)

      implicit none

      type (io_output_object), intent(in) :: output_obj
      type (field2dReal), intent(inout) :: field

      include 'netcdf.inc'

      integer :: nferr
      integer :: varID
      integer, dimension(2) :: start2, count2

      start2(1) = field % ioinfo % start(1)
      start2(2) = field % ioinfo % start(2)
      count2(1) = field % ioinfo % count(1)
      count2(2) = field % ioinfo % count(2)

#include "output_field2dreal.inc"

#if (RKIND == 8)
      nferr = nf_put_vara_double(output_obj % wr_ncid, varID, start2, count2, field % array)
#else
      nferr = nf_put_vara_real(output_obj % wr_ncid, varID, start2, count2, field % array)
#endif
 
      nferr = nf_sync(output_obj % wr_ncid)

   end subroutine io_output_field2dReal
 
 
   subroutine io_output_field3dReal(output_obj, field)

      implicit none

      type (io_output_object), intent(in) :: output_obj
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

#include "output_field3dreal.inc"

#if (RKIND == 8)
      nferr = nf_put_vara_double(output_obj % wr_ncid, varID, start3, count3, field % array)
#else
      nferr = nf_put_vara_real(output_obj % wr_ncid, varID, start3, count3, field % array)
#endif
 
      nferr = nf_sync(output_obj % wr_ncid)

   end subroutine io_output_field3dReal


   subroutine io_output_field0dReal_time(output_obj, field)

      implicit none

      type (io_output_object), intent(in) :: output_obj
      type (field0dReal), intent(inout) :: field

      include 'netcdf.inc'

      integer :: nferr
      integer :: varID
      integer, dimension(1) :: start1, count1

      start1(1) = output_obj % time
      count1(1) = 1

#include "output_field0dreal_time.inc"

#if (RKIND == 8)
      nferr = nf_put_vara_double(output_obj % wr_ncid, varID, start1, count1, field % scalar)
#else
      nferr = nf_put_vara_real(output_obj % wr_ncid, varID, start1, count1, field % scalar)
#endif
 
      nferr = nf_sync(output_obj % wr_ncid)

   end subroutine io_output_field0dReal_time


   subroutine io_output_field1dReal_time(output_obj, field)

      implicit none

      type (io_output_object), intent(in) :: output_obj
      type (field1dReal), intent(inout) :: field

      include 'netcdf.inc'

      integer :: nferr
      integer :: varID
      integer, dimension(2) :: start2, count2

      start2(1) = field % ioinfo % start(1)
      start2(2) = output_obj % time
      count2(1) = field % ioinfo % count(1)
      count2(2) = 1

#include "output_field1dreal_time.inc"

#if (RKIND == 8)
      nferr = nf_put_vara_double(output_obj % wr_ncid, varID, start2, count2, field % array)
#else
      nferr = nf_put_vara_real(output_obj % wr_ncid, varID, start2, count2, field % array)
#endif
 
      nferr = nf_sync(output_obj % wr_ncid)

   end subroutine io_output_field1dReal_time


   subroutine io_output_field2dReal_time(output_obj, field)

      implicit none

      type (io_output_object), intent(in) :: output_obj
      type (field2dReal), intent(inout) :: field

      include 'netcdf.inc'

      integer :: nferr
      integer :: varID
      integer, dimension(3) :: start3, count3

      start3(1) = field % ioinfo % start(1)
      start3(2) = field % ioinfo % start(2)
      start3(3) = output_obj % time
      count3(1) = field % ioinfo % count(1)
      count3(2) = field % ioinfo % count(2)
      count3(3) = 1

#include "output_field2dreal_time.inc"

#if (RKIND == 8)
      nferr = nf_put_vara_double(output_obj % wr_ncid, varID, start3, count3, field % array)
#else
      nferr = nf_put_vara_real(output_obj % wr_ncid, varID, start3, count3, field % array)
#endif
 
      nferr = nf_sync(output_obj % wr_ncid)

   end subroutine io_output_field2dReal_time


   subroutine io_output_field3dReal_time(output_obj, field)

      implicit none

      type (io_output_object), intent(in) :: output_obj
      type (field3dReal), intent(inout) :: field

      include 'netcdf.inc'

      integer :: nferr
      integer :: varID
      integer, dimension(4) :: start4, count4

      start4(1) = field % ioinfo % start(1)
      start4(2) = field % ioinfo % start(2)
      start4(3) = field % ioinfo % start(3)
      start4(4) = output_obj % time
      count4(1) = field % ioinfo % count(1)
      count4(2) = field % ioinfo % count(2)
      count4(3) = field % ioinfo % count(3)
      count4(4) = 1

#include "output_field3dreal_time.inc"

#if (RKIND == 8)
      nferr = nf_put_vara_double(output_obj % wr_ncid, varID, start4, count4, field % array)
#else
      nferr = nf_put_vara_real(output_obj % wr_ncid, varID, start4, count4, field % array)
#endif
 
      nferr = nf_sync(output_obj % wr_ncid)

   end subroutine io_output_field3dReal_time


   subroutine io_output_field1dInteger(output_obj, field)

      implicit none

      type (io_output_object), intent(in) :: output_obj
      type (field1dInteger), intent(inout) :: field

      include 'netcdf.inc'

      integer :: nferr
      integer :: varID
      integer, dimension(1) :: start1, count1

      start1(1) = field % ioinfo % start(1)
      count1(1) = field % ioinfo % count(1)

#include "output_field1dinteger.inc"

      nferr = nf_put_vara_int(output_obj % wr_ncid, varID, start1, count1, field % array)
 
      nferr = nf_sync(output_obj % wr_ncid)

   end subroutine io_output_field1dInteger


   subroutine io_output_field2dInteger(output_obj, field)

      implicit none

      type (io_output_object), intent(in) :: output_obj
      type (field2dInteger), intent(inout) :: field

      include 'netcdf.inc'

      integer :: nferr
      integer :: varID
      integer, dimension(2) :: start2, count2

      start2(1) = field % ioinfo % start(1)
      start2(2) = field % ioinfo % start(2)
      count2(1) = field % ioinfo % count(1)
      count2(2) = field % ioinfo % count(2)

#include "output_field2dinteger.inc"

      nferr = nf_put_vara_int(output_obj % wr_ncid, varID, start2, count2, field % array)
 
      nferr = nf_sync(output_obj % wr_ncid)

   end subroutine io_output_field2dInteger


   subroutine io_output_finalize(output_obj, dminfo)
 
      implicit none
 
      include 'netcdf.inc'

      type (io_output_object), intent(inout) :: output_obj
      type (dm_info), intent(in) :: dminfo
 
      integer :: nferr
 
      if (dminfo % my_proc_id == 0) then
      nferr = nf_close(output_obj % wr_ncid)
      end if
 
   end subroutine io_output_finalize
 
end module io_output
 
