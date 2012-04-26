module mpas_io_output

   use mpas_grid_types
   use mpas_dmpar
   use mpas_sort
   use mpas_configure
   use mpas_io_streams

   integer, parameter :: OUTPUT = 1
   integer, parameter :: RESTART = 2
   integer, parameter :: SFC = 3
 
   type io_output_object
      character (len=StrKIND) :: filename
      integer :: stream

      integer :: time

      type (MPAS_Stream_type) :: io_stream
   end type io_output_object


   contains

 
   subroutine mpas_output_state_init(output_obj, domain, stream, outputSuffix)

      implicit none

      type (io_output_object), intent(inout) :: output_obj
      type (domain_type), intent(in) :: domain
      character (len=*) :: stream
      character (len=*), optional :: outputSuffix

      character (len=StrKIND) :: tempfilename

      type (block_type), pointer :: block_ptr

      block_ptr => domain % blocklist

      if (trim(stream) == 'OUTPUT') then
         if(present(outputSuffix)) then
            call mpas_insert_string_suffix(config_output_name, outputSuffix, tempfilename)
         else
            tempfilename = config_output_name
         end if
         output_obj % filename = trim(tempfilename)
         output_obj % stream = OUTPUT
      else if (trim(stream) == 'RESTART') then
         if(present(outputSuffix)) then
            call mpas_insert_string_suffix(config_restart_name, outputSuffix, tempfilename)
         else
            tempfilename = config_restart_name
         end if
         output_obj % filename = trim(tempfilename)
         output_obj % stream = RESTART
      else if (trim(stream) == 'SFC') then
         ! Keep filename as whatever was set by the user
         output_obj % stream = SFC
      end if

      ! For now, we assume that a domain consists only of one block,
      !   although in future, work needs to be done to write model state
      !   from many distributed blocks
      call mpas_io_output_init(domain, output_obj, domain % dminfo, &
                          block_ptr % mesh &
                         )

   end subroutine mpas_output_state_init


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


   subroutine mpas_output_state_for_domain(output_obj, domain, itime)
   
      implicit none
   
      type (io_output_object), intent(inout) :: output_obj
      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itime

      integer :: ierr
      integer :: i, j
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
      type (field0dChar) :: char0d
      type (field1dChar) :: char1d

      output_obj % time = itime

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

      call MPAS_writeStream(output_obj % io_stream, output_obj % time, ierr)

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

   end subroutine mpas_output_state_for_domain


   subroutine mpas_output_state_finalize(output_obj, dminfo)

      implicit none

      type (io_output_object), intent(inout) :: output_obj
      type (dm_info), intent(in) :: dminfo

      call mpas_io_output_finalize(output_obj, dminfo)

   end subroutine mpas_output_state_finalize


   subroutine mpas_io_output_init( domain, output_obj, &
                              dminfo, &
                              mesh &
                            )
 
      implicit none
 
      type (domain_type), intent(in) :: domain
      type (io_output_object), intent(inout) :: output_obj
      type (dm_info), intent(in) :: dminfo
      type (mesh_type), intent(in) :: mesh
 
      integer :: nferr, ierr
      integer, dimension(10) :: dimlist
 
      call MPAS_createStream(output_obj % io_stream, trim(output_obj % filename), MPAS_IO_PNETCDF, MPAS_IO_WRITE, 1, nferr)

#include "add_output_fields.inc"

      if (mesh % on_a_sphere) then
         call MPAS_writeStreamAtt(output_obj % io_stream, 'on_a_sphere', 'YES             ', nferr)
      else
         call MPAS_writeStreamAtt(output_obj % io_stream, 'on_a_sphere', 'NO              ', nferr)
      end if
      call MPAS_writeStreamAtt(output_obj % io_stream, 'sphere_radius', mesh % sphere_radius, nferr)

#include "add_output_atts.inc"
 
   end subroutine mpas_io_output_init


   subroutine mpas_io_output_finalize(output_obj, dminfo)
 
      implicit none
 
      type (io_output_object), intent(inout) :: output_obj
      type (dm_info), intent(in) :: dminfo
 
      integer :: nferr
 
      call MPAS_closeStream(output_obj % io_stream, nferr)
 
   end subroutine mpas_io_output_finalize
 
end module mpas_io_output
