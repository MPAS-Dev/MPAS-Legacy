module global_diagnostics

   use grid_types
   use configure
   use constants
   use dmpar

   implicit none
   save
   public

   contains

   subroutine computeGlobalDiagnostics(dminfo, state, grid, timeIndex, dt)

      ! Note: this routine assumes that there is only one block per processor. No looping
      ! is preformed over blocks.
      ! dminfo is the domain info needed for global communication
      ! state contains the state variables needed to compute global diagnostics
      ! grid conains the meta data about the grid
      ! timeIndex is the current time step counter
      ! dt is the duration of each time step

      ! Sums of variables at vertices are not weighted by thickness (since h is not known at
      !    vertices as it is at cell centers and at edges).

      implicit none

      type (dm_info), intent(in) :: dminfo
      type (grid_state), intent(inout) :: state
      type (grid_meta), intent(in) :: grid
      integer, intent(in) :: timeIndex
      real (kind=RKIND), intent(in) :: dt

      integer :: nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve, nCellsGlobal, nEdgesGlobal, nVerticesGlobal
      real (kind=RKIND) ::  areaCellGlobal, areaEdgeGlobal, areaTriangleGlobal
      real (kind=RKIND), dimension(:), pointer ::  areaCell, dcEdge, dvEdge, areaTriangle, areaEdge
      real (kind=RKIND), dimension(:,:), pointer :: h, u, v, h_edge, circulation, vorticity, ke, pv_edge, pv_vertex, &
         pv_cell, gradPVn, gradPVt, zmid, zbot, pmid, pbot, MontPot
      real (kind=RKIND) :: volumeCellGlobal, volumeEdgeGlobal, CFLNumberGlobal
      real (kind=RKIND) ::  localCFL, localSum
      integer :: elementIndex, variableIndex, nVariables, nSums, nMaxes, nMins
      integer :: timeLevel

      integer, parameter :: kMaxVariables = 1024 ! this must be a little more than double the number of variables to be reduced

      real (kind=RKIND), dimension(kMaxVariables) :: sums, mins, maxes, averages, verticalSumMins, verticalSumMaxes, reductions

      integer :: fileID

      nVertLevels = grid % nVertLevels
      nCellsSolve = grid % nCellsSolve
      nEdgesSolve = grid % nEdgesSolve
      nVerticesSolve = grid % nVerticesSolve

      areaCell => grid % areaCell % array
      dcEdge => grid % dcEdge % array
      dvEdge => grid % dvEdge % array
      areaTriangle => grid % areaTriangle % array
      allocate(areaEdge(1:nEdgesSolve))
      areaEdge = dcEdge(1:nEdgesSolve)*dvEdge(1:nEdgesSolve)

      h => state % h % array
      u => state % u % array
      v => state % v % array
      h_edge => state % h_edge % array
      circulation => state % circulation % array
      vorticity => state % vorticity % array
      ke => state % ke % array
      pv_edge => state % pv_edge % array
      pv_vertex => state % pv_vertex % array
      pv_cell => state % pv_cell % array
      gradPVn => state % gradPVn % array
      gradPVt => state % gradPVt % array
      zmid => state % zmid % array
      zbot => state % zbot % array
      pmid => state % pmid % array
      pbot => state % pbot % array
      MontPot => state % MontPot % array

      variableIndex = 0
      ! h
      variableIndex = variableIndex + 1
      call computeFieldAreaWeightedLocalStats(dminfo, nVertLevels, nCellsSolve, areaCell(1:nCellsSolve), h(:,1:nCellsSolve), &
        sums(variableIndex), mins(variableIndex), maxes(variableIndex), verticalSumMins(variableIndex), verticalSumMaxes(variableIndex))

      ! u
      variableIndex = variableIndex + 1
      call computeFieldVolumeWeightedLocalStats(dminfo, nVertLevels, nEdgesSolve, areaEdge(1:nEdgesSolve), h_edge(:,1:nEdgesSolve), &
        u(:,1:nEdgesSolve), sums(variableIndex), mins(variableIndex), maxes(variableIndex), verticalSumMins(variableIndex), &
        verticalSumMaxes(variableIndex))

      ! v
      variableIndex = variableIndex + 1
      call computeFieldVolumeWeightedLocalStats(dminfo, nVertLevels, nEdgesSolve, areaEdge(1:nEdgesSolve), h_edge(:,1:nEdgesSolve), &
        v(:,1:nEdgesSolve), sums(variableIndex), mins(variableIndex), maxes(variableIndex), verticalSumMins(variableIndex), &
        verticalSumMaxes(variableIndex))

      ! h_edge
      variableIndex = variableIndex + 1
      call computeFieldAreaWeightedLocalStats(dminfo, nVertLevels, nEdgesSolve, areaEdge(1:nEdgesSolve), h_edge(:,1:nEdgesSolve), &
        sums(variableIndex), mins(variableIndex), maxes(variableIndex), verticalSumMins(variableIndex), verticalSumMaxes(variableIndex))

      ! circulation
      variableIndex = variableIndex + 1
      call computeFieldLocalStats(dminfo, nVertLevels, nVerticesSolve, circulation(:,1:nVerticesSolve), &
        sums(variableIndex), mins(variableIndex), maxes(variableIndex), verticalSumMins(variableIndex), verticalSumMaxes(variableIndex))

      ! vorticity
      variableIndex = variableIndex + 1
      call computeFieldAreaWeightedLocalStats(dminfo, nVertLevels, nVerticesSolve, areaTriangle(1:nVerticesSolve), &
        vorticity(:,1:nVerticesSolve), sums(variableIndex), mins(variableIndex), maxes(variableIndex), &
        verticalSumMins(variableIndex), verticalSumMaxes(variableIndex))

      ! ke
      variableIndex = variableIndex + 1
      call computeFieldVolumeWeightedLocalStats(dminfo, nVertLevels, nCellsSolve, areaCell(1:nCellsSolve), h(:,1:nCellsSolve), &
        ke(:,1:nCellsSolve), sums(variableIndex), mins(variableIndex), maxes(variableIndex), verticalSumMins(variableIndex), &
        verticalSumMaxes(variableIndex))

      ! pv_edge
      variableIndex = variableIndex + 1
      call computeFieldVolumeWeightedLocalStats(dminfo, nVertLevels, nEdgesSolve, areaEdge(1:nEdgesSolve), h_edge(:,1:nEdgesSolve), &
        pv_edge(:,1:nEdgesSolve), sums(variableIndex), mins(variableIndex), maxes(variableIndex), verticalSumMins(variableIndex), &
        verticalSumMaxes(variableIndex))

      ! pv_vertex
      variableIndex = variableIndex + 1
      call computeFieldAreaWeightedLocalStats(dminfo, nVertLevels, nVerticesSolve, areaTriangle(1:nVerticesSolve), &
        pv_vertex(:,1:nVerticesSolve), sums(variableIndex), mins(variableIndex), maxes(variableIndex), &
        verticalSumMins(variableIndex), verticalSumMaxes(variableIndex))

      ! pv_cell
      variableIndex = variableIndex + 1
      call computeFieldVolumeWeightedLocalStats(dminfo, nVertLevels, nCellsSolve, areaCell(1:nCellsSolve), h(:,1:nCellsSolve), &
        pv_cell(:,1:nCellsSolve), sums(variableIndex), mins(variableIndex), maxes(variableIndex), verticalSumMins(variableIndex), &
        verticalSumMaxes(variableIndex))

      ! gradPVn
      variableIndex = variableIndex + 1
      call computeFieldVolumeWeightedLocalStats(dminfo, nVertLevels, nEdgesSolve, areaEdge(1:nEdgesSolve), h_edge(:,1:nEdgesSolve), &
        gradPVn(:,1:nEdgesSolve), sums(variableIndex), mins(variableIndex), maxes(variableIndex), verticalSumMins(variableIndex), &
        verticalSumMaxes(variableIndex))

      ! gradPVt
      variableIndex = variableIndex + 1
      call computeFieldVolumeWeightedLocalStats(dminfo, nVertLevels, nEdgesSolve, areaEdge(1:nEdgesSolve), h_edge(:,1:nEdgesSolve), &
        gradPVt(:,1:nEdgesSolve), sums(variableIndex), mins(variableIndex), maxes(variableIndex), verticalSumMins(variableIndex), &
        verticalSumMaxes(variableIndex))

      ! zmid
      variableIndex = variableIndex + 1
      call computeFieldAreaWeightedLocalStats(dminfo, nVertLevels, nCellsSolve, areaCell(1:nCellsSolve), &
        zmid(:,1:nCellsSolve), sums(variableIndex), mins(variableIndex), maxes(variableIndex), verticalSumMins(variableIndex), &
        verticalSumMaxes(variableIndex))

      ! zbot
      variableIndex = variableIndex + 1
      call computeFieldAreaWeightedLocalStats(dminfo, nVertLevels, nCellsSolve, areaCell(1:nCellsSolve), &
        zbot(:,1:nCellsSolve), sums(variableIndex), mins(variableIndex), maxes(variableIndex), verticalSumMins(variableIndex), &
        verticalSumMaxes(variableIndex))

      ! pmid
      variableIndex = variableIndex + 1
      call computeFieldVolumeWeightedLocalStats(dminfo, nVertLevels, nCellsSolve, areaCell(1:nCellsSolve), h(:,1:nCellsSolve), &
        pmid(:,1:nCellsSolve), sums(variableIndex), mins(variableIndex), maxes(variableIndex), verticalSumMins(variableIndex), &
        verticalSumMaxes(variableIndex))

      ! pbot
      variableIndex = variableIndex + 1
      call computeFieldVolumeWeightedLocalStats(dminfo, nVertLevels, nCellsSolve, areaCell(1:nCellsSolve), h(:,1:nCellsSolve), &
        pbot(:,1:nCellsSolve), sums(variableIndex), mins(variableIndex), maxes(variableIndex), verticalSumMins(variableIndex), &
        verticalSumMaxes(variableIndex))

      ! MontPot
      variableIndex = variableIndex + 1
      call computeFieldVolumeWeightedLocalStats(dminfo, nVertLevels, nCellsSolve, areaCell(1:nCellsSolve), h(:,1:nCellsSolve), &
        MontPot(:,1:nCellsSolve), sums(variableIndex), mins(variableIndex), maxes(variableIndex), verticalSumMins(variableIndex), &
        verticalSumMaxes(variableIndex))

      nVariables = variableIndex
      nSums = nVariables
      nMins = nVariables
      nMaxes = nVariables

      nSums = nSums + 1
      sums(nSums) = sum(areaCell(1:nCellsSolve))

      nSums = nSums + 1
      sums(nSums) = sum(dcEdge(1:nEdgesSolve)*dvEdge(1:nEdgesSolve))

      nSums = nSums + 1
      sums(nSums) = sum(areaTriangle(1:nVerticesSolve))

      nSums = nSums + 1
      sums(nSums) = nCellsSolve

      nSums = nSums + 1
      sums(nSums) = nEdgesSolve

      nSums = nSums + 1
      sums(nSums) = nVerticesSolve

      localCFL = 0.0
      do elementIndex = 1,nEdgesSolve
         localCFL = max(localCFL, maxval(dt*u(:,elementIndex)/dcEdge(elementIndex)))
      end do
      nMaxes = nMaxes + 1
      maxes(nMaxes) = localCFL

      mins(nMins+1:nMins+nVariables) = verticalSumMins(1:nVariables)
      nMins = nMins + nVariables
      maxes(nMaxes+1:nMaxes+nVariables) = verticalSumMaxes(1:nVariables)
      nMaxes = nMaxes + nVariables

      ! global reduction of the 5 arrays (packed into 3 to minimize global communication)
      call dmpar_sum_real_array(dminfo, nSums, sums(1:nSums), reductions(1:nSums))
      sums(1:nVariables) = reductions(1:nVariables)
      areaCellGlobal = reductions(nVariables+1)
      areaEdgeGlobal = reductions(nVariables+2)
      areaTriangleGlobal = reductions(nVariables+3)
      nCellsGlobal = int(reductions(nVariables+4))
      nEdgesGlobal = int(reductions(nVariables+5))
      nVerticesGlobal = int(reductions(nVariables+6))

      call dmpar_min_real_array(dminfo, nMins, mins(1:nMins), reductions(1:nMins))
      mins(1:nVariables) = reductions(1:nVariables)
      verticalSumMins(1:nVariables) = reductions(nMins-nVariables+1:nMins)

      call dmpar_max_real_array(dminfo, nMaxes, maxes(1:nMaxes), reductions(1:nMaxes))
      maxes(1:nVariables) = reductions(1:nVariables)
      CFLNumberGlobal = reductions(nVariables+1)
      verticalSumMaxes(1:nVariables) = reductions(nMaxes-nVariables+1:nMaxes)

      volumeCellGlobal = sums(1)
      volumeEdgeGlobal = sums(4)
      ! compute the averages (slightly different depending on how the sum was computed)
      variableIndex = 0
      ! h
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaCellGlobal*nVertLevels)

      ! u
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeEdgeGlobal

      ! v
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeEdgeGlobal

      ! h_edge
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaEdgeGlobal*nVertLevels)

      ! circulation
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(nVerticesGlobal*nVertLevels)

      ! vorticity
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaTriangleGlobal*nVertLevels)

      ! ke
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeCellGlobal

      ! pv_edge
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeEdgeGlobal

      ! pv_vertex
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaTriangleGlobal*nVertLevels)

      ! pv_cell
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeCellGlobal

      ! gradPVn
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeEdgeGlobal

      ! gradPVt
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeEdgeGlobal

      ! zmid
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaCellGlobal*nVertLevels)

      ! zbot
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaCellGlobal*nVertLevels)

      ! pmid
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeCellGlobal

      ! pbot
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeCellGlobal

      ! MontPot
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeCellGlobal

      ! write out the data to files
      if (dminfo % my_proc_id == IO_NODE) then
         fileID = getFreeUnit()
         open(fileID,file='stats_min.txt',ACCESS='append')
            write (fileID,'(100es24.16)') mins(1:nVariables)
         close (fileID)
         open(fileID,file='stats_max.txt',ACCESS='append')
            write (fileID,'(100es24.16)') maxes(1:nVariables)
         close (fileID)
         open(fileID,file='stats_sum.txt',ACCESS='append')
            write (fileID,'(100es24.16)') sums(1:nVariables)
         close (fileID)
         open(fileID,file='stats_avg.txt',ACCESS='append')
            write (fileID,'(100es24.16)') averages(1:nVariables)
         close (fileID)
         open(fileID,file='stats_time.txt',ACCESS='append')
            write (fileID,'(i5,100es24.16)') timeIndex, &
               state % xtime % scalar, dt, &
               CFLNumberGlobal
         close (fileID)
         open(fileID,file='stats_colmin.txt',ACCESS='append')
            write (fileID,'(100es24.16)') verticalSumMins(1:nVariables)
         close (fileID)
         open(fileID,file='stats_colmax.txt',ACCESS='append')
            write (fileID,'(100es24.16)') verticalSumMaxes(1:nVariables)
         close (fileID)
      end if

      state % areaCellGlobal % scalar = areaCellGlobal
      state % areaEdgeGlobal % scalar = areaEdgeGlobal
      state % areaTriangleGlobal % scalar = areaTriangleGlobal

      state % volumeCellGlobal % scalar = volumeCellGlobal
      state % volumeEdgeGlobal % scalar = volumeEdgeGlobal
      state % CFLNumberGlobal % scalar = CFLNumberGlobal
      deallocate(areaEdge)

   end subroutine computeGlobalDiagnostics

   integer function getFreeUnit()
      implicit none

      integer :: index
      logical :: isOpened

      getFreeUnit = 0
      do index = 1,99
         if((index /= 5) .and. (index /= 6)) then
            inquire(unit = index, opened = isOpened)
            if( .not. isOpened) then
               getFreeUnit = index
               return
            end if
         end if
      end do
   end function getFreeUnit

   subroutine computeFieldLocalStats(dminfo, nVertLevels, nElements, field, localSum, localMin, localMax, localVertSumMin, &
      localVertSumMax)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: localSum, localMin, localMax, localVertSumMin, &
      localVertSumMax

      localSum = sum(field)
      localMin = minval(field)
      localMax = maxval(field)
      localVertSumMin = minval(sum(field,1))
      localVertSumMax = maxval(sum(field,1))

   end subroutine computeFieldLocalStats

   subroutine computeFieldAreaWeightedLocalStats(dminfo, nVertLevels, nElements, areas, field, localSum, localMin, &
      localMax, localVertSumMin, localVertSumMax)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nElements), intent(in) :: areas
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: localSum, localMin, localMax, localVertSumMin, &
      localVertSumMax

      integer :: elementIndex

      localSum = 0.0
      do elementIndex = 1, nElements
        localSum = localSum + areas(elementIndex) * sum(field(:,elementIndex))
      end do

      localMin = minval(field)
      localMax = maxval(field)
      localVertSumMin = minval(sum(field,1))
      localVertSumMax = maxval(sum(field,1))

   end subroutine computeFieldAreaWeightedLocalStats

   subroutine computeFieldThicknessWeightedLocalStats(dminfo, nVertLevels, nElements, h, field, &
      localSum, localMin, localMax, localVertSumMin, localVertSumMax)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: h
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: localSum, localMin, localMax, localVertSumMin, &
      localVertSumMax

      real (kind=RKIND), dimension(nVertLevels, nElements) :: hTimesField

      integer :: elementIndex

      localSum = sum(h*field)
      localMin = minval(field)
      localMax = maxval(field)
      localVertSumMin = minval(sum(h*field,1))
      localVertSumMax = maxval(sum(h*field,1))

   end subroutine computeFieldThicknessWeightedLocalStats

   subroutine computeFieldVolumeWeightedLocalStats(dminfo, nVertLevels, nElements, areas, h, field, &
      localSum, localMin, localMax, localVertSumMin, localVertSumMax)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nElements), intent(in) :: areas
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: h
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: localSum, localMin, localMax, localVertSumMin, &
      localVertSumMax

      real (kind=RKIND), dimension(nVertLevels, nElements) :: hTimesField

      integer :: elementIndex

      localSum = 0.0
      do elementIndex = 1, nElements
        localSum = localSum + areas(elementIndex) * sum(h(:,elementIndex)*field(:,elementIndex))
      end do

      localMin = minval(field)
      localMax = maxval(field)
      localVertSumMin = minval(sum(h*field,1))
      localVertSumMax = maxval(sum(h*field,1))

   end subroutine computeFieldVolumeWeightedLocalStats


   subroutine computeGlobalSum(dminfo, nVertLevels, nElements, field, globalSum)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: globalSum

      real (kind=RKIND) :: localSum

      localSum = sum(field)
      call dmpar_sum_real(dminfo, localSum, globalSum)

   end subroutine computeGlobalSum

   subroutine computeAreaWeightedGlobalSum(dminfo, nVertLevels, nElements, areas, field, globalSum)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nElements), intent(in) :: areas
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: globalSum
      
      integer :: elementIndex
      real (kind=RKIND) :: localSum

      localSum = 0.
      do elementIndex = 1, nElements
        localSum = localSum + areas(elementIndex) * sum(field(:,elementIndex))
      end do
   
      call dmpar_sum_real(dminfo, localSum, globalSum)
       
   end subroutine computeAreaWeightedGlobalSum

   subroutine computeVolumeWeightedGlobalSum(dminfo, nVertLevels, nElements, areas, h, field, globalSum)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nElements), intent(in) :: areas
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: h
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: globalSum

      real (kind=RKIND), dimension(nVertLevels, nElements) :: hTimesField

      hTimesField = h*field

      call computeAreaWeightedGlobalSum(dminfo, nVertLevels, nElements, areas, hTimesField, globalSum)

   end subroutine computeVolumeWeightedGlobalSum

   subroutine computeGlobalMin(dminfo, nVertLevels, nElements, field, globalMin)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: globalMin

      real (kind=RKIND) :: localMin

      localMin = minval(field)
      call dmpar_min_real(dminfo, localMin, globalMin)

   end subroutine computeGlobalMin

   subroutine computeGlobalMax(dminfo, nVertLevels, nElements, field, globalMax)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: globalMax

      real (kind=RKIND) :: localMax

      localMax = maxval(field)
      call dmpar_max_real(dminfo, localMax, globalMax)

   end subroutine computeGlobalMax

   subroutine computeGlobalVertSumHorizMin(dminfo, nVertLevels, nElements, field, globalMin)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: globalMin

      real (kind=RKIND) :: localMin

      localMin = minval(sum(field,1))
      call dmpar_min_real(dminfo, localMin, globalMin)

   end subroutine computeGlobalVertSumHorizMin

   subroutine computeGlobalVertSumHorizMax(dminfo, nVertLevels, nElements, field, globalMax)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: globalMax

      real (kind=RKIND) :: localMax

      localMax = maxval(sum(field,1))
      call dmpar_max_real(dminfo, localMax, globalMax)

   end subroutine computeGlobalVertSumHorizMax

   subroutine computeGlobalVertThicknessWeightedSumHorizMin(dminfo, nVertLevels, nElements, h, field, globalMin)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: h, field
      real (kind=RKIND), intent(out) :: globalMin

      real (kind=RKIND) :: localMin

      localMin = minval(sum(h*field,1))
      call dmpar_min_real(dminfo, localMin, globalMin)

   end subroutine computeGlobalVertThicknessWeightedSumHorizMin

   subroutine computeGlobalVertThicknessWeightedSumHorizMax(dminfo, nVertLevels, nElements, h, field, globalMax)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: h, field
      real (kind=RKIND), intent(out) :: globalMax

      real (kind=RKIND) :: localMax

      localMax = maxval(sum(h*field,1))
      call dmpar_max_real(dminfo, localMax, globalMax)

   end subroutine computeGlobalVertThicknessWeightedSumHorizMax

end module global_diagnostics
