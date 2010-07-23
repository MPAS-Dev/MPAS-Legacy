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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                            INSTRUCTIONS                               !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! To add a new Diagnostic as a Global Stat, follow these steps.
      ! 1. Define the array to integrate, and the variable for the value above.
      ! 2. Allocate the array with the correct dimensions.
      ! 3. Fill the array with the data to be integrated.
      !     eg. G_h = Sum(h dA)/Sum(dA), See below for array filling
      ! 4. Call Function to compute Global Stat that you want.
      ! 5. Finish computing the global stat/integral
      ! 6. Write out your global stat to the file
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (dm_info), intent(in) :: dminfo
      type (grid_state), intent(inout) :: state
      type (grid_meta), intent(in) :: grid
      integer, intent(in) :: timeIndex
      real (kind=RKIND), intent(in) :: dt

      integer :: nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve, nCellsGlobal, nEdgesGlobal, nVerticesGlobal, iTracer
      integer :: nCells

      ! Step 1
      ! 1. Define the array to integrate, and the variable for the value above.
      real (kind=RKIND) ::  areaCellGlobal, areaEdgeGlobal, areaTriangleGlobal
      real (kind=RKIND), dimension(:), pointer ::  areaCell, dcEdge, dvEdge, areaTriangle, areaEdge, h_s
      real (kind=RKIND), dimension(:,:), pointer :: h, u, v, h_edge, pv_edge, pv_vertex, pv_cell, h_vertex

      real (kind=RKIND), dimension(:), allocatable :: eres, h_ref
      real (kind=RKIND), dimension(:,:), allocatable :: pe_vertex, u2, h2, hb, pv_temp, pve_bot, h_top, h_bot, cor_energy, sshr, ke_sink, pe_sink
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers

      real (kind=RKIND) ::  global_temp, gh, gpv, gpe, ge, ger, e_cor, sshr_val, h_ref_val, ke_sink_val, pe_sink_val, ke
      real (kind=RKIND) :: volumeCellGlobal, volumeEdgeGlobal, CFLNumberGlobal
      integer :: elementIndex, variableIndex, nVariables, nSums, nMaxes, nMins
      integer :: timeLevel,k,i

      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell

      integer, parameter :: kMaxVariables = 1024 ! this must be a little more than double the number of variables to be reduced
      real (kind=RKIND), dimension(kMaxVariables) :: sums, mins, maxes, averages, verticalSumMins, verticalSumMaxes, reductions

      integer :: fileID, iCell1, iCell2, iEdge
      
      cellsOnEdge => grid % cellsOnEdge % array
      edgesOnCell => grid % edgesOnCell % array

      nVertLevels = grid % nVertLevels
      nCellsSolve = grid % nCellsSolve
      nEdgesSolve = grid % nEdgesSolve
      nVerticesSolve = grid % nVerticesSolve
      nCells = grid % nCells

      h_s => grid % h_s % array
      areaCell => grid % areaCell % array
      dcEdge => grid % dcEdge % array
      dvEdge => grid % dvEdge % array
      areaTriangle => grid % areaTriangle % array
      allocate(areaEdge(1:nEdgesSolve))
      areaEdge = dcEdge(1:nEdgesSolve)*dvEdge(1:nEdgesSolve)

      h => state % h % array
      u => state % u % array
      v => state % v % array
      tracers => state % tracers % array
      h_edge => state % h_edge % array
      h_vertex => state % h_vertex % array
      pv_edge => state % pv_edge % array
      pv_vertex => state % pv_vertex % array
      pv_cell => state % pv_cell % array

      ! Step 2
      ! 2. Allocate the array with the correct dimensions.
      allocate(h_top(nVertLevels,nCellsSolve))
      allocate(h_bot(nVertLevels,nCellsSolve))
      allocate(h2(nVertLevels,nCellsSolve))
      allocate(hb(nVertLevels,nCellsSolve))
      allocate(sshr(nVertLevels,nCellsSolve))

      allocate(h_ref(nCellsSolve))

      allocate(pv_temp(nVertLevels,nVerticesSolve))
      allocate(pe_vertex(nVertLevels,nVerticesSolve))
      allocate(pve_bot(nVertLevels,nVerticesSolve))

      allocate(u2(nVertLevels,nEdgesSolve))
      allocate(cor_energy(nVertLevels,nEdgesSolve))

      allocate(ke_sink(nVertLevels,nEdgesSolve))
      allocate(pe_sink(nVertLevels,nCells))

      allocate(eres(nCellsSolve))

      pe_sink = 0

      ! Build Arrays for Global Integrals
      ! Step 3
      ! 3. Fill the array with the data to be integrated.
      !     eg. G_h = Sum(h dA)/Sum(dA), See below for array filling
      do i = 1,nVertLevels
        ! G_h Top = Sum(h dA)
        h_top(i,:) = h(i,1:nCellsSolve)*areaCell(1:nCellsSolve)
        ! G_h Bot = Sum(dA)
        h_bot(i,:) = areaCell(1:nCellsSolve)
        pv_temp(i,:) = pv_vertex(i,1:nVerticesSolve)*h_vertex(i,1:nVerticesSolve)*areaTriangle(1:nVerticesSolve) 
        pe_vertex(i,:) = pv_vertex(i,1:nVerticesSolve)*pv_vertex(i,1:nVerticesSolve)*h_vertex(i,1:nVerticesSolve)*areaTriangle(1:nVerticesSolve)
        pve_bot(i,:) = h_vertex(i,1:nVerticesSolve)*areaTriangle(1:nVerticesSolve)
        u2(i,:) = u(i,1:nEdgesSolve)*u(i,1:nEdgesSolve)*h_edge(i,1:nEdgesSolve)*areaEdge(1:nEdgesSolve)*0.5
        h2(i,:) = gravity*h(i,1:nCellsSolve)*h(i,1:nCellsSolve)*areaCell(1:nCellsSolve)*0.5
        hb(i,:) = gravity*h(i,1:nCellsSolve)*h_s(1:nCellsSolve)*areaCell(1:nCellsSolve)
        cor_energy(i,:) = h_edge(i,1:nEdgesSolve)*u(i,1:nEdgesSolve)*pv_edge(i,1:nEdgesSolve)*areaEdge(1:nEdgesSolve)
        sshr(i,:) = areaCell(1:nCellsSolve)*(h(i,1:nCellsSolve)+h_s(1:nCellsSolve))

        do k = 1,nEdgesSolve
            iCell1 = cellsOnEdge(k,1)
            iCell2 = cellsOnEdge(k,2)
            ke_sink(i,k) = areaEdge(k)*h_edge(i,k)*u(i,k)*gravity*(h(i,iCell2)+h_s(iCell2) - h(i,iCell1)-h_s(iCell1))/dcEdge(k)
            pe_sink(i,iCell1) = pe_sink(i,iCell1) + h_edge(i,k)*u(i,k)*dvEdge(k)
            pe_sink(i,iCell2) = pe_sink(i,iCell2) - h_edge(i,k)*u(i,k)*dvEdge(k)
        end do
        pe_sink(i,:) = pe_sink(i,1:nCellsSolve)*gravity*(h(i,1:nCellsSolve)+h_s(1:nCellsSolve))
      end do

      ! Step 4
      ! 4. Call Function to compute Global Stat that you want.
      ! KE and PE Sink Terms
      call computeGlobalSum(dminfo, nVertLevels, nEdgesSolve, ke_sink, ke_sink_val)
      call computeGlobalSum(dminfo, nVertLevels, nCellsSolve, pe_sink, pe_sink_val)

      call computeGlobalSum(dminfo, 1, nCellsSolve, h_ref, h_ref_val)

      ! Global Area Average of h
      call computeGlobalSum(dminfo, nVertLevels, nCellsSolve, h_top, gh)
      call computeGlobalSum(dminfo, nVertLevels, nCellsSolve, h_bot, global_temp)

      ! Step 5
      ! 5. Finish computing the global stat/integral
      gh = gh/global_temp

      ! Mean SSH for PE Res
      call computeGlobalSum(dminfo, nVertLevels, nCellsSolve, sshr, sshr_val)

      h_ref(:) = (sshr_val/global_temp)-h_s(1:nCellsSolve)


      ! Gloabl Area-Thickness Average of pv and pe
      call computeGlobalSum(dminfo, nVertLevels, nVerticesSolve, pv_temp, gpv)
      call computeGlobalSum(dminfo, nVertLevels, nVerticesSolve, pe_vertex, gpe)
      call computeGlobalSum(dminfo, nVertLevels, nVerticesSolve, pve_bot, global_temp)

      gpv = gpv/global_temp
      gpe = gpe/global_temp

      ! Global Integral of Total Energy with subtracting reference res
      call computeGlobalSum(dminfo, nVertLevels, nEdgesSolve, u2, ke)
      call computeGlobalSum(dminfo, nVertLevels, nCellsSolve, h2, global_temp)
      ge = ke + global_temp
      call computeGlobalSum(dminfo, nVertLevels, nCellsSolve, hb, global_temp)
      ge = ge + global_temp

      eres(1:nCellsSolve) = areaCell(1:nCellsSolve)*h_ref*h_ref*gravity*0.5
      call computeGlobalSum(dminfo, nVertLevels, nCellsSolve, eres, global_temp)
      ger = global_temp
      eres(1:nCellsSolve) = areaCell(1:nCellsSolve)*h_ref*h_s(1:nCellsSolve)*gravity
      call computeGlobalSum(dminfo, nVertLevels, nCellsSolve, eres, global_temp)
      ger = ger + global_temp
      ge = ge - ger

      ! Global Integral of spurious Coriolis energy for Time to KE doubling
      call computeGlobalSum(dminfo, nVertLevels, nEdgesSolve, cor_energy, e_cor)

      ! Step 6
      ! 6. Write out your global stat to the file
      if (dminfo % my_proc_id == IO_NODE) then
         fileID = getFreeUnit()
         if (timeIndex/config_stats_interval == 1) then
             open(fileID, file='GlobalIntegrals.txt',STATUS='unknown')
         else
             open(fileID, file='GlobalIntegrals.txt',ACCESS='append')
         endif 
         write(fileID,'(100es24.16)') gh, gpv, gpe, ge, e_cor, abs(ke/e_cor), ke_sink_val+pe_sink_val,abs(ke/(ke_sink_val+pe_sink_val))
         close(fileID)
      end if

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

end module global_diagnostics
