!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_advection
!
!> \brief MPAS ocean tracer advection driver
!> \author Doug Jacobsen
!> \date   03/09/12
!> \version SVN:$Id:$
!> \details
!>  This module contains driver routine for tracer advection tendencys
!>  as well as the routines for setting up advection coefficients and 
!>  initialization of the advection routines.
!
!-----------------------------------------------------------------------

module mpas_ocn_tracer_advection

   use mpas_kind_types
   use mpas_grid_types
   use mpas_configure

   use mpas_ocn_tracer_advection_std
   use mpas_ocn_tracer_advection_mono
     
   implicit none
   private
   save

   public :: mpas_ocn_tracer_advection_init,         &
             mpas_ocn_tracer_advection_coefficients, &
             mpas_ocn_tracer_advection_tend

   logical :: monotonicOn

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_ocn_tracer_advection_coefficients
!
!> \brief MPAS ocean tracer advection coefficients
!> \author Doug Jacobsen
!> \date   03/09/12
!> \version SVN:$Id:$
!> \details
!>  This routine precomputes the advection coefficients for horizontal
!>  advection of tracers.
!
!-----------------------------------------------------------------------
   subroutine mpas_ocn_tracer_advection_coefficients( grid, err )!{{{

      implicit none
      type (mesh_type) :: grid !< Input: Grid information
      integer, intent(out) :: err

      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_2nd, adv_coefs_3rd
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, advCellsForEdge, highOrderAdvectionMask, lowOrderAdvectionMask, boundaryCell
      integer, dimension(:), pointer :: nEdgesOnCell, nAdvCellsForEdge, maxLevelCell

      integer, dimension(:), pointer :: cell_list, ordered_cell_list
      integer :: cell1, cell2, iEdge, n, i, j, j_in, iCell, k, nVertLevels
      logical :: addcell, highOrderAdvection

      deriv_two => grid % deriv_two % array
      adv_coefs => grid % adv_coefs % array
      adv_coefs_2nd => grid % adv_coefs_2nd % array
      adv_coefs_3rd => grid % adv_coefs_3rd % array
      cellsOnCell => grid % cellsOnCell % array
      cellsOnEdge => grid % cellsOnEdge % array
      advCellsForEdge => grid % advCellsForEdge % array
      boundaryCell => grid % boundaryCell % array
      highOrderAdvectionMask => grid % highOrderAdvectionMask % array
      lowOrderAdvectionMask => grid % lowOrderAdvectionMask % array
      nEdgesOnCell => grid % nEdgesOnCell % array
      maxLevelCell => grid % maxLevelCell % array
      nAdvCellsForEdge => grid % nAdvCellsForEdge % array

      nVertLevels = grid % nVertLevels

      allocate(cell_list(grid % maxEdges2 + 2))
      allocate(ordered_cell_list(grid % maxEdges2 + 2))

      err = 0

      highOrderAdvectionMask = 0
      lowOrderAdvectionMask = 0
      if(config_horiz_tracer_adv_order == 2) then
        
      end if

      do iEdge = 1, grid % nEdges
        nAdvCellsForEdge(iEdge) = 0
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        

        do k = 1, nVertLevels
          if (boundaryCell(k, cell1) == 1 .or. boundaryCell(k, cell2) == 1) then
            highOrderAdvectionMask(k, iEdge) = 0
            lowOrderAdvectionMask(k, iEdge) = 1
          else
            highOrderAdvectionMask(k, iEdge) = 1
            lowOrderAdvectionMask(k, iEdge) = 0
          end if
        end do

        !
        ! do only if this edge flux is needed to update owned cells
        !
        if (cell1 <= grid%nCells .or. cell2 <= grid%nCells) then

          cell_list(1) = cell1
          cell_list(2) = cell2
          n = 2 

        !  add cells surrounding cell 1.  n is number of cells currently in list
          do i = 1, nEdgesOnCell(cell1)
            if(cellsOnCell(i,cell1) /= cell2) then
              n = n + 1
              cell_list(n) = cellsOnCell(i,cell1)
            end if
          end do

        !  add cells surrounding cell 2 (brute force approach)
          do iCell = 1, nEdgesOnCell(cell2)
            addcell = .true.
            do i=1,n
              if(cell_list(i) == cellsOnCell(iCell,cell2)) addcell = .false.
            end do
            if(addcell) then
              n = n+1
              cell_list(n) = cellsOnCell(iCell,cell2)
            end if
          end do

        ! order the list by increasing cell number (brute force approach)

          do i=1,n
            ordered_cell_list(i) = grid % nCells + 2
            j_in = 1
            do j=1,n
              if(ordered_cell_list(i) > cell_list(j) ) then
                j_in = j
                ordered_cell_list(i) = cell_list(j)
              end if
            end do
!           ordered_cell_list(i) = cell_list(j_in)
            cell_list(j_in) = grid % nCells + 3
          end do

          nAdvCellsForEdge(iEdge) = n
          do iCell = 1, nAdvCellsForEdge(iEdge)
            advCellsForEdge(iCell,iEdge) = ordered_cell_list(iCell)
          end do

        ! we have the ordered list, now construct coefficients

          adv_coefs(:,iEdge) = 0.
          adv_coefs_2nd(:,iEdge) = 0.
          adv_coefs_3rd(:,iEdge) = 0.
        
        ! pull together third and fourth order contributions to the flux
        ! first from cell1

          j_in = 0
          do j=1, n
            if( ordered_cell_list(j) == cell1 ) j_in = j
          end do
          adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(1,1,iEdge)
          adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) + deriv_two(1,1,iEdge)

          do iCell = 1, nEdgesOnCell(cell1)
            j_in = 0
            do j=1, n
              if( ordered_cell_list(j) == cellsOnCell(iCell,cell1) ) j_in = j
            end do
            adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(iCell+1,1,iEdge)
            adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) + deriv_two(iCell+1,1,iEdge)
          end do

        ! pull together third and fourth order contributions to the flux
        ! now from cell2

          j_in = 0
          do j=1, n
            if( ordered_cell_list(j) == cell2 ) j_in = j
          enddo
          adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(1,2,iEdge)
          adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) - deriv_two(1,2,iEdge)

          do iCell = 1, nEdgesOnCell(cell2)
            j_in = 0
            do j=1, n
              if( ordered_cell_list(j) == cellsOnCell(iCell,cell2) ) j_in = j
            enddo
            adv_coefs    (j_in,iEdge) = adv_coefs    (j_in,iEdge) + deriv_two(iCell+1,2,iEdge)
            adv_coefs_3rd(j_in,iEdge) = adv_coefs_3rd(j_in,iEdge) - deriv_two(iCell+1,2,iEdge)
          end do

          do j = 1,n
            adv_coefs    (j,iEdge) = - (grid % dcEdge % array (iEdge) **2) * adv_coefs    (j,iEdge) / 12.
            adv_coefs_3rd(j,iEdge) = - (grid % dcEdge % array (iEdge) **2) * adv_coefs_3rd(j,iEdge) / 12.
          end do

        ! 2nd order centered contribution - place this in the main flux weights

          j_in = 0
          do j=1, n
            if( ordered_cell_list(j) == cell1 ) j_in = j
          enddo
          adv_coefs(j_in,iEdge) = adv_coefs(j_in,iEdge) + 0.5
          adv_coefs_2nd(j_in,iEdge) = adv_coefs_2nd(j_in,iEdge) + 0.5

          j_in = 0
          do j=1, n
            if( ordered_cell_list(j) == cell2 ) j_in = j
          enddo
          adv_coefs(j_in,iEdge) = adv_coefs(j_in,iEdge) + 0.5
          adv_coefs_2nd(j_in,iEdge) = adv_coefs_2nd(j_in,iEdge) + 0.5

        !  multiply by edge length - thus the flux is just dt*ru times the results of the vector-vector multiply

          do j=1,n
            adv_coefs    (j,iEdge) = grid % dvEdge % array(iEdge) * adv_coefs    (j,iEdge)
            adv_coefs_2nd(j,iEdge) = grid % dvEdge % array(iEdge) * adv_coefs_2nd(j,iEdge)
            adv_coefs_3rd(j,iEdge) = grid % dvEdge % array(iEdge) * adv_coefs_3rd(j,iEdge)
          end do

        end if  ! only do for edges of owned-cells
        
      end do ! end loop over edges

      deallocate(cell_list)
      deallocate(ordered_cell_list)

      ! If 2nd order advection, set masks appropriately.
      if(config_horiz_tracer_adv_order == 2) then
        lowOrderAdvectionMask = 1
        highOrderAdvectionMask = 0
      end if

      if (maxval(highOrderAdvectionMask+lowOrderAdvectionMask) > 1) then
        write(*,*) "Masks don't sum to 1."
        err = 1
      endif

   end subroutine mpas_ocn_tracer_advection_coefficients!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_ocn_tracer_advection_tend
!
!> \brief MPAS ocean tracer advection tendency
!> \author Doug Jacobsen
!> \date   03/09/12
!> \version SVN:$Id:$
!> \details
!>  This routine is the driver routine for computing the tendency for
!>  advection of tracers.
!
!-----------------------------------------------------------------------
   subroutine mpas_ocn_tracer_advection_tend(tracers, uh, w, h, verticalCellSize, dt, grid, tend_h, tend)!{{{

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tend !< Input/Output: tracer tendency
      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers !< Input/Output: tracer values
      real (kind=RKIND), dimension(:,:), intent(in) :: uh !< Input: Thickness weighted horizontal velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: w  !< Input: Vertical velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: h  !< Input: Thickness field
      real (kind=RKIND), dimension(:,:), intent(in) :: verticalCellSize !< Input: Distance between vertical interfaces of a cell
      real (kind=RKIND), intent(in) :: dt !< Input: Time step
      type (mesh_type), intent(in) :: grid !< Input: grid information
      real (kind=RKIND), dimension(:,:), intent(in) :: tend_h !< Input: Thickness tendency information

      if(monotonicOn) then
         call mpas_ocn_tracer_advection_mono_tend(tracers, uh, w, h, verticalCellSize, dt, grid, tend_h, tend)
      else
         call mpas_ocn_tracer_advection_std_tend(tracers, uh, w, verticalCellSize, grid, tend)
      endif
   end subroutine mpas_ocn_tracer_advection_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_ocn_tracer_advection_init
!
!> \brief MPAS ocean tracer advection tendency
!> \author Doug Jacobsen
!> \date   03/09/12
!> \version SVN:$Id:$
!> \details
!>  This routine is the driver routine for initialization of 
!>  the tracer advection routines.
!
!-----------------------------------------------------------------------
   subroutine mpas_ocn_tracer_advection_init(err)!{{{

      integer, intent(inout) :: err !< Input/Output: Error flag

      integer :: err_tmp

      err = 0

      call mpas_ocn_tracer_advection_std_init(err_tmp)
      call mpas_ocn_tracer_advection_mono_init(err_tmp)

      err = ior(err, err_tmp)

      monotonicOn = .false.

      if(config_monotonic) then
         monotonicOn = .true.
      endif

   end subroutine mpas_ocn_tracer_advection_init!}}}

end module mpas_ocn_tracer_advection
