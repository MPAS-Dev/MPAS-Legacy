!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_vadv_spline3
!
!> \brief MPAS ocean vertical tracer advection driver
!> \author Doug Jacobsen
!> \date   16 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the main driver routine for computing 
!>  vertical advection tendencies.  
!
!-----------------------------------------------------------------------

module ocn_tracer_vadv_spline3

   use mpas_grid_types
   use mpas_configure
   use mpas_timer
   use mpas_spline_interpolation

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_tracer_vadv_spline3_tend, &
             ocn_tracer_vadv_spline3_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: spline3On


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_vadv_spline3_tend
!
!> \brief   Computes tendency term for vertical tracer advection 3rd order spline
!> \author  Doug Jacobsen
!> \date    16 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the vertical advection tendency for tracer
!>  based on current state using a 3rd order spline.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vadv_spline3_tend(grid, wTop, tracers, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         wTop    !< Input: vertical tracer in top layer

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers     !< Input: tracers

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend          !< Input/Output: tracer tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, nCells, nCellsSolve, k, iTracer, num_tracers, nVertLevels

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND), dimension(:), pointer :: hRatioZLevelK, &
            hRatioZLevelKm1, zTopZLevel, zMidZLevel

      real (kind=RKIND), dimension(:), allocatable :: tracer2ndDer,  &
            tracersIn, tracersOut, posZMidZLevel, posZTopZLevel
      real (kind=RKIND), dimension(:,:,:), allocatable :: tracerTop

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.spline3On) return
      ! Compute tracerTop using linear interpolation.

      call mpas_timer_start("compute_scalar_tend-vert adv spline 3")

      nCells = grid % nCells
      nCellsSolve = grid % nCellsSolve
      nVertLevels = grid % nVertLevels
      num_tracers = size(tracers, 1)
      maxLevelCell => grid % maxLevelCell % array

      hRatioZLevelK => grid % hRatioZLevelK % array
      hRatioZLevelKm1 => grid % hRatioZLevelKm1 % array
      zMidZLevel => grid % zMidZLevel % array
      zTopZLevel => grid % zTopZLevel % array

      allocate(tracerTop(num_tracers,nVertLevels+1,nCells))

      ! Compute tracerTop using cubic spline interpolation.

      allocate(tracer2ndDer(nVertLevels))
      allocate(tracersIn(nVertLevels),tracersOut(nVertLevels), &
            posZMidZLevel(nVertLevels), posZTopZLevel(nVertLevels-1))

      ! For the ocean, zlevel coordinates are negative and decreasing, 
      ! but spline functions assume increasing, so flip to positive.

      posZMidZLevel = -zMidZLevel(1:nVertLevels)
      posZTopZLevel = -zTopZLevel(2:nVertLevels)

      do iCell=1,nCellsSolve 
         ! mrp 110201 efficiency note: push tracer loop down
         ! into spline subroutines to improve efficiency
         do iTracer=1,num_tracers

            ! Place data in arrays to avoid creating new temporary arrays for every 
            ! subroutine call.  
            tracersIn(1:maxLevelCell(iCell))=tracers(iTracer,1:maxLevelCell(iCell),iCell)

            call mpas_cubic_spline_coefficients(posZMidZLevel, &
               tracersIn, maxLevelCell(iCell), tracer2ndDer)

            call mpas_interpolate_cubic_spline( &
               posZMidZLevel, tracersIn, tracer2ndDer, maxLevelCell(iCell), &
               posZTopZLevel, tracersOut, maxLevelCell(iCell)-1 )

            tracerTop(iTracer,2:maxLevelCell(iCell),iCell) = tracersOut(1:maxLevelCell(iCell)-1)

         end do
      end do

      do iCell=1,nCellsSolve 
         do k=1,maxLevelCell(iCell)  
            do iTracer=1,num_tracers
               tend(iTracer,k,iCell) = tend(iTracer,k,iCell) &
                  - (   wTop(k  ,iCell)*tracerTop(iTracer,k  ,iCell) &
                      - wTop(k+1,iCell)*tracerTop(iTracer,k+1,iCell))
            end do
         end do
      end do

      deallocate(tracer2ndDer)
      deallocate(tracersIn,tracersOut, posZMidZLevel, posZTopZLevel)
      deallocate(tracerTop)

      call mpas_timer_stop("compute_scalar_tend-vert adv spline 3")
   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vadv_spline3_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_vadv_spline3_init
!
!> \brief   Initializes ocean tracer vertical advection quantities
!> \author  Doug Jacobsen
!> \date    16 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  3rd order spline based vertical tracer advection in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vadv_spline3_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      spline3On = .false.

      if(config_vert_tracer_adv_order.eq.3) then
        spline3On = .true.
      endif


   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vadv_spline3_init!}}}

!***********************************************************************

end module ocn_tracer_vadv_spline3

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
