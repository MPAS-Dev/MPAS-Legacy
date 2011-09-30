!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_vadv_stencil4
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

module ocn_tracer_vadv_stencil4

   use grid_types
   use configure
   use timer

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

   public :: ocn_tracer_vadv_stencil4_tend, &
             ocn_tracer_vadv_stencil4_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: stencil4On


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_vadv_stencil4_tend
!
!> \brief   Computes tendency term for vertical tracer advection 4th order stencil
!> \author  Phil Jones, Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the vertical advection tendency for tracer
!>  based on current state using a 4th order stencil.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vadv_stencil4_tend(grid, wTop, tracers, tend, err)!{{{

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

      integer :: nCellsSolve, iCell, k, iTracer, num_tracers, nVertLevels
      integer :: nCells

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND) :: cSingWTop, flux3Coef
      real (kind=RKIND), dimension(:), pointer :: hRatioZLevelK, hRatioZLevelKm1
      real (kind=RKIND), dimension(:,:,:), allocatable :: tracerTop


      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not. Stencil4On) return

      call timer_start("compute_scalar_tend-vert adv stencil 4")

      nCells = grid % nCells
      nCellsSolve = grid % nCellsSolve
      num_tracers = size(tracers, 1)
      nVertLevels = grid % nVertLevels
      maxLevelCell => grid % maxLevelCell % array
      hRatioZLevelK => grid % hRatioZLevelK % array
      hRatioZLevelKm1 => grid % hRatioZLevelKm1 % array

      allocate(tracerTop(num_tracers,nVertLevels+1,nCells))

      ! Compute tracerTop using 4rd order stencil [-1 7 7 -1]

      do iCell=1,nCellsSolve 
         k=2
            do iTracer=1,num_tracers
              tracerTop(iTracer,k,iCell) = &
                   hRatioZLevelK(k)  *tracers(iTracer,k-1,iCell) &
                 + hRatioZLevelKm1(k)*tracers(iTracer,k  ,iCell)
            end do
         do k=3,maxLevelCell(iCell)-1
            do iTracer=1,num_tracers
               tracerTop(iTracer,k,iCell) = &
                  (-   tracers(iTracer,k-2,iCell) &
                   +7.*tracers(iTracer,k-1,iCell) &
                   +7.*tracers(iTracer,k  ,iCell) &
                   -   tracers(iTracer,k+1,iCell) &
                  )/12.
            end do
         end do
         k=maxLevelCell(iCell)
            do iTracer=1,num_tracers
              tracerTop(iTracer,k,iCell) = &
                   hRatioZLevelK(k)  *tracers(iTracer,k-1,iCell) &
                 + hRatioZLevelKm1(k)*tracers(iTracer,k  ,iCell)
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

      deallocate(tracerTop)
      call timer_stop("compute_scalar_tend-vert adv stencil 4")

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vadv_stencil4_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_vadv_stencil4_init
!
!> \brief   Initializes ocean tracer vertical advection quantities
!> \author  Phil Jones, Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  4th order stencil based vertical tracer advection in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vadv_stencil4_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0
      stencil4On = .false.

      if(config_vert_tracer_adv_order.eq.4) then
          stencil4On = .true.
      endif


   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vadv_stencil4_init!}}}

!***********************************************************************

end module ocn_tracer_vadv_stencil4

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
