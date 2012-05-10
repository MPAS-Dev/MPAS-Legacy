!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vmix
!
!> \brief MPAS ocean vertical mixing driver
!> \author Doug Jacobsen
!> \date   19 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module is the main driver for 
!>  vertical mixing in the ocean. 
!>
!
!-----------------------------------------------------------------------

module ocn_vmix

   use mpas_grid_types
   use mpas_configure
   use mpas_timer

   use ocn_vmix_coefs_const
   use ocn_vmix_coefs_tanh
   use ocn_vmix_coefs_rich

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

   private :: tridiagonal_solve, &
              tridiagonal_solve_mult

   public :: ocn_vmix_coefs, &
             ocn_vel_vmix_tend_explicit, &
             ocn_tracer_vmix_tend_explicit, &
             ocn_vel_vmix_tend_implicit, &
             ocn_tracer_vmix_tend_implicit, &
             ocn_vmix_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: explicitOn, implicitOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vmix_coefs
!
!> \brief   Computes coefficients for vertical mixing
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the vertical mixing coefficients for momentum
!>  and tracers based user choices of mixing parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_coefs(grid, s, d, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (state_type), intent(inout) :: &
         s             !< Input/Output: state information

      type (diagnostics_type), intent(inout) :: &
         d             !< Input/Output: diagnostic information

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

      integer :: err1, err2, err3

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing coefficients
      !
      !-----------------------------------------------------------------

      err = 0

      call ocn_vmix_coefs_const_build(grid, s, d, err1)
      call ocn_vmix_coefs_tanh_build(grid, s, d, err2)
      call ocn_vmix_coefs_rich_build(grid, s, d, err3)

      err = ior(err1, ior(err2, err3))

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs!}}}

!***********************************************************************
!
!  routine ocn_vel_vmix_tendExplict
!
!> \brief   Computes tendencies for explict momentum vertical mixing
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the tendencies for explicit vertical mixing for momentum
!>  using computed coefficients.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vmix_tend_explicit(grid, u, h_edge, vertViscTopOfEdge, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         u             !< Input: velocity

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         h_edge        !< Input: thickness at edge

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertViscTopOfEdge !< Input: vertical mixing coefficients

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: tendency information

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

      integer :: iEdge, nEdgesSolve, k, nVertLevels

      integer, dimension(:), pointer :: maxLevelEdgeTop

      real (kind=RKIND), dimension(:), allocatable :: fluxVertTop

      err = 0

      if(implicitOn) return

      nEdgessolve = grid % nEdgesSolve
      nVertLevels = grid % nVertLevels
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array

      allocate(fluxVertTop(nVertLevels+1))
      fluxVertTop(1) = 0.0
      do iEdge=1,nEdgesSolve
         do k=2,maxLevelEdgeTop(iEdge)
           fluxVertTop(k) = vertViscTopOfEdge(k,iEdge) &
              * ( u(k-1,iEdge) - u(k,iEdge) ) &
              * 2 / (h_edge(k-1,iEdge) + h_edge(k,iEdge))
         enddo
         fluxVertTop(maxLevelEdgeTop(iEdge)+1) = 0.0

         do k=1,maxLevelEdgeTop(iEdge)
           tend(k,iEdge) = tend(k,iEdge) &
             + (fluxVertTop(k) - fluxVertTop(k+1)) &
             / h_edge(k,iEdge)
         enddo

      end do
      deallocate(fluxVertTop)
   !--------------------------------------------------------------------

   end subroutine ocn_vel_vmix_tend_explicit!}}}

!***********************************************************************
!
!  routine ocn_vel_vmix_tend_implicit
!
!> \brief   Computes tendencies for implicit momentum vertical mixing
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the tendencies for implicit vertical mixing for momentum
!>  using computed coefficients.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vmix_tend_implicit(grid, dt, ke_edge, vertViscTopOfEdge, h, h_edge, u, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         ke_edge        !< Input: kinetic energy at edge

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertViscTopOfEdge !< Input: vertical mixing coefficients

      real (kind=RKIND), intent(in) :: &
         dt            !< Input: time step

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         h             !< Input: thickness at cell center

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         u             !< Input: velocity

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         h_edge        !< Input: thickness at edge

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

      integer :: iEdge, nEdges, k, cell1, cell2, nVertLevels, N

      integer, dimension(:), pointer :: maxLevelEdgeTop

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(:), allocatable :: A, B, C, uTemp

      err = 0

      if(explicitOn) return

      nEdges = grid % nEdges
      nVertLevels = grid % nVertLevels
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      cellsOnEdge => grid % cellsOnEdge % array

      allocate(A(nVertLevels),B(nVertLevels),C(nVertLevels),uTemp(nVertLevels)) 
      A(1)=0

      do iEdge=1,nEdges
        N=maxLevelEdgeTop(iEdge)
        if (N.gt.0) then

         ! Compute A(k), B(k), C(k)
         ! h_edge is computed in compute_solve_diag, and is not available yet,
         ! so recompute h_edge here.
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,N
            h_edge(k,iEdge) = 0.5 * (h(k,cell1) + h(k,cell2))
         end do

         ! A is lower diagonal term
         do k=2,N
            A(k) = -2.0*dt*vertViscTopOfEdge(k,iEdge) &
               / (h_edge(k-1,iEdge) + h_edge(k,iEdge)) &
               / h_edge(k,iEdge)
         enddo

         ! C is upper diagonal term
         do k=1,N-1
            C(k) = -2.0*dt*vertViscTopOfEdge(k+1,iEdge) &
               / (h_edge(k,iEdge) + h_edge(k+1,iEdge)) &
               / h_edge(k,iEdge)
         enddo

         ! B is diagonal term
         B(1) = 1 - C(1)
         do k=2,N-1
            B(k) = 1 - A(k) - C(k)
         enddo

         ! Apply bottom drag boundary condition on the viscous term
         B(N) = 1 - A(N) + dt*config_bottom_drag_coeff  &
            *sqrt(2.0*ke_edge(k,iEdge))/h_edge(k,iEdge)

         call tridiagonal_solve(A(2:N),B,C(1:N-1),u(:,iEdge),uTemp,N)

         u(1:N,iEdge) = uTemp(1:N)
         u(N+1:nVertLevels,iEdge) = 0.0

        end if
      end do

      deallocate(A,B,C,uTemp)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vmix_tend_implicit!}}}

!***********************************************************************
!
!  routine ocn_tracer_vmix_tendExplict
!
!> \brief   Computes tendencies for explict tracer vertical mixing
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the tendencies for explicit vertical mixing for
!>  tracers using computed coefficients.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vmix_tend_explicit(grid, h, vertDiffTopOfCell, tracers, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         h        !< Input: thickness at cell center

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertDiffTopOfCell !< Input: vertical mixing coefficients

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers             !< Input: tracers

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend          !< Input/Output: tendency information

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

      integer :: iCell, nCellsSolve, k, iTracer, num_tracers, nVertLevels

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND), dimension(:,:), allocatable :: fluxVertTop

      err = 0

      if(implicitOn) return

      nCellsSolve = grid % nCellsSolve
      nVertLevels = grid % nVertLevels
      num_tracers = size(tracers, dim=1)

      maxLevelCell => grid % maxLevelCell % array

      allocate(fluxVertTop(num_tracers,nVertLevels+1))
      fluxVertTop(:,1) = 0.0
      do iCell=1,nCellsSolve 

         do k=2,maxLevelCell(iCell)
           do iTracer=1,num_tracers
             ! compute \kappa_v d\phi/dz
             fluxVertTop(iTracer,k) = vertDiffTopOfCell(k,iCell) &
                * (   tracers(iTracer,k-1,iCell)    &
                    - tracers(iTracer,k  ,iCell) )  &
                * 2 / (h(k-1,iCell) + h(k,iCell))

           enddo
         enddo
         fluxVertTop(:,maxLevelCell(iCell)+1) = 0.0

         do k=1,maxLevelCell(iCell)
           do iTracer=1,num_tracers
             ! This is h d/dz( fluxVertTop) but h and dz cancel, so 
             ! reduces to delta( fluxVertTop)
             tend(iTracer,k,iCell) = tend(iTracer,k,iCell) &
               + fluxVertTop(iTracer,k) - fluxVertTop(iTracer,k+1)
           enddo
         enddo

      enddo ! iCell loop
      deallocate(fluxVertTop)
   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vmix_tend_explicit!}}}

!***********************************************************************
!
!  routine ocn_tracer_vmix_tend_implicit
!
!> \brief   Computes tendencies for implicit tracer vertical mixing
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the tendencies for implicit vertical mixing for
!>  tracers using computed coefficients.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vmix_tend_implicit(grid, dt, vertDiffTopOfCell, h, tracers, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertDiffTopOfCell !< Input: vertical mixing coefficients

      real (kind=RKIND), intent(in) :: &
         dt            !< Input: time step

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         h             !< Input: thickness at cell center

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tracers        !< Input: tracers

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

      integer :: iCell, nCells, k, nVertLevels, num_tracers, N

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND), dimension(:), allocatable :: A,B,C
      real (kind=RKIND), dimension(:,:), allocatable :: tracersTemp

      err = 0

      if(explicitOn) return

      nCells = grid % nCells
      nVertLevels = grid % nVertLevels
      num_tracers = size(tracers, dim=1)
      maxLevelCell => grid % maxLevelCell % array

      allocate(A(nVertLevels),B(nVertLevels),C(nVertLevels),tracersTemp(num_tracers,nVertLevels))

      do iCell=1,nCells
         ! Compute A(k), B(k), C(k) for tracers
         N = maxLevelCell(iCell)

         ! A is lower diagonal term
         A(1)=0
         do k=2,N
            A(k) = -2.0*dt*vertDiffTopOfCell(k,iCell) &
                 / (h(k-1,iCell) + h(k,iCell)) / h(k,iCell)
         enddo

         ! C is upper diagonal term
         do k=1,N-1
            C(k) = -2.0*dt*vertDiffTopOfCell(k+1,iCell) &
                 / (h(k,iCell) + h(k+1,iCell)) / h(k,iCell)
         enddo
         C(N) = 0.0

         ! B is diagonal term
         do k=1,N
            B(k) = 1 - A(k) - C(k)
         enddo

         call tridiagonal_solve_mult(A(2:N),B,C(1:N-1),tracers(:,:,iCell), &
              tracersTemp, N, nVertLevels,num_tracers)

         tracers(:,1:N,iCell) = tracersTemp(:,1:N)
         tracers(:,N+1:nVertLevels,iCell) = -1e34
      end do

      deallocate(A,B,C,tracersTemp)

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vmix_tend_implicit!}}}

!***********************************************************************
!
!  routine ocn_vmix_init
!
!> \brief   Initializes ocean vertical mixing quantities
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  vertical mixing in the ocean. This primarily determines if
!>  explicit or implicit vertical mixing is to be used.
!
!-----------------------------------------------------------------------


   subroutine ocn_vmix_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      integer :: err1, err2, err3

      err = 0

      explicitOn = .true.
      implicitOn = .false.

      if(config_implicit_vertical_mix) then
          explicitOn = .false.
          implicitOn = .true.
      end if

      call ocn_vmix_coefs_const_init(err1)
      call ocn_vmix_coefs_tanh_init(err2)
      call ocn_vmix_coefs_rich_init(err3)

      err = ior(err1, ior(err2, err3))

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_init!}}}

subroutine tridiagonal_solve(a,b,c,r,x,n)!{{{
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve the matrix equation Ax=r for x, where A is tridiagonal.
! A is an nxn matrix, with:
!   a sub-diagonal, filled from 1:n-1 (a(1) appears on row 2)
!   b diagonal, filled from 1:n
!   c sup-diagonal, filled from 1:n-1  (c(1) apears on row 1)
!
! Input: a,b,c,r,n
!
! Output: x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   implicit none

   integer,intent(in) :: n
   real (KIND=RKIND), dimension(n), intent(in) :: a,b,c,r
   real (KIND=RKIND), dimension(n), intent(out) :: x
   real (KIND=RKIND), dimension(n) :: bTemp,rTemp
   real (KIND=RKIND) :: m
   integer i

   call mpas_timer_start("tridiagonal_solve")
 
   ! Use work variables for b and r
   bTemp(1) = b(1)
   rTemp(1) = r(1)
 
   ! First pass: set the coefficients
   do i = 2,n
      m = a(i-1)/bTemp(i-1)
      bTemp(i) = b(i) - m*c(i-1)
      rTemp(i) = r(i) - m*rTemp(i-1)
   end do 
 
   x(n) = rTemp(n)/bTemp(n)
   ! Second pass: back-substition
   do i = n-1, 1, -1
      x(i) = (rTemp(i) - c(i)*x(i+1))/bTemp(i)
   end do

   call mpas_timer_stop("tridiagonal_solve")
 
end subroutine tridiagonal_solve!}}}

subroutine tridiagonal_solve_mult(a,b,c,r,x,n,nDim,nSystems)!{{{
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve the matrix equation Ax=r for x, where A is tridiagonal.
! A is an nxn matrix, with:
!   a sub-diagonal, filled from 1:n-1 (a(1) appears on row 2)
!   b diagonal, filled from 1:n
!   c sup-diagonal, filled from 1:n-1  (c(1) apears on row 1)
!
! Input: a,b,c,r,n
!
! Output: x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   implicit none

   integer,intent(in) :: n, nDim, nSystems
   real (KIND=RKIND), dimension(n), intent(in) :: a,b,c
   real (KIND=RKIND), dimension(nSystems,nDim), intent(in) :: r
   real (KIND=RKIND), dimension(nSystems,nDim), intent(out) :: x
   real (KIND=RKIND), dimension(n) :: bTemp
   real (KIND=RKIND), dimension(nSystems,n) :: rTemp
   real (KIND=RKIND) :: m
   integer i,j

   call mpas_timer_start("tridiagonal_solve_mult")
 
   ! Use work variables for b and r
   bTemp(1) = b(1)
   do j = 1,nSystems
      rTemp(j,1) = r(j,1)
   end do
 
   ! First pass: set the coefficients
   do i = 2,n
      m = a(i-1)/bTemp(i-1)
      bTemp(i) = b(i) - m*c(i-1)
      do j = 1,nSystems
         rTemp(j,i) = r(j,i) - m*rTemp(j,i-1)
      end do 
   end do 
 
   do j = 1,nSystems
      x(j,n) = rTemp(j,n)/bTemp(n)
   end do
   ! Second pass: back-substition
   do i = n-1, 1, -1
      do j = 1,nSystems
         x(j,i) = (rTemp(j,i) - c(i)*x(j,i+1))/bTemp(i)
      end do
   end do
 
   call mpas_timer_stop("tridiagonal_solve_mult")

end subroutine tridiagonal_solve_mult!}}}

!***********************************************************************

end module ocn_vmix

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
