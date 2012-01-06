!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tendency
!
!> \brief MPAS ocean tendency driver
!> \author Doug Jacobsen
!> \date   23 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the routines for computing
!>  various tendencies for the ocean. As well as routines
!>  for computing diagnostic variables, and other quantities
!>  such as wTop.
!
!-----------------------------------------------------------------------

module ocn_tendency

   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_timer

   use ocn_thick_hadv
   use ocn_thick_vadv

   use ocn_vel_coriolis
   use ocn_vel_pressure_grad
   use ocn_vel_vadv
   use ocn_vel_hmix
   use ocn_vel_forcing

   use ocn_tracer_hadv
   use ocn_tracer_vadv
   use ocn_tracer_hmix
   use ocn_restoring

   use ocn_equation_of_state
   use ocn_vmix

   use ocn_time_average

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

   public :: ocn_tend_h, &
             ocn_tend_u, &
             ocn_tend_scalar, &
             ocn_diagnostic_solve, &
             ocn_wtop, &
             ocn_fuperp

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tend_h
!
!> \brief   Computes thickness tendency
!> \author  Doug Jacobsen
!> \date    23 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the thickness tendency for the ocean
!
!-----------------------------------------------------------------------

   subroutine ocn_tend_h(tend, s, d, grid)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute height and normal wind tendencies, as well as diagnostic variables
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed tendencies for prognostic variables
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (tend_type), intent(inout) :: tend
      type (state_type), intent(in) :: s
      type (diagnostics_type), intent(in) :: d
      type (mesh_type), intent(in) :: grid

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, &
        vertex1, vertex2, eoe, i, j, err

! mrp 110512 I just split compute_tend into compute_tend_u and ocn_tend_h.
!  Most of these variables can be removed, but at a later time.
      integer :: nCells, nEdges, nVertices, nVertLevels, nEdgesSolve
      real (kind=RKIND) :: flux, vorticity_abs, h_vertex, workpv, q, &
        upstream_bias, wTopEdge, rho0Inv, r
      real (kind=RKIND), dimension(:), pointer :: &
        h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle, &
        zMidZLevel, zTopZLevel 
      real (kind=RKIND), dimension(:,:), pointer :: &
        weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, pressure, &
        tend_h, circulation, vorticity, ke, ke_edge, pv_edge, &
        MontPot, wTop, divergence, vertViscTopOfEdge
      type (dm_info) :: dminfo

      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, &
        maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot
      integer, dimension(:,:), pointer :: &
        cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, &
        edgesOnEdge, edgesOnVertex
      real (kind=RKIND) :: u_diffusion
      real (kind=RKIND), dimension(:), allocatable:: fluxVertTop,w_dudzTopEdge

      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_divergence
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_u
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_circulation, delsq_vorticity

      call mpas_timer_start("ocn_tend_h")

      h           => s % h % array
      u           => s % u % array
      v           => s % v % array
      wTop        => s % wTop % array
      h_edge      => s % h_edge % array
      circulation => s % circulation % array
      vorticity   => s % vorticity % array
      divergence  => s % divergence % array
      ke          => s % ke % array
      ke_edge     => s % ke_edge % array
      pv_edge     => s % pv_edge % array
      MontPot     => s % MontPot % array
      pressure    => s % pressure % array
      vertViscTopOfEdge => d % vertViscTopOfEdge % array

      weightsOnEdge     => grid % weightsOnEdge % array
      kiteAreasOnVertex => grid % kiteAreasOnVertex % array
      cellsOnEdge       => grid % cellsOnEdge % array
      cellsOnVertex     => grid % cellsOnVertex % array
      verticesOnEdge    => grid % verticesOnEdge % array
      nEdgesOnCell      => grid % nEdgesOnCell % array
      edgesOnCell       => grid % edgesOnCell % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array
      edgesOnVertex     => grid % edgesOnVertex % array
      dcEdge            => grid % dcEdge % array
      dvEdge            => grid % dvEdge % array
      areaCell          => grid % areaCell % array
      areaTriangle      => grid % areaTriangle % array
      h_s               => grid % h_s % array
      fVertex           => grid % fVertex % array
      fEdge             => grid % fEdge % array
      zMidZLevel        => grid % zMidZLevel % array
      zTopZLevel        => grid % zTopZLevel % array
      maxLevelCell      => grid % maxLevelCell % array
      maxLevelEdgeTop      => grid % maxLevelEdgeTop % array
      maxLevelVertexBot    => grid % maxLevelVertexBot % array

      tend_h      => tend % h % array
                  
      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nEdgesSolve = grid % nEdgesSolve
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels

      !
      ! height tendency: start accumulating tendency terms
      !
      tend_h = 0.0

      !
      ! height tendency: horizontal advection term -\nabla\cdot ( hu)
      !
      ! See Ringler et al. (2010) jcp paper, eqn 19, 21, and fig. 3. 
      ! for explanation of divergence operator.
      !
      ! for z-level, only compute height tendency for top layer.

      call mpas_timer_start("ocn_tend_h-horiz adv")

      call ocn_thick_hadv_tend(grid, u, h_edge, tend_h, err)

      call mpas_timer_stop("ocn_tend_h-horiz adv")

      !
      ! height tendency: vertical advection term -d/dz(hw)
      !
      ! Vertical advection computed for top layer of a z grid only.
      call mpas_timer_start("ocn_tend_h-vert adv")

      call ocn_thick_vadv_tend(grid, wtop, tend_h, err)

      call mpas_timer_stop("ocn_tend_h-vert adv")
      call mpas_timer_stop("ocn_tend_h")
   
   end subroutine ocn_tend_h!}}}

!***********************************************************************
!
!  routine ocn_tend_u
!
!> \brief   Computes velocity tendency
!> \author  Doug Jacobsen
!> \date    23 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the velocity tendency for the ocean
!
!-----------------------------------------------------------------------

   subroutine ocn_tend_u(tend, s, d, grid)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute height and normal wind tendencies, as well as diagnostic variables
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed tendencies for prognostic variables
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (tend_type), intent(inout) :: tend
      type (state_type), intent(in) :: s
      type (diagnostics_type), intent(in) :: d
      type (mesh_type), intent(in) :: grid

! mrp 110512 I just split compute_tend into ocn_tend_u and compute_tend_h.
!  Some of these variables can be removed, but at a later time.
      integer :: iEdge, iCell, iVertex, k, cell1, cell2, &
        vertex1, vertex2, eoe, i, j

      integer :: nCells, nEdges, nVertices, nVertLevels, nEdgesSolve, err
      real (kind=RKIND) :: flux, vorticity_abs, h_vertex, workpv, q, &
        upstream_bias, wTopEdge, rho0Inv, r, visc_vorticity_coef
      real (kind=RKIND), dimension(:), pointer :: &
        h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle, &
        zMidZLevel, zTopZLevel, meshScalingDel2, meshScalingDel4
      real (kind=RKIND), dimension(:,:), pointer :: &
        weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, pressure, &
        tend_u, circulation, vorticity, ke, ke_edge, pv_edge, &
        MontPot, wTop, divergence, vertViscTopOfEdge
      type (dm_info) :: dminfo

      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, &
        maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot
      integer, dimension(:,:), pointer :: &
        cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, &
        edgesOnEdge, edgesOnVertex
      real (kind=RKIND) :: u_diffusion
      real (kind=RKIND), dimension(:), allocatable:: fluxVertTop,w_dudzTopEdge

      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_divergence
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_u
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_circulation, delsq_vorticity


      real (kind=RKIND), dimension(:,:), pointer :: u_src
      real (kind=RKIND), parameter :: rho_ref = 1000.0

      call mpas_timer_start("ocn_tend_u")

      h           => s % h % array
      u           => s % u % array
      v           => s % v % array
      wTop        => s % wTop % array
      h_edge      => s % h_edge % array
      circulation => s % circulation % array
      vorticity   => s % vorticity % array
      divergence  => s % divergence % array
      ke          => s % ke % array
      ke_edge     => s % ke_edge % array
      pv_edge     => s % pv_edge % array
      MontPot     => s % MontPot % array
      pressure    => s % pressure % array
      vertViscTopOfEdge => d % vertViscTopOfEdge % array

      weightsOnEdge     => grid % weightsOnEdge % array
      kiteAreasOnVertex => grid % kiteAreasOnVertex % array
      cellsOnEdge       => grid % cellsOnEdge % array
      cellsOnVertex     => grid % cellsOnVertex % array
      verticesOnEdge    => grid % verticesOnEdge % array
      nEdgesOnCell      => grid % nEdgesOnCell % array
      edgesOnCell       => grid % edgesOnCell % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array
      edgesOnVertex     => grid % edgesOnVertex % array
      dcEdge            => grid % dcEdge % array
      dvEdge            => grid % dvEdge % array
      areaCell          => grid % areaCell % array
      areaTriangle      => grid % areaTriangle % array
      h_s               => grid % h_s % array
! mrp 110516 cleanup fvertex fedge not used in this subroutine
      fVertex           => grid % fVertex % array
      fEdge             => grid % fEdge % array
      zMidZLevel        => grid % zMidZLevel % array
      zTopZLevel        => grid % zTopZLevel % array
      maxLevelCell      => grid % maxLevelCell % array
      maxLevelEdgeTop      => grid % maxLevelEdgeTop % array
      maxLevelVertexBot    => grid % maxLevelVertexBot % array

      tend_u      => tend % u % array
                  
      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nEdgesSolve = grid % nEdgesSolve
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels

      u_src => grid % u_src % array

      meshScalingDel2 => grid % meshScalingDel2 % array
      meshScalingDel4 => grid % meshScalingDel4 % array

      !
      ! velocity tendency: start accumulating tendency terms
      !
      ! mrp 110516 efficiency: could remove next line and have first tend_u operation not be additive
      tend_u(:,:) = 0.0

      !
      ! velocity tendency: nonlinear Coriolis term and grad of kinetic energy
      !

      call mpas_timer_start("ocn_tend_u-coriolis")

      call ocn_vel_coriolis_tend(grid, pv_edge, h_edge, u, ke, tend_u, err)

      call mpas_timer_stop("ocn_tend_u-coriolis")

      !
      ! velocity tendency: vertical advection term -w du/dz
      !
      call mpas_timer_start("ocn_tend_u-vert adv")

      call ocn_vel_vadv_tend(grid, u, wtop, tend_u, err)

      call mpas_timer_stop("ocn_tend_u-vert adv")

      !
      ! velocity tendency: pressure gradient
      !
      call mpas_timer_start("ocn_tend_u-pressure grad")

      if (config_vert_grid_type.eq.'isopycnal') then
          call ocn_vel_pressure_grad_tend(grid, MontPot, tend_u, err)
      elseif (config_vert_grid_type.eq.'zlevel') then
          call ocn_vel_pressure_grad_tend(grid, pressure, tend_u, err)
      end if

      call mpas_timer_stop("ocn_tend_u-pressure grad")

      !
      ! velocity tendency: del2 dissipation, \nu_2 \nabla^2 u
      !   computed as \nu( \nabla divergence + k \times \nabla vorticity )
      !   strictly only valid for config_h_mom_eddy_visc2 == constant
      !
      call mpas_timer_start("ocn_tend_u-horiz mix")

      call ocn_vel_hmix_tend(grid, divergence, vorticity, tend_u, err)

      call mpas_timer_stop("ocn_tend_u-horiz mix")

      !
      ! velocity tendency: forcing and bottom drag
      !
      ! mrp 101115 note: in order to include flux boundary conditions, we will need to 
      ! know the bottom edge with nonzero velocity and place the drag there.

      call mpas_timer_start("ocn_tend_u-forcings")

      call ocn_vel_forcing_tend(grid, u, u_src, ke_edge, h_edge, tend_u, err)

      call mpas_timer_stop("ocn_tend_u-forcings")

      !
      ! velocity tendency: vertical mixing d/dz( nu_v du/dz))
      !
      if (.not.config_implicit_vertical_mix) then
          call mpas_timer_start("ocn_tend_u-explicit vert mix")

          call ocn_vel_vmix_tend_explicit(grid, u, h_edge, vertvisctopofedge, tend_u, err)

          call mpas_timer_stop("ocn_tend_u-explicit vert mix")
      endif
      call mpas_timer_stop("ocn_tend_u")

   end subroutine ocn_tend_u!}}}

!***********************************************************************
!
!  routine ocn_tendSalar
!
!> \brief   Computes scalar tendency
!> \author  Doug Jacobsen
!> \date    23 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the scalar (tracer) tendency for the ocean
!
!-----------------------------------------------------------------------

   subroutine ocn_tend_scalar(tend, s, d, grid)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !        note: the variable s % tracers really contains the tracers, 
   !              not tracers*h
   !
   ! Output: tend - computed scalar tendencies
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (tend_type), intent(inout) :: tend
      type (state_type), intent(in) :: s
      type (diagnostics_type), intent(in) :: d
      type (mesh_type), intent(in) :: grid

      integer :: i, k, iCell, iEdge, iTracer, cell1, cell2, upwindCell,&
        nEdges, nCells, nCellsSolve, nVertLevels, num_tracers, err
      real (kind=RKIND) :: invAreaCell1, invAreaCell2, tracer_turb_flux
      real (kind=RKIND) :: flux, tracer_edge, r
      real (kind=RKIND), dimension(:), pointer :: &
        h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: &
        u,h,wTop, h_edge, vertDiffTopOfCell
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        tracers, tend_tr
      integer, dimension(:,:), pointer :: boundaryEdge
      type (dm_info) :: dminfo

      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, &
        maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot
      integer, dimension(:,:), pointer :: cellsOnEdge, boundaryCell
      real (kind=RKIND), dimension(:), pointer :: zTopZLevel,zMidZLevel, &
         hRatioZLevelK, hRatioZLevelKm1, meshScalingDel2, meshScalingDel4
      real (kind=RKIND), dimension(:), allocatable:: tracer2ndDer, tracersIn, tracersOut, posZMidZLevel, &
            posZTopZLevel
      real (kind=RKIND), dimension(:,:), allocatable:: fluxVertTop, boundaryMask
      real (kind=RKIND), dimension(:,:,:), allocatable::tr_flux, tr_div, delsq_tracer, tracerTop


      real (kind=RKIND) :: d2fdx2_cell1, d2fdx2_cell2
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND) :: coef_3rd_order, flux3Coef, cSignWTop

      integer :: index_temperature, index_salinity, rrr
      real (kind=RKIND), dimension(:), pointer :: temperatureRestore, salinityRestore

      call mpas_timer_start("ocn_tend_scalar")

      u           => s % u % array
      h           => s % h % array
      boundaryCell=> grid % boundaryCell % array
      wTop        => s % wTop % array
      tracers     => s % tracers % array
      h_edge      => s % h_edge % array
      vertDiffTopOfCell => d % vertDiffTopOfCell % array

      tend_tr     => tend % tracers % array
                  
      areaCell          => grid % areaCell % array
      cellsOnEdge       => grid % cellsOnEdge % array
      dvEdge            => grid % dvEdge % array
      dcEdge            => grid % dcEdge % array
      zTopZLevel        => grid % zTopZLevel % array
      zMidZLevel        => grid % zMidZLevel % array
      hRatioZLevelK    => grid % hRatioZLevelK % array
      hRatioZLevelKm1    => grid % hRatioZLevelKm1 % array
      boundaryEdge      => grid % boundaryEdge % array
      maxLevelCell      => grid % maxLevelCell % array
      maxLevelEdgeTop      => grid % maxLevelEdgeTop % array
      maxLevelVertexBot    => grid % maxLevelVertexBot % array

      nEdges      = grid % nEdges
      nCells      = grid % nCells
      nCellsSolve = grid % nCellsSolve
      nVertLevels = grid % nVertLevels
      num_tracers = s % num_tracers

      meshScalingDel2 => grid % meshScalingDel2 % array
      meshScalingDel4 => grid % meshScalingDel4 % array


      deriv_two   => grid % deriv_two % array

      !
      ! initialize tracer tendency (RHS of tracer equation) to zero.
      !
      tend_tr(:,:,:) = 0.0

      !
      ! tracer tendency: horizontal advection term -div( h \phi u)
      !
      ! mrp 101115 note: in order to include flux boundary conditions, we will need to 
      ! assign h_edge for maxLevelEdgeTop:maxLevelEdgeBot in the compute_solve_diagnostics
      ! and then change maxLevelEdgeTop to maxLevelEdgeBot in the following section.
      ! tracer_edge at the boundary will also need to be defined for flux boundaries.

      call mpas_timer_start("ocn_tend_scalar-horiz adv")

      call ocn_tracer_hadv_tend(grid, u, h_edge, tracers, tend_tr, err)

      call mpas_timer_stop("ocn_tend_scalar-horiz adv")


      !
      ! tracer tendency: vertical advection term -d/dz( h \phi w)
      !

      call mpas_timer_start("ocn_tend_scalar-vert adv")

      call ocn_tracer_vadv_tend(grid, wtop, tracers, tend_tr, err)

      call mpas_timer_stop("ocn_tend_scalar-vert adv")

      !
      ! tracer tendency: del2 horizontal tracer diffusion, div(h \kappa_2 \nabla \phi)
      !
      call mpas_timer_start("ocn_tend_scalar-horiz diff")

      call ocn_tracer_hmix_tend(grid, h_edge, tracers, tend_tr, err)

      call mpas_timer_stop("ocn_tend_scalar-horiz diff")

! mrp 110516 printing
!print *, 'tend_tr 1',minval(tend_tr(3,1,1:nCells)),&
!                   maxval(tend_tr(3,1,1:nCells))
!print *, 'tracer  1',minval(tracers(3,1,1:nCells)),&
!                   maxval(tracers(3,1,1:nCells))
! mrp 110516 printing end

      !
      ! tracer tendency: vertical diffusion h d/dz( \kappa_v d\phi/dz)
      !
      if (.not.config_implicit_vertical_mix) then
         call mpas_timer_start("ocn_tend_scalar-explicit vert diff")

         call ocn_tracer_vmix_tend_explicit(grid, h, vertdifftopofcell, tracers, tend_tr, err)

         call mpas_timer_stop("ocn_tend_scalar-explicit vert diff")
      endif

! mrp 110516 printing
!print *, 'tend_tr 2',minval(tend_tr(3,1,1:nCells)),&
!                   maxval(tend_tr(3,1,1:nCells))
! mrp 110516 printing end

      !
      ! add restoring to T and S in top model layer
      !
      call mpas_timer_start("ocn_tend_scalar-restoring")

      call ocn_restoring_tend(grid, h, s%index_temperature, s%index_salinity, tracers, tend_tr, err)

      call mpas_timer_stop("ocn_tend_scalar-restoring")

 10   format(2i8,10e20.10)
      call mpas_timer_stop("ocn_tend_scalar")

   end subroutine ocn_tend_scalar!}}}

!***********************************************************************
!
!  routine ocn_diagnostic_solve
!
!> \brief   Computes diagnostic variables
!> \author  Doug Jacobsen
!> \date    23 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the diagnostic variables for the ocean
!
!-----------------------------------------------------------------------

   subroutine ocn_diagnostic_solve(dt, s, grid)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields used in the tendency computations
   !
   ! Input: grid - grid metadata
   !
   ! Output: s - computed diagnostics
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), intent(in) :: dt
      type (state_type), intent(inout) :: s
      type (mesh_type), intent(in) :: grid


      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j, cov
      real (kind=RKIND) :: flux, vorticity_abs, h_vertex, workpv, rho0Inv

      integer :: nCells, nEdges, nVertices, nVertLevels, vertexDegree, fCoef, err


      real (kind=RKIND), dimension(:), pointer :: &
        h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle, &
        hZLevel
      real (kind=RKIND), dimension(:,:), pointer :: &
        weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, w, pressure,&
        circulation, vorticity, ke, ke_edge, MontPot, wTop, &
        pv_edge, pv_vertex, pv_cell, gradPVn, gradPVt, divergence, &
        rho, temperature, salinity, kev, kevc
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      real (kind=RKIND), dimension(:), allocatable:: pTop
      real (kind=RKIND), dimension(:,:), allocatable:: div_u
      character :: c1*6

      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, &
        verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex, &
        boundaryEdge, boundaryCell
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, &
        maxLevelCell, maxLevelEdgeTop, maxLevelEdgeBot, &
        maxLevelVertexBot,  maxLevelVertexTop
      real (kind=RKIND) :: d2fdx2_cell1, d2fdx2_cell2
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND) :: coef_3rd_order
      real (kind=RKIND) :: r, h1, h2

      call mpas_timer_start("ocn_diagnostic_solve")

      h           => s % h % array
      u           => s % u % array
      v           => s % v % array
      wTop        => s % wTop % array
      h_edge      => s % h_edge % array
      circulation => s % circulation % array
      vorticity   => s % vorticity % array
      divergence  => s % divergence % array
      ke          => s % ke % array
      kev         => s % kev % array
      kevc        => s % kevc % array
      ke_edge     => s % ke_edge % array
      pv_edge     => s % pv_edge % array
      pv_vertex   => s % pv_vertex % array
      pv_cell     => s % pv_cell % array
      gradPVn     => s % gradPVn % array
      gradPVt     => s % gradPVt % array
      rho         => s % rho % array
      tracers     => s % tracers % array
      MontPot     => s % MontPot % array
      pressure    => s % pressure % array

      weightsOnEdge     => grid % weightsOnEdge % array
      kiteAreasOnVertex => grid % kiteAreasOnVertex % array
      cellsOnEdge       => grid % cellsOnEdge % array
      cellsOnVertex     => grid % cellsOnVertex % array
      verticesOnEdge    => grid % verticesOnEdge % array
      nEdgesOnCell      => grid % nEdgesOnCell % array
      edgesOnCell       => grid % edgesOnCell % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array
      edgesOnVertex     => grid % edgesOnVertex % array
      dcEdge            => grid % dcEdge % array
      dvEdge            => grid % dvEdge % array
      areaCell          => grid % areaCell % array
      areaTriangle      => grid % areaTriangle % array
      h_s               => grid % h_s % array
      fVertex           => grid % fVertex % array
      fEdge             => grid % fEdge % array
      hZLevel           => grid % hZLevel % array
      deriv_two         => grid % deriv_two % array
      maxLevelCell      => grid % maxLevelCell % array
      maxLevelEdgeTop   => grid % maxLevelEdgeTop % array
      maxLevelEdgeBot   => grid % maxLevelEdgeBot % array
      maxLevelVertexBot => grid % maxLevelVertexBot % array
      maxLevelVertexTop => grid % maxLevelVertexTop % array
                  
      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels
      vertexDegree = grid % vertexDegree

      boundaryEdge => grid % boundaryEdge % array
      boundaryCell => grid % boundaryCell % array

      !
      ! Compute height on cell edges at velocity locations
      !   Namelist options control the order of accuracy of the reconstructed h_edge value
      !
      ! mrp 101115 note: in order to include flux boundary conditions, we will need to 
      ! assign h_edge for maxLevelEdgeTop:maxLevelEdgeBot in the following section

      ! mrp 110516 efficiency note: For z-level, only do this on level 1.  h_edge for all
      ! lower levels is defined by hZlevel.

      call mpas_timer_start("ocn_diagnostic_solve-hEdge")

      coef_3rd_order = 0.
      if (config_thickness_adv_order == 3) coef_3rd_order = 1.0
      if (config_thickness_adv_order == 3 .and. config_monotonic) coef_3rd_order = 0.25

      if (config_thickness_adv_order == 2) then
          call mpas_timer_start("ocn_diagnostic_solve-hEdge 2")

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1,maxLevelEdgeTop(iEdge)
               h_edge(k,iEdge) = 0.5 * (h(k,cell1) + h(k,cell2))
            end do
         end do
          call mpas_timer_stop("ocn_diagnostic_solve-hEdge 2")

      else if (config_thickness_adv_order == 3) then
          call mpas_timer_start("ocn_diagnostic_solve-hEdge 3")

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            do k=1,maxLevelEdgeTop(iEdge)

               d2fdx2_cell1 = 0.0
               d2fdx2_cell2 = 0.0

               !-- if not a boundary cell
               if(boundaryCell(k,cell1).eq.0.and.boundaryCell(k,cell2).eq.0) then

                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * h(k,cell1)
                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * h(k,cell2)

                  !-- all edges of cell 1
                  do i=1, grid % nEdgesOnCell % array (cell1)
                          d2fdx2_cell1 = d2fdx2_cell1 + &
                          deriv_two(i+1,1,iEdge) * h(k,grid % CellsOnCell % array (i,cell1))
                  end do

                  !-- all edges of cell 2
                  do i=1, grid % nEdgesOnCell % array (cell2)
                          d2fdx2_cell2 = d2fdx2_cell2 + &
                          deriv_two(i+1,2,iEdge) * h(k,grid % CellsOnCell % array (i,cell2))
                  end do

               endif

               !-- if u > 0:
               if (u(k,iEdge) > 0) then
                  h_edge(k,iEdge) =     &
                       0.5*(h(k,cell1) + h(k,cell2))      &
                       -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                       -(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12.
               !-- else u <= 0:
               else
                  h_edge(k,iEdge) =     &
                       0.5*(h(k,cell1) + h(k,cell2))      &
                       -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                       +(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12.
               end if

            end do   ! do k
         end do         ! do iEdge

          call mpas_timer_stop("ocn_diagnostic_solve-hEdge 3")
      else  if (config_thickness_adv_order == 4) then
          call mpas_timer_start("ocn_diagnostic_solve-hEdge 4")

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            do k=1,maxLevelEdgeTop(iEdge)

               d2fdx2_cell1 = 0.0
               d2fdx2_cell2 = 0.0

               !-- if not a boundary cell
               if(boundaryCell(k,cell1).eq.0.and.boundaryCell(k,cell2).eq.0) then

                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * h(k,cell1)
                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * h(k,cell2)

                  !-- all edges of cell 1
                  do i=1, grid % nEdgesOnCell % array (cell1)
                          d2fdx2_cell1 = d2fdx2_cell1 + &
                          deriv_two(i+1,1,iEdge) * h(k,grid % CellsOnCell % array (i,cell1))
                  end do

                  !-- all edges of cell 2
                  do i=1, grid % nEdgesOnCell % array (cell2)
                          d2fdx2_cell2 = d2fdx2_cell2 + &
                          deriv_two(i+1,2,iEdge) * h(k,grid % CellsOnCell % array (i,cell2))
                  end do

               endif

               h_edge(k,iEdge) =   &
                    0.5*(h(k,cell1) + h(k,cell2))      &
                       -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.

            end do   ! do k
         end do         ! do iEdge

         call mpas_timer_stop("ocn_diagnostic_solve-hEdge 4")
      endif   ! if(config_thickness_adv_order == 2)
      call mpas_timer_stop("ocn_diagnostic_solve-hEdge")

      !
      ! set the velocity and height at dummy address
      !    used -1e34 so error clearly occurs if these values are used.
      !
!mrp 110516 change to zero, change back later:
      u(:,nEdges+1) = -1e34
      h(:,nCells+1) = -1e34
      tracers(s % index_temperature,:,nCells+1) = -1e34
      tracers(s % index_salinity,:,nCells+1) = -1e34

      !
      ! Compute circulation and relative vorticity at each vertex
      !
      circulation(:,:) = 0.0
      do iEdge=1,nEdges
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)
         do k=1,maxLevelEdgeBot(iEdge)
            circulation(k,vertex1) = circulation(k,vertex1) - dcEdge(iEdge) * u(k,iEdge)
            circulation(k,vertex2) = circulation(k,vertex2) + dcEdge(iEdge) * u(k,iEdge)
         end do
      end do
      do iVertex=1,nVertices
         do k=1,maxLevelVertexBot(iVertex)
            vorticity(k,iVertex) = circulation(k,iVertex) / areaTriangle(iVertex)
         end do
      end do

      !
      ! Compute the divergence at each cell center
      !
      divergence(:,:) = 0.0
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,maxLevelEdgeBot(iEdge)
             divergence(k,cell1) = divergence(k,cell1) + u(k,iEdge)*dvEdge(iEdge)
             divergence(k,cell2) = divergence(k,cell2) - u(k,iEdge)*dvEdge(iEdge)
         enddo
      end do
      do iCell = 1,nCells
         r = 1.0 / areaCell(iCell)
         do k = 1,maxLevelCell(iCell)
            divergence(k,iCell) = divergence(k,iCell) * r
         enddo
      enddo

      !
      ! Compute kinetic energy in each cell
      !
      ke(:,:) = 0.0
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,maxLevelEdgeBot(iEdge)
              ke(k,cell1) = ke(k,cell1) + 0.25 * dcEdge(iEdge) * dvEdge(iEdge) * u(k,iEdge)**2.0
              ke(k,cell2) = ke(k,cell2) + 0.25 * dcEdge(iEdge) * dvEdge(iEdge) * u(k,iEdge)**2.0
         enddo
      end do
      do iCell = 1,nCells
         do k = 1,maxLevelCell(iCell)
            ke(k,iCell) = ke(k,iCell) / areaCell(iCell)
         enddo
      enddo

      !
      ! Compute kinetic energy in each vertex
      !
      kev(:,:) = 0.0; kevc(:,:) = 0.0
      do iEdge=1,nEdges
         do k=1,nVertLevels
            kev(k,verticesOnEdge(1,iEdge)) = kev(k,verticesOnEdge(1,iEdge)) + dcEdge(iEdge) * dvEdge(iEdge) * u(k,iEdge)**2
            kev(k,verticesOnEdge(2,iEdge)) = kev(k,verticesOnEdge(2,iEdge)) + dcEdge(iEdge) * dvEdge(iEdge) * u(k,iEdge)**2
         end do
      end do
      do iVertex = 1,nVertices
         do k=1,nVertLevels
           kev(k,iVertex) = kev(k,iVertex) / areaTriangle(iVertex) / 4.0
         enddo
      enddo
      do iVertex = 1, nVertices
       do i=1,grid % vertexDegree
         iCell = cellsOnVertex(i,iVertex)
         do k=1,nVertLevels
           kevc(k,iCell) = kevc(k,iCell) + kiteAreasOnVertex(i, iVertex) * kev(k, iVertex) / areaCell(iCell)
         enddo
       enddo
      enddo

      !
      ! Compute kinetic energy in each cell by blending ke and kevc
      !
      if(config_include_KE_vertex) then
      do iCell=1,nCells
         do k=1,nVertLevels
            ke(k,iCell) = 5.0/8.0*ke(k,iCell) + 3.0/8.0*kevc(k,iCell)
         end do
      end do
      endif

      !
      ! Compute v (tangential) velocities
      !
      v(:,:) = 0.0
      do iEdge = 1,nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            ! mrp 101115 note: in order to include flux boundary conditions,
            ! the following loop may need to change to maxLevelEdgeBot
            do k = 1,maxLevelEdgeTop(iEdge) 
               v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
            end do
         end do
      end do

      !
      ! Compute ke on cell edges at velocity locations for quadratic bottom drag. 
      !
      ! mrp 101025 efficiency note: we could get rid of ke_edge completely by 
      ! using sqrt(u(k,iEdge)**2 + v(k,iEdge)**2) in its place elsewhere.
      ke_edge = 0.0  !mrp remove 0 for efficiency
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,maxLevelEdgeTop(iEdge)
            ke_edge(k,iEdge) = 0.5 * (ke(k,cell1) + ke(k,cell2))
         end do
      end do

      !
      ! Compute height at vertices, pv at vertices, and average pv to edge locations
      !  ( this computes pv_vertex at all vertices bounding real cells and distance-1 ghost cells )
      !
      if (trim(config_time_integration) == 'RK4') then
         ! for RK4, PV is really PV = (eta+f)/h
         fCoef = 1
      elseif (trim(config_time_integration) == 'split_explicit' &
          .or.trim(config_time_integration) == 'unsplit_explicit') then
         ! for split explicit, PV is eta/h because f is added separately to the momentum forcing.
! mrp temp, new should be:
         fCoef = 0
! old, for testing:
!         fCoef = 1
      end if

      do iVertex = 1,nVertices
         do k=1,maxLevelVertexBot(iVertex)
            h_vertex = 0.0
            do i=1,vertexDegree
               h_vertex = h_vertex + h(k,cellsOnVertex(i,iVertex)) * kiteAreasOnVertex(i,iVertex)
            end do
            h_vertex = h_vertex / areaTriangle(iVertex)

            pv_vertex(k,iVertex) = (fCoef*fVertex(iVertex) + vorticity(k,iVertex)) / h_vertex
         end do
      end do

      !
      ! Compute pv at cell centers
      !    ( this computes pv_cell for all real cells and distance-1 ghost cells )
      !
      pv_cell(:,:) = 0.0
      do iVertex = 1,nVertices
         do i=1,vertexDegree
            iCell = cellsOnVertex(i,iVertex)
            do k = 1,maxLevelCell(iCell)
               pv_cell(k,iCell) = pv_cell(k,iCell)  &
                  + kiteAreasOnVertex(i, iVertex) * pv_vertex(k, iVertex) &
                    / areaCell(iCell)
            enddo
         enddo
      enddo

      !
      ! Compute pv at the edges
      !   ( this computes pv_edge at all edges bounding real cells )
      !
      pv_edge(:,:) = 0.0
      do iVertex = 1,nVertices
         do i=1,vertexDegree
            iEdge = edgesOnVertex(i,iVertex)
            do k=1,maxLevelEdgeBot(iEdge)
               pv_edge(k,iEdge) =  pv_edge(k,iEdge) + 0.5 * pv_vertex(k,iVertex)
            enddo
        end do
      end do

      !
      ! Compute gradient of PV in normal direction
      !   ( this computes gradPVn for all edges bounding real cells )
      !
      gradPVn(:,:) = 0.0
      do iEdge = 1,nEdges
         do k=1,maxLevelEdgeTop(iEdge)
            gradPVn(k,iEdge) = (  pv_cell(k,cellsOnEdge(2,iEdge)) &
                                - pv_cell(k,cellsOnEdge(1,iEdge))) &
                               / dcEdge(iEdge)
         enddo
      enddo

      !
      ! Compute gradient of PV in the tangent direction
      !   ( this computes gradPVt at all edges bounding real cells and distance-1 ghost cells )
      !
      do iEdge = 1,nEdges
         do k = 1,maxLevelEdgeBot(iEdge)
           gradPVt(k,iEdge) = (  pv_vertex(k,verticesOnEdge(2,iEdge)) &
                               - pv_vertex(k,verticesOnEdge(1,iEdge))) &
                                 /dvEdge(iEdge)
         enddo
      enddo

      !
      ! Modify PV edge with upstream bias.
      !
      do iEdge = 1,nEdges
         do k = 1,maxLevelEdgeBot(iEdge)
           pv_edge(k,iEdge) = pv_edge(k,iEdge) &
             - 0.5 * dt* (  u(k,iEdge) * gradPVn(k,iEdge) &
                          + v(k,iEdge) * gradPVt(k,iEdge) )
         enddo
      enddo

      !
      ! equation of state
      !
      ! For an isopycnal model, density should remain constant.
      ! For zlevel, calculate in-situ density
      if (config_vert_grid_type.eq.'zlevel') then
         call ocn_equation_of_state_rho(s, grid, 0, 'relative', err)
      ! mrp 110324 In order to visualize rhoDisplaced, include the following
         call ocn_equation_of_state_rho(s, grid, 1, 'relative', err)
      endif

      !
      ! Pressure
      ! This section must be after computing rho
      !
      if (config_vert_grid_type.eq.'isopycnal') then

        ! For Isopycnal model.
        ! Compute pressure at top of each layer, and then
        ! Montgomery Potential.
        allocate(pTop(nVertLevels))
        do iCell=1,nCells

           ! assume atmospheric pressure at the surface is zero for now.
           pTop(1) = 0.0
           ! For isopycnal mode, p is the Montgomery Potential.
           ! At top layer it is g*SSH, where SSH may be off by a 
           ! constant (ie, h_s can be relative to top or bottom)
           MontPot(1,iCell) = gravity &
              * (h_s(iCell) + sum(h(1:nVertLevels,iCell)))

           do k=2,nVertLevels
              pTop(k) = pTop(k-1) + rho(k-1,iCell)*gravity* h(k-1,iCell)

              ! from delta M = p delta / rho
              MontPot(k,iCell) = MontPot(k-1,iCell) &
                 + pTop(k)*(1.0/rho(k,iCell) - 1.0/rho(k-1,iCell)) 
           end do

        end do
        deallocate(pTop)

      elseif (config_vert_grid_type.eq.'zlevel') then

        ! For z-level model.
        ! Compute pressure at middle of each level.  
        ! At k=1, where p is pressure at a depth of hZLevel(1)/2, not
        ! pressure at middle of layer including SSH.

        do iCell=1,nCells
           ! compute pressure for z-level coordinates
           ! assume atmospheric pressure at the surface is zero for now.

           pressure(1,iCell) = rho(1,iCell)*gravity &
              * (h(1,iCell)-0.5*hZLevel(1)) 

           do k=2,maxLevelCell(iCell)
              pressure(k,iCell) = pressure(k-1,iCell)  &
                + 0.5*gravity*(  rho(k-1,iCell)*hZLevel(k-1) &
                               + rho(k  ,iCell)*hZLevel(k  ))
           end do

        end do

      endif

      call ocn_wtop(s,grid)

      call mpas_timer_stop("ocn_diagnostic_solve")

   end subroutine ocn_diagnostic_solve!}}}

!***********************************************************************
!
!  routine ocn_wtop
!
!> \brief   Computes vertical velocity
!> \author  Doug Jacobsen
!> \date    23 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the vertical velocity in the top layer for the ocean
!
!-----------------------------------------------------------------------

   subroutine ocn_wtop(s, grid)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields used in the tendency computations
   !
   ! Input: grid - grid metadata
   !
   ! Output: s - computed diagnostics
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (state_type), intent(inout) :: s
      type (mesh_type), intent(in) :: grid

      ! mrp 110512 could clean this out, remove pointers?
      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j, cov
      real (kind=RKIND) :: flux, vorticity_abs, h_vertex, workpv, rho0Inv

      integer :: nCells, nEdges, nVertices, nVertLevels, vertexDegree


      real (kind=RKIND), dimension(:), pointer :: &
        h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle, &
        hZLevel
      real (kind=RKIND), dimension(:,:), pointer :: u,wTop
      real (kind=RKIND), dimension(:,:), allocatable:: div_u

      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, &
        verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex, &
        boundaryEdge, boundaryCell
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, &
        maxLevelCell, maxLevelEdgeTop, maxLevelEdgeBot, &
        maxLevelVertexBot,  maxLevelVertexTop

        call mpas_timer_start("wTop")

      u           => s % u % array
      wTop        => s % wTop % array

      areaCell          => grid % areaCell % array
      cellsOnEdge       => grid % cellsOnEdge % array
      hZLevel           => grid % hZLevel % array
      maxLevelCell      => grid % maxLevelCell % array
      maxLevelEdgeBot   => grid % maxLevelEdgeBot % array
      dvEdge            => grid % dvEdge % array

      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertLevels = grid % nVertLevels

      !
      ! vertical velocity through layer interface
      !
      if (config_vert_grid_type.eq.'isopycnal') then
        ! set vertical velocity to zero in isopycnal case
        wTop=0.0  

      elseif (config_vert_grid_type.eq.'zlevel') then

        !
        ! Compute div(u) for each cell
        ! See Ringler et al. (2010) jcp paper, eqn 19, 21, and fig. 3.
        !
        allocate(div_u(nVertLevels,nCells+1))
        div_u(:,:) = 0.0
        do iEdge=1,nEdges
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           do k=2,maxLevelEdgeBot(iEdge)
              flux = u(k,iEdge) * dvEdge(iEdge) 
              div_u(k,cell1) = div_u(k,cell1) + flux
              div_u(k,cell2) = div_u(k,cell2) - flux
           end do 
        end do 

        do iCell=1,nCells
           ! Vertical velocity through layer interface at top and 
           ! bottom is zero.
           wTop(1,iCell) = 0.0
           wTop(maxLevelCell(iCell)+1,iCell) = 0.0
           do k=maxLevelCell(iCell),2,-1
              wTop(k,iCell) = wTop(k+1,iCell) &
                 - div_u(k,iCell)/areaCell(iCell)*hZLevel(k)
           end do
        end do
        deallocate(div_u)

      endif

      call mpas_timer_stop("wTop")

   end subroutine ocn_wtop!}}}

!***********************************************************************
!
!  routine ocn_fuperp
!
!> \brief   Computes f u_perp
!> \author  Doug Jacobsen
!> \date    23 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes f u_perp for the ocean
!
!-----------------------------------------------------------------------

   subroutine ocn_fuperp(s, grid)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Put f*uBcl^{perp} in u as a work variable
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed tendencies for prognostic variables
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (state_type), intent(inout) :: s
      type (mesh_type), intent(in) :: grid

! mrp 110512 I just split compute_tend into compute_tend_u and compute_tend_h.
!  Some of these variables can be removed, but at a later time.
      integer :: iEdge, iCell, iVertex, k, cell1, cell2, &
        vertex1, vertex2, eoe, i, j

      integer :: nCells, nEdges, nVertices, nVertLevels, nEdgesSolve
      real (kind=RKIND) :: flux, vorticity_abs, h_vertex, workpv, q, &
        upstream_bias, wTopEdge, rho0Inv, r
      real (kind=RKIND), dimension(:), pointer :: &
        h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle, &
        zMidZLevel, zTopZLevel 
      real (kind=RKIND), dimension(:,:), pointer :: &
        weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, uBcl, v, pressure, &
        tend_u, circulation, vorticity, ke, ke_edge, pv_edge, &
        MontPot, wTop, divergence, vertViscTopOfEdge
      type (dm_info) :: dminfo

      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, &
        maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot
      integer, dimension(:,:), pointer :: &
        cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, &
        edgesOnEdge, edgesOnVertex
      real (kind=RKIND) :: u_diffusion
      real (kind=RKIND), dimension(:), allocatable:: fluxVertTop,w_dudzTopEdge

      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_divergence
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_u
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_circulation, delsq_vorticity


      real (kind=RKIND), dimension(:,:), pointer :: u_src
      real (kind=RKIND), parameter :: rho_ref = 1000.0

      call mpas_timer_start("ocn_fuperp")

      h           => s % h % array
      u           => s % u % array
      uBcl        => s % uBcl % array
      v           => s % v % array
      wTop        => s % wTop % array
      h_edge      => s % h_edge % array
      circulation => s % circulation % array
      vorticity   => s % vorticity % array
      divergence  => s % divergence % array
      ke          => s % ke % array
      ke_edge     => s % ke_edge % array
      pv_edge     => s % pv_edge % array
      MontPot     => s % MontPot % array
      pressure    => s % pressure % array

      weightsOnEdge     => grid % weightsOnEdge % array
      kiteAreasOnVertex => grid % kiteAreasOnVertex % array
      cellsOnEdge       => grid % cellsOnEdge % array
      cellsOnVertex     => grid % cellsOnVertex % array
      verticesOnEdge    => grid % verticesOnEdge % array
      nEdgesOnCell      => grid % nEdgesOnCell % array
      edgesOnCell       => grid % edgesOnCell % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array
      edgesOnVertex     => grid % edgesOnVertex % array
      dcEdge            => grid % dcEdge % array
      dvEdge            => grid % dvEdge % array
      areaCell          => grid % areaCell % array
      areaTriangle      => grid % areaTriangle % array
      h_s               => grid % h_s % array
      fVertex           => grid % fVertex % array
      fEdge             => grid % fEdge % array
      zMidZLevel        => grid % zMidZLevel % array
      zTopZLevel        => grid % zTopZLevel % array
      maxLevelCell      => grid % maxLevelCell % array
      maxLevelEdgeTop      => grid % maxLevelEdgeTop % array
      maxLevelVertexBot    => grid % maxLevelVertexBot % array

                  
      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nEdgesSolve = grid % nEdgesSolve
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels

      !
      ! Put f*uBcl^{perp} in u as a work variable
      !
      do iEdge=1,grid % nEdgesSolve
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k=1,maxLevelEdgeTop(iEdge)

            u(k,iEdge) = 0.0
            do j = 1,nEdgesOnEdge(iEdge)
               eoe = edgesOnEdge(j,iEdge)
               u(k,iEdge) = u(k,iEdge) + weightsOnEdge(j,iEdge) * uBcl(k,eoe) * fEdge(eoe) 
            end do
         end do
      end do

      call mpas_timer_stop("ocn_fuperp")

   end subroutine ocn_fuperp!}}}

!***********************************************************************

end module ocn_tendency

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
