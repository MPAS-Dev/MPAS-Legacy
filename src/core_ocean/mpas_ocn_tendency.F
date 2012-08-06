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

   use mpas_ocn_tracer_advection

   use ocn_thick_hadv
   use ocn_thick_vadv
   use ocn_gm

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

   type (timer_node), pointer :: diagEOSTimer
   type (timer_node), pointer :: thickHadvTimer, thickVadvTimer
   type (timer_node), pointer :: velCorTimer, velVadvTimer, velPgradTimer, velHmixTimer, velForceTimer, velExpVmixTimer
   type (timer_node), pointer :: tracerHadvTimer, tracerVadvTimer, tracerHmixTimer, tracerExpVmixTimer, tracerRestoringTimer

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
             ocn_fuperp, &
             ocn_tendency_init, &
             ocn_filter_btr_mode_u, &
             ocn_filter_btr_mode_tend_u

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: hadv2nd, hadv3rd, hadv4th
   integer :: ke_cell_flag, ke_vertex_flag
   real (kind=RKIND) ::  coef_3rd_order, fCoef


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

   subroutine ocn_tend_h(tend, s, grid)!{{{
      implicit none

      type (tend_type), intent(inout) :: tend !< Input/Output: Tendency structure
      type (state_type), intent(in) :: s !< Input: State information
      type (mesh_type), intent(in) :: grid !< Input: Grid information

      real (kind=RKIND), dimension(:,:), pointer :: h_edge, wTop, tend_h, uTransport

      integer :: err

      call mpas_timer_start("ocn_tend_h")

      uTransport  => s % uTransport % array
      wTop        => s % wTop % array
      h_edge      => s % h_edge % array

      tend_h      => tend % h % array
                  
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
      ! QC Comment (3/15/12): need to make sure that uTranport is the right
      ! transport velocity here.
      call mpas_timer_start("hadv", .false., thickHadvTimer)
      call ocn_thick_hadv_tend(grid, uTransport, h_edge, tend_h, err)
      call mpas_timer_stop("hadv", thickHadvTimer)

      !
      ! height tendency: vertical advection term -d/dz(hw)
      !
      call mpas_timer_start("vadv", .false., thickVadvTimer)
      call ocn_thick_vadv_tend(grid, wtop, tend_h, err)
      call mpas_timer_stop("vadv", thickVadvTimer)

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
      implicit none

      type (tend_type), intent(inout) :: tend !< Input/Output: Tendency structure
      type (state_type), intent(in) :: s !< Input: State information
      type (diagnostics_type), intent(in) :: d !< Input: Diagnostic information
      type (mesh_type), intent(in) :: grid !< Input: Grid information

      real (kind=RKIND), dimension(:,:), pointer :: &
        h_edge, h, u, rho, zMid, pressure, &
        tend_u, circulation, vorticity, ke, ke_edge, Vor_edge, &
        MontPot, wTop, divergence, vertViscTopOfEdge

      real (kind=RKIND), dimension(:,:), pointer :: u_src

      integer :: err

      call mpas_timer_start("ocn_tend_u")

      u           => s % u % array
      rho         => s % rho % array
      wTop        => s % wTop % array
      zMid        => s % zMid % array
      h_edge      => s % h_edge % array
      vorticity   => s % vorticity % array
      divergence  => s % divergence % array
      ke          => s % ke % array
      ke_edge     => s % ke_edge % array
      Vor_edge     => s % Vor_edge % array
      MontPot     => s % MontPot % array
      pressure    => s % pressure % array
      vertViscTopOfEdge => d % vertViscTopOfEdge % array

      tend_u      => tend % u % array
                  
      u_src => grid % u_src % array

      !
      ! velocity tendency: start accumulating tendency terms
      !
      ! mrp 110516 efficiency: could remove next line and have first tend_u operation not be additive
      tend_u(:,:) = 0.0

      !
      ! velocity tendency: nonlinear Coriolis term and grad of kinetic energy
      !

      call mpas_timer_start("coriolis", .false., velCorTimer)
      call ocn_vel_coriolis_tend(grid, Vor_edge, h_edge, u, ke, tend_u, err)
      call mpas_timer_stop("coriolis", velCorTimer)

      !
      ! velocity tendency: vertical advection term -w du/dz
      !
      call mpas_timer_start("vadv", .false., velVadvTimer)
      call ocn_vel_vadv_tend(grid, u, h_edge, wtop, tend_u, err)
      call mpas_timer_stop("vadv", velVadvTimer)

      !
      ! velocity tendency: pressure gradient
      !
      call mpas_timer_start("pressure grad", .false., velPgradTimer)
      if (config_pressure_type.eq.'MontgomeryPotential') then
          call ocn_vel_pressure_grad_tend(grid, MontPot,  zMid, rho, tend_u, err)
      else
          call ocn_vel_pressure_grad_tend(grid, pressure, zMid, rho, tend_u, err)
      end if
      call mpas_timer_stop("pressure grad", velPgradTimer)

      !
      ! velocity tendency: del2 dissipation, \nu_2 \nabla^2 u
      !   computed as \nu( \nabla divergence + k \times \nabla vorticity )
      !   strictly only valid for config_h_mom_eddy_visc2 == constant
      !
      call mpas_timer_start("hmix", .false., velHmixTimer)
      call ocn_vel_hmix_tend(grid, divergence, vorticity, tend_u, err)
      call mpas_timer_stop("hmix", velHmixTimer)

      !
      ! velocity tendency: forcing and bottom drag
      !
      ! mrp 101115 note: in order to include flux boundary conditions, we will need to 
      ! know the bottom edge with nonzero velocity and place the drag there.

      call mpas_timer_start("forcings", .false., velForceTimer)
      call ocn_vel_forcing_tend(grid, u, u_src, ke_edge, h_edge, tend_u, err)
      call mpas_timer_stop("forcings", velForceTimer)

      !
      ! velocity tendency: vertical mixing d/dz( nu_v du/dz))
      !
      if (.not.config_implicit_vertical_mix) then
          call mpas_timer_start("explicit vmix", .false., velExpVmixTimer)
          call ocn_vel_vmix_tend_explicit(grid, u, h_edge, vertvisctopofedge, tend_u, err)
          call mpas_timer_stop("explicit vmix", velExpVmixTimer)
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
   subroutine ocn_tend_scalar(tend, s, d, grid, dt)!{{{
      implicit none

      type (tend_type), intent(inout) :: tend !< Input/Output: Tendency structure
      type (state_type), intent(in) :: s !< Input: State information
      type (diagnostics_type), intent(in) :: d !< Input: Diagnostic information
      type (mesh_type), intent(in) :: grid !< Input: Grid information
      real (kind=RKIND), intent(in) :: dt !< Input: Time step

      real (kind=RKIND), dimension(:,:), pointer :: &
        uTransport, h,wTop, h_edge, vertDiffTopOfCell, tend_h, uh
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        tracers, tend_tr

      integer :: err, iEdge, k

      call mpas_timer_start("ocn_tend_scalar")

      uTransport  => s % uTransport % array
      h           => s % h % array
      wTop        => s % wTop % array
      tracers     => s % tracers % array
      h_edge      => s % h_edge % array
      vertDiffTopOfCell => d % vertDiffTopOfCell % array

      tend_tr     => tend % tracers % array
      tend_h      => tend % h % array

      allocate(uh(grid % nVertLevels, grid % nEdges+1))
      !
      ! QC Comment (3/15/12): need to make sure that uTransport is the right
      ! transport velocity for the tracer.
      do iEdge = 1, grid % nEdges
         do k = 1, grid % nVertLevels
            uh(k, iEdge) = uTransport(k, iEdge) * h_edge(k, iEdge)
         end do
      end do

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

      ! Monotonoic Advection, or standard advection
      call mpas_timer_start("adv", .false., tracerHadvTimer)
      call mpas_ocn_tracer_advection_tend(tracers, uh, wTop, h, h, dt, grid, tend_h, tend_tr)
      call mpas_timer_stop("adv", tracerHadvTimer)

      !
      ! tracer tendency: del2 horizontal tracer diffusion, div(h \kappa_2 \nabla \phi)
      !
      call mpas_timer_start("hmix", .false., tracerHmixTimer)
      call ocn_tracer_hmix_tend(grid, h_edge, tracers, tend_tr, err)
      call mpas_timer_stop("hmix", tracerHmixTimer)

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
         call mpas_timer_start("explicit vmix", .false., tracerExpVmixTimer)

         call ocn_tracer_vmix_tend_explicit(grid, h, vertdifftopofcell, tracers, tend_tr, err)

         call mpas_timer_stop("explicit vmix", tracerExpVmixTimer)
      endif

! mrp 110516 printing
!print *, 'tend_tr 2',minval(tend_tr(3,1,1:nCells)),&
!                   maxval(tend_tr(3,1,1:nCells))
! mrp 110516 printing end

      !
      ! add restoring to T and S in top model layer
      !
      call mpas_timer_start("restoring", .false., tracerRestoringTimer)

      call ocn_restoring_tend(grid, h, s%index_temperature, s%index_salinity, tracers, tend_tr, err)

      call mpas_timer_stop("restoring", tracerRestoringTimer)

 10   format(2i8,10e20.10)
      call mpas_timer_stop("ocn_tend_scalar")

      deallocate(uh)

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
      implicit none

      real (kind=RKIND), intent(in) :: dt !< Input: Time step
      type (state_type), intent(inout) :: s !< Input/Output: State information
      type (mesh_type), intent(in) :: grid !< Input: Grid information


      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j
      integer :: boundaryMask, velMask, nCells, nEdges, nVertices, nVertLevels, vertexDegree, err

      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, &
        maxLevelCell, maxLevelEdgeTop, maxLevelEdgeBot, &
        maxLevelVertexBot
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, &
        verticesOnEdge, edgesOnEdge, edgesOnVertex,boundaryCell

      real (kind=RKIND) :: d2fdx2_cell1, d2fdx2_cell2, coef_3rd_order, r_tmp, invAreaCell1, invAreaCell2, invAreaTri1, invAreaTri2, invLength, h_vertex

      real (kind=RKIND), dimension(:), allocatable:: pTop

      real (kind=RKIND), dimension(:), pointer :: &
        h_s, fVertex, dvEdge, dcEdge, areaCell, areaTriangle, &
        referenceBottomDepth, ssh
      real (kind=RKIND), dimension(:,:), pointer :: &
        weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, pressure,&
        circulation, vorticity, ke, ke_edge, MontPot, wTop, zMid, &
        Vor_edge, Vor_vertex, Vor_cell, gradVor_n, gradVor_t, divergence, &
        rho, temperature, salinity, kev, kevc, uBolusGM, uTransport
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers, deriv_two
      real (kind=RKIND), dimension(:,:), allocatable:: div_u
      character :: c1*6

      h           => s % h % array
      u           => s % u % array
      uTransport  => s % uTransport % array
      uBolusGM    => s % uBolusGM % array
      v           => s % v % array
      h_edge      => s % h_edge % array
      circulation => s % circulation % array
      vorticity   => s % vorticity % array
      divergence  => s % divergence % array
      ke          => s % ke % array
      kev         => s % kev % array
      kevc        => s % kevc % array
      ke_edge     => s % ke_edge % array
      Vor_edge     => s % Vor_edge % array
      Vor_vertex   => s % Vor_vertex % array
      Vor_cell     => s % Vor_cell % array
      gradVor_n     => s % gradVor_n % array
      gradVor_t     => s % gradVor_t % array
      rho         => s % rho % array
      MontPot     => s % MontPot % array
      pressure    => s % pressure % array
      zMid        => s % zMid % array
      ssh         => s % ssh % array
      tracers     => s % tracers % array

      weightsOnEdge     => grid % weightsOnEdge % array
      kiteAreasOnVertex => grid % kiteAreasOnVertex % array
      cellsOnEdge       => grid % cellsOnEdge % array
      cellsOnVertex     => grid % cellsOnVertex % array
      verticesOnEdge    => grid % verticesOnEdge % array
      nEdgesOnCell      => grid % nEdgesOnCell % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array
      edgesOnVertex     => grid % edgesOnVertex % array
      dcEdge            => grid % dcEdge % array
      dvEdge            => grid % dvEdge % array
      areaCell          => grid % areaCell % array
      areaTriangle      => grid % areaTriangle % array
      h_s               => grid % h_s % array
      fVertex           => grid % fVertex % array
      referenceBottomDepth        => grid % referenceBottomDepth % array
      deriv_two         => grid % deriv_two % array
      maxLevelCell      => grid % maxLevelCell % array
      maxLevelEdgeTop   => grid % maxLevelEdgeTop % array
      maxLevelEdgeBot   => grid % maxLevelEdgeBot % array
      maxLevelVertexBot => grid % maxLevelVertexBot % array
                  
      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels
      vertexDegree = grid % vertexDegree

      boundaryCell => grid % boundaryCell % array


      !
      ! Compute height on cell edges at velocity locations
      !   Namelist options control the order of accuracy of the reconstructed h_edge value
      !
      ! mrp 101115 note: in order to include flux boundary conditions, we will need to 
      ! assign h_edge for maxLevelEdgeTop:maxLevelEdgeBot in the following section

      ! initialize h_edge to avoid divide by zero and NaN problems.
      h_edge = -1.0e34
      coef_3rd_order = config_coef_3rd_order

      do iEdge=1,nEdges*hadv2nd
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,maxLevelEdgeTop(iEdge)
            h_edge(k,iEdge) = 0.5 * (h(k,cell1) + h(k,cell2))
         end do
      end do

      do iEdge=1,nEdges*hadv3rd
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k=1,maxLevelEdgeTop(iEdge)

            d2fdx2_cell1 = 0.0
            d2fdx2_cell2 = 0.0

            boundaryMask = abs(transfer(.not.(boundaryCell(k,cell1) == 0 .and. boundaryCell(k,cell2) == 0), boundaryMask))

            d2fdx2_cell1 = deriv_two(1,1,iEdge) * h(k,cell1) * boundaryMask
            d2fdx2_cell2 = deriv_two(1,2,iEdge) * h(k,cell2) * boundaryMask

            !-- all edges of cell 1
            do i=1, nEdgesOnCell(cell1) * boundaryMask
               d2fdx2_cell1 = d2fdx2_cell1 + &
               deriv_two(i+1,1,iEdge) * h(k,grid % CellsOnCell % array (i,cell1))
            end do

            !-- all edges of cell 2
            do i=1, nEdgesOnCell(cell2) * boundaryMask
               d2fdx2_cell2 = d2fdx2_cell2 + &
               deriv_two(i+1,2,iEdge) * h(k,grid % CellsOnCell % array (i,cell2))
            end do

            velMask = 2*(abs(transfer(u(k,iEdge) <= 0, velMask))) - 1

            h_edge(k,iEdge) = 0.5*(h(k,cell1) + h(k,cell2)) - (dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                            + velMask * (dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12.

         end do   ! do k
      end do         ! do iEdge

      do iEdge=1,nEdges*hadv4th
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k=1,maxLevelEdgeTop(iEdge)

            d2fdx2_cell1 = 0.0
            d2fdx2_cell2 = 0.0

            boundaryMask = abs(transfer(.not.(boundaryCell(k,cell1) == 0 .and. boundaryCell(k,cell2) == 0), boundaryMask))

            d2fdx2_cell1 = deriv_two(1,1,iEdge) * h(k,cell1) * boundaryMask
            d2fdx2_cell2 = deriv_two(1,2,iEdge) * h(k,cell2) * boundaryMask

            !-- all edges of cell 1
            do i=1, nEdgesOnCell(cell1) * boundaryMask
               d2fdx2_cell1 = d2fdx2_cell1 + &
               deriv_two(i+1,1,iEdge) * h(k,grid % CellsOnCell % array (i,cell1))
            end do

            !-- all edges of cell 2
            do i=1, nEdgesOnCell(cell2) * boundaryMask
               d2fdx2_cell2 = d2fdx2_cell2 + &
               deriv_two(i+1,2,iEdge) * h(k,grid % CellsOnCell % array (i,cell2))
            end do

            h_edge(k,iEdge) =   &
                 0.5*(h(k,cell1) + h(k,cell2))      &
                    -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.

         end do   ! do k
      end do         ! do iEdge

      !
      ! set the velocity and height at dummy address
      !    used -1e34 so error clearly occurs if these values are used.
      !
      u(:,nEdges+1) = -1e34
      h(:,nCells+1) = -1e34
      tracers(s % index_temperature,:,nCells+1) = -1e34
      tracers(s % index_salinity,:,nCells+1) = -1e34

      circulation(:,:) = 0.0
      vorticity(:,:) = 0.0
      divergence(:,:) = 0.0
      ke(:,:) = 0.0
      v(:,:) = 0.0
      do iEdge=1,nEdges
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         invAreaTri1 = 1.0 / areaTriangle(vertex1)
         invAreaTri2 = 1.0 / areaTriangle(vertex2)

         !dwj: 02/23/12 arraCell(nCells+1) is still 0, this is a temporary fix
         invAreaCell1 = 1.0 / max(areaCell(cell1), 1.0)
         invAreaCell2 = 1.0 / max(areaCell(cell2), 1.0)

         do k=1,maxLevelEdgeBot(iEdge)
            ! Compute circulation and relative vorticity at each vertex
            r_tmp = dcEdge(iEdge) * u(k,iEdge)
            circulation(k,vertex1) = circulation(k,vertex1) - r_tmp
            circulation(k,vertex2) = circulation(k,vertex2) + r_tmp

            vorticity(k, vertex1) = vorticity(k, vertex1) - r_tmp * invAreaTri1
            vorticity(k, vertex2) = vorticity(k, vertex2) + r_tmp * invAreaTri2

            ! Compute the divergence at each cell center
            r_tmp = dvEdge(iEdge) * u(k, iEdge)
            divergence(k,cell1) = divergence(k,cell1) + r_tmp * invAreaCell1
            divergence(k,cell2) = divergence(k,cell2) - r_tmp * invAreaCell2

            ! Compute kinetic energy in each cell
            r_tmp = r_tmp * dcEdge(iEdge) * u(k,iEdge)
            ke(k,cell1) = ke(k,cell1) + 0.25 * r_tmp * invAreaCell1
            ke(k,cell2) = ke(k,cell2) + 0.25 * r_tmp * invAreaCell2
         end do

         ! Compute v (tangential) velocities
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
      ! Compute kinetic energy in each vertex
      !
      kev(:,:) = 0.0; kevc(:,:) = 0.0
      do iEdge=1,nEdges*ke_vertex_flag
         do k=1,nVertLevels
            r_tmp = dcEdge(iEdge) * dvEdge(iEdge) * u(k, iEdge)**2
            kev(k,verticesOnEdge(1,iEdge)) = kev(k,verticesOnEdge(1,iEdge)) + r_tmp
            kev(k,verticesOnEdge(2,iEdge)) = kev(k,verticesOnEdge(2,iEdge)) + r_tmp
         end do
      end do
      do iVertex = 1,nVertices*ke_vertex_flag
         do k=1,nVertLevels
           kev(k,iVertex) = kev(k,iVertex) / areaTriangle(iVertex) * 0.25
         enddo
      enddo
      do iVertex = 1, nVertices*ke_vertex_flag
       do i=1,grid % vertexDegree
         iCell = cellsOnVertex(i,iVertex)
         !dwj: 02/23/12 arraCell(nCells+1) is still 0, this is a temporary fix
         invAreaCell1 = 1.0 / max(areaCell(iCell), 1.0)
         do k=1,nVertLevels
           kevc(k,iCell) = kevc(k,iCell) + kiteAreasOnVertex(i, iVertex) * kev(k, iVertex) * invAreaCell1
         enddo
       enddo
      enddo

      !
      ! Compute kinetic energy in each cell by blending ke and kevc
      !
      do iCell=1,nCells*ke_vertex_flag
         do k=1,nVertLevels
            ke(k,iCell) = 5.0/8.0*ke(k,iCell) + 3.0/8.0*kevc(k,iCell)
         end do
      end do

      !
      ! Compute ke on cell edges at velocity locations for quadratic bottom drag. 
      !
      ! mrp 101025 efficiency note: we could get rid of ke_edge completely by 
      ! using sqrt(u(k,iEdge)**2 + v(k,iEdge)**2) in its place elsewhere.
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,maxLevelEdgeTop(iEdge)
            ke_edge(k,iEdge) = 0.5 * (ke(k,cell1) + ke(k,cell2))
         end do
      end do

      !
      ! Compute height at vertices, pv at vertices, and average pv to edge locations
      !  ( this computes Vor_vertex at all vertices bounding real cells and distance-1 ghost cells )
      !
      do iVertex = 1,nVertices
         invAreaTri1 = 1.0 / areaTriangle(iVertex)
         do k=1,maxLevelVertexBot(iVertex)
            h_vertex = 0.0
            do i=1,vertexDegree
               h_vertex = h_vertex + h(k,cellsOnVertex(i,iVertex)) * kiteAreasOnVertex(i,iVertex)
            end do
            h_vertex = h_vertex * invAreaTri1

            Vor_vertex(k,iVertex) = (fCoef*fVertex(iVertex) + vorticity(k,iVertex)) / h_vertex
         end do
      end do

      Vor_cell(:,:) = 0.0
      Vor_edge(:,:) = 0.0
      do iVertex = 1,nVertices
         do i=1,vertexDegree
            iCell = cellsOnVertex(i,iVertex)
            iEdge = edgesOnVertex(i,iVertex)

            !dwj: 02/23/12 arraCell(nCells+1) is still 0, this is a temporary fix
            invAreaCell1 = 1.0 / max(areaCell(iCell), 1.0)

            ! Compute pv at cell centers
            !    ( this computes Vor_cell for all real cells and distance-1 ghost cells )
            do k = 1,maxLevelCell(iCell)
               Vor_cell(k,iCell) = Vor_cell(k,iCell) + kiteAreasOnVertex(i, iVertex) * Vor_vertex(k, iVertex) * invAreaCell1
            enddo

            ! Compute pv at the edges
            !   ( this computes Vor_edge at all edges bounding real cells )
            do k=1,maxLevelEdgeBot(iEdge)
               Vor_edge(k,iEdge) = Vor_edge(k,iEdge) + 0.5 * Vor_vertex(k,iVertex)
            enddo
         enddo
      enddo

!     gradVor_n(:,:) = 0.0
!     gradVor_t(:,:) = 0.0
      do iEdge = 1,nEdges
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)
         vertex1 = verticesOnedge(1, iEdge)
         vertex2 = verticesOnedge(2, iEdge)

         invLength = 1.0 / dcEdge(iEdge)
         ! Compute gradient of PV in normal direction
         !   ( this computes gradVor_n for all edges bounding real cells )
         do k=1,maxLevelEdgeTop(iEdge)
            gradVor_n(k,iEdge) = (Vor_cell(k,cell2) - Vor_cell(k,cell1)) * invLength
         enddo

         invLength = 1.0 / dvEdge(iEdge)
         ! Compute gradient of PV in the tangent direction
         !   ( this computes gradVor_t at all edges bounding real cells and distance-1 ghost cells )
         do k = 1,maxLevelEdgeBot(iEdge)
           gradVor_t(k,iEdge) = (Vor_vertex(k,vertex2) - Vor_vertex(k,vertex1)) * invLength
         enddo

      enddo

      !
      ! Modify PV edge with upstream bias.
      !
      do iEdge = 1,nEdges
         do k = 1,maxLevelEdgeBot(iEdge)
           Vor_edge(k,iEdge) = Vor_edge(k,iEdge) &
             - config_apvm_scale_factor * dt* (  u(k,iEdge) * gradVor_n(k,iEdge) &
                          + v(k,iEdge) * gradVor_t(k,iEdge) )
         enddo
      enddo

      !
      ! equation of state
      !
      ! For an isopycnal model, density should remain constant.
      ! For zlevel, calculate in-situ density
      if (config_vert_grid_type.ne.'isopycnal') then
         call mpas_timer_start("equation of state", .false., diagEOSTimer)
         call ocn_equation_of_state_rho(s, grid, 0, 'relative', err)
      ! mrp 110324 In order to visualize rhoDisplaced, include the following
         call ocn_equation_of_state_rho(s, grid, 1, 'relative', err)
         call mpas_timer_stop("equation of state", diagEOSTimer)
      endif

      !
      ! Pressure
      ! This section must be after computing rho
      !
      ! dwj: 10/25/2011 - Need to explore isopycnal vs zlevel flags
      if (config_pressure_type.eq.'MontgomeryPotential') then

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

      else

        do iCell=1,nCells
           ! pressure for generalized coordinates
           ! assume atmospheric pressure at the surface is zero for now.
           pressure(1,iCell) = rho(1,iCell)*gravity &
              * 0.5*h(1,iCell)

           do k=2,maxLevelCell(iCell)
              pressure(k,iCell) = pressure(k-1,iCell)  &
                + 0.5*gravity*(  rho(k-1,iCell)*h(k-1,iCell) &
                               + rho(k  ,iCell)*h(k  ,iCell))
           end do

           ! Compute zMid, the z-coordinate of the middle of the layer.
           ! This is used for the rho g grad z momentum term.
           ! Note the negative sign, since referenceBottomDepth is positive
           ! and z-coordinates are negative below the surface.
           k = maxLevelCell(iCell)
           zMid(k:nVertLevels,iCell) = -referenceBottomDepth(k) + 0.5*h(k,iCell)

           do k=maxLevelCell(iCell)-1, 1, -1
              zMid(k,iCell) = zMid(k+1,iCell)  &
                + 0.5*(  h(k+1,iCell) &
                       + h(k  ,iCell))
           end do

        end do

      endif

      !
      ! Sea Surface Height
      !
      do iCell=1,nCells
         ! Start at the bottom where we know the depth, and go up.
         ! The bottom depth for this cell is 
         ! referenceBottomDepth(maxLevelCell(iCell)).
         ! Note the negative sign, since referenceBottomDepth is positive
         ! and z-coordinates are negative below the surface.

         ssh(iCell) = -referenceBottomDepth(maxLevelCell(iCell)) &
           + sum(h(1:maxLevelCell(iCell),iCell))

      end do

      !
      ! Apply the GM closure as a bolus velocity
      !
      if (config_h_kappa .GE. epsilon(0D0)) then
         call ocn_gm_compute_uBolus(s,grid)
      else
         ! mrp efficiency note: if uBolusGM is guaranteed to be zero, this can be removed.
         uBolusGM = 0.0
      end if

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
   subroutine ocn_wtop(s1,s2, grid)!{{{
      implicit none

      type (state_type), intent(inout) :: s1 !< Input/Output: State 1 information
      type (state_type), intent(inout) :: s2 !< Input/Output: State 2 information
      type (mesh_type), intent(in) :: grid !< Input: Grid information

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j, cov
      real (kind=RKIND) :: flux, vorticity_abs, h_vertex, workpv, rho0Inv, hSum

      integer :: nCells, nEdges, nVertices, nVertLevels, vertexDegree


      real (kind=RKIND), dimension(:), pointer :: &
        h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle, zstarWeight
      real (kind=RKIND), dimension(:,:), pointer :: uTransport,h,wTop, h_edge
      real (kind=RKIND), dimension(:,:), allocatable:: div_hu
      real (kind=RKIND), dimension(:), allocatable:: div_hu_btr, h_tend_col

      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, &
        verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex, &
        boundaryEdge, boundaryCell
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, &
        maxLevelCell, maxLevelEdgeTop, maxLevelEdgeBot, &
        maxLevelVertexBot,  maxLevelVertexTop

      h           => s1 % h % array
      h_edge      => s1 % h_edge % array
      uTransport  => s2 % uTransport % array
      wTop        => s2 % wTop % array

      areaCell          => grid % areaCell % array
      cellsOnEdge       => grid % cellsOnEdge % array
      maxLevelCell      => grid % maxLevelCell % array
      maxLevelEdgeBot   => grid % maxLevelEdgeBot % array
      dvEdge            => grid % dvEdge % array
      zstarWeight       => grid % zstarWeight % array

      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertLevels = grid % nVertLevels

      allocate(div_hu(nVertLevels,nCells+1), div_hu_btr(nCells+1), &
          h_tend_col(nVertLevels))

      !
      ! Compute div(h^{edge} u) for each cell
      ! See Ringler et al. (2010) jcp paper, eqn 19, 21, and fig. 3.
      !
      div_hu(:,:) = 0.0
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,maxLevelEdgeBot(iEdge)
            flux = uTransport(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge) 
            div_hu(k,cell1) = div_hu(k,cell1) + flux
            div_hu(k,cell2) = div_hu(k,cell2) - flux
         end do 
      end do 

      do iCell=1,nCells
         div_hu_btr(iCell) = 0.0
         do k=1,maxLevelCell(iCell)
            div_hu(k,iCell) = div_hu(k,iCell) / areaCell(iCell)
            div_hu_btr(iCell) = div_hu_btr(iCell) + div_hu(k,iCell)
         end do
      end do

      !
      ! vertical velocity through layer interface
      !
      !dwj: 10/25/2011 - Need to explore isopycnal vs zlevel flags
      if (config_vert_grid_type.eq.'isopycnal') then
        ! set vertical velocity to zero in isopycnal case
        wTop=0.0  

      else  ! zlevel or zstar type vertical grid

        do iCell=1,nCells

           hSum = 0.0
           do k=1,maxLevelCell(iCell)
              h_tend_col(k) = - zstarWeight(k)*h(k,iCell)*div_hu_btr(iCell)
              hSum = hSum + zstarWeight(k)*h(k,iCell)
           end do
           if(hSum > 0.0) then
             h_tend_col = h_tend_col / hSum
           else
           end if

           ! Vertical velocity through layer interface at top and 
           ! bottom is zero.
           wTop(1,iCell) = 0.0
           wTop(maxLevelCell(iCell)+1,iCell) = 0.0
           do k=maxLevelCell(iCell),2,-1
              wTop(k,iCell) = wTop(k+1,iCell) - div_hu(k,iCell) - h_tend_col(k)
           end do
        end do

      endif

      deallocate(div_hu, div_hu_btr, h_tend_col)

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
      implicit none

      type (state_type), intent(inout) :: s !< Input/Output: State information
      type (mesh_type), intent(in) :: grid !< Input: Grid information

! mrp 110512 I just split compute_tend into compute_tend_u and compute_tend_h.
!  Some of these variables can be removed, but at a later time.
      integer :: iEdge, cell1, cell2, eoe, i, j, k

      integer :: nEdgesSolve
      real (kind=RKIND), dimension(:), pointer :: fEdge
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge, u, uBcl
      type (dm_info) :: dminfo

      integer, dimension(:), pointer :: maxLevelEdgeTop, nEdgesOnEdge
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnEdge

      call mpas_timer_start("ocn_fuperp")

      u           => s % u % array
      uBcl        => s % uBcl % array
      weightsOnEdge     => grid % weightsOnEdge % array
      fEdge             => grid % fEdge % array
      maxLevelEdgeTop      => grid % maxLevelEdgeTop % array
      cellsOnEdge       => grid % cellsOnEdge % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array

      fEdge       => grid % fEdge % array

      nEdgesSolve = grid % nEdgesSolve

      !
      ! Put f*uBcl^{perp} in u as a work variable
      !
      do iEdge=1,nEdgesSolve
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
!
!  routine ocn_filter_btr_mode_u
!
!> \brief   filters barotropic mode out of the velocity variable.
!> \author  Mark Petersen
!> \date    23 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine filters barotropic mode out of the velocity variable.
!
!-----------------------------------------------------------------------
   subroutine ocn_filter_btr_mode_u(s, grid)!{{{
      implicit none

      type (state_type), intent(inout) :: s
      type (mesh_type), intent(in) :: grid

      integer :: iEdge, k, nEdges
      real (kind=RKIND) :: vertSum, uhSum, hSum
      real (kind=RKIND), dimension(:,:), pointer :: h_edge, u
      integer, dimension(:), pointer :: maxLevelEdgeTop

      call mpas_timer_start("ocn_filter_btr_mode_u")

      u           => s % u % array
      h_edge      => s % h_edge % array
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      nEdges      = grid % nEdges

      do iEdge=1,nEdges

        ! hSum is initialized outside the loop because on land boundaries 
        ! maxLevelEdgeTop=0, but I want to initialize hSum with a 
        ! nonzero value to avoid a NaN.
        uhSum = h_edge(1,iEdge) * u(1,iEdge)
        hSum  = h_edge(1,iEdge)

        do k=2,maxLevelEdgeTop(iEdge)
          uhSum = uhSum + h_edge(k,iEdge) * u(k,iEdge)
          hSum  =  hSum + h_edge(k,iEdge)
        enddo

        vertSum = uhSum/hSum
        do k=1,maxLevelEdgeTop(iEdge)
          u(k,iEdge) = u(k,iEdge) - vertSum
        enddo
      enddo ! iEdge

      call mpas_timer_stop("ocn_filter_btr_mode_u")

   end subroutine ocn_filter_btr_mode_u!}}}

!***********************************************************************
!
!  routine ocn_filter_btr_mode_tend_u
!
!> \brief   ocn_filters barotropic mode out of the u tendency
!> \author  Mark Petersen
!> \date    23 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine filters barotropic mode out of the u tendency.
!
!-----------------------------------------------------------------------
   subroutine ocn_filter_btr_mode_tend_u(tend, s, grid)!{{{
      implicit none

      type (tend_type), intent(inout) :: tend
      type (state_type), intent(in) :: s
      type (mesh_type), intent(in) :: grid

      integer :: iEdge, k, nEdges
      real (kind=RKIND) :: vertSum, uhSum, hSum
      real (kind=RKIND), dimension(:,:), pointer :: h_edge, tend_u

      integer, dimension(:), pointer :: maxLevelEdgeTop

      call mpas_timer_start("ocn_filter_btr_mode_tend_u")

      tend_u      => tend % u % array
      h_edge      => s % h_edge % array
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      nEdges      = grid % nEdges

      do iEdge=1,nEdges

        ! hSum is initialized outside the loop because on land boundaries 
        ! maxLevelEdgeTop=0, but I want to initialize hSum with a 
        ! nonzero value to avoid a NaN.
        uhSum = h_edge(1,iEdge) * tend_u(1,iEdge)
        hSum  = h_edge(1,iEdge)

        do k=2,maxLevelEdgeTop(iEdge)
          uhSum = uhSum + h_edge(k,iEdge) * tend_u(k,iEdge)
          hSum  =  hSum + h_edge(k,iEdge)
        enddo

        vertSum = uhSum/hSum
        do k=1,maxLevelEdgeTop(iEdge)
          tend_u(k,iEdge) = tend_u(k,iEdge) - vertSum
        enddo
      enddo ! iEdge

      call mpas_timer_stop("ocn_filter_btr_mode_tend_u")

   end subroutine ocn_filter_btr_mode_tend_u!}}}

!***********************************************************************
!
!  routine ocn_tendency_init
!
!> \brief   Initializes flags used within tendency routines.
!> \author  Doug Jacobsen
!> \date    4 November 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes flags related to quantities computed within
!>  other tendency routines.
!
!-----------------------------------------------------------------------
    subroutine ocn_tendency_init(err)!{{{
        integer, intent(out) :: err !< Output: Error flag

        err = 0

        coef_3rd_order = 0.

        if (config_thickness_adv_order == 2) then
            hadv2nd = 1
            hadv3rd = 0
            hadv4th = 0
        else if (config_thickness_adv_order == 3) then
            hadv2nd = 0
            hadv3rd = 1
            hadv4th = 0

            if(config_monotonic) then
                coef_3rd_order = 0.25
            else
                coef_3rd_order = 1.0
            endif
        else if (config_thickness_adv_order == 4) then
            hadv2nd = 0
            hadv3rd = 0
            hadv4th = 1
        end if


        if(config_include_KE_vertex) then
            ke_vertex_flag = 1
            ke_cell_flag = 0
        else
            ke_vertex_flag = 0
            ke_cell_flag = 1
        endif

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

    end subroutine ocn_tendency_init!}}}

!***********************************************************************

end module ocn_tendency

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
