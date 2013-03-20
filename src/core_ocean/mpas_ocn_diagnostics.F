!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_diagnostics
!
!> \brief MPAS ocean diagnostics driver
!> \author Mark Petersen
!> \date   23 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the routines for computing
!>  diagnostic variables, and other quantities such as wTop.
!
!-----------------------------------------------------------------------

module ocn_diagnostics

   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_timer

   use ocn_gm
   use ocn_equation_of_state

   implicit none
   private
   save

   type (timer_node), pointer :: diagEOSTimer

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

   public :: ocn_diagnostic_solve, &
             ocn_wtop, &
             ocn_fuperp, &
             ocn_filter_btr_mode_u, &
             ocn_filter_btr_mode_tend_u, &
             ocn_diagnostics_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: ke_cell_flag, ke_vertex_flag
   real (kind=RKIND) ::  coef_3rd_order, fCoef

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_diagnostic_solve
!
!> \brief   Computes diagnostic variables
!> \author  Mark Petersen
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
        verticesOnEdge, edgesOnEdge, edgesOnVertex,boundaryCell, kiteIndexOnCell, &
        verticesOnCell, edgeSignOnVertex, edgeSignOnCell, edgesOnCell

      real (kind=RKIND) :: d2fdx2_cell1, d2fdx2_cell2, coef_3rd_order, r_tmp, &
        invAreaCell1, invAreaCell2, invAreaTri1, invAreaTri2, invLength, h_vertex, coef

      real (kind=RKIND), dimension(:), allocatable:: pTop, div_hu

      real (kind=RKIND), dimension(:), pointer :: &
        bottomDepth, fVertex, dvEdge, dcEdge, areaCell, areaTriangle, ssh, seaSurfacePressure
      real (kind=RKIND), dimension(:,:), pointer :: &
        weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, pressure,&
        circulation, vorticity, ke, ke_edge, MontPot, wTop, zMid, &
        Vor_edge, Vor_vertex, Vor_cell, gradVor_n, gradVor_t, divergence, &
        rho, rhoDisplaced, temperature, salinity, kev, kevc, uBolusGM, uTransport, &
        vertVelocityTop, BruntVaisalaFreqTop
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers, deriv_two
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
      Vor_edge    => s % Vor_edge % array
      Vor_vertex  => s % Vor_vertex % array
      Vor_cell    => s % Vor_cell % array
      gradVor_n   => s % gradVor_n % array
      gradVor_t   => s % gradVor_t % array
      rho         => s % rho % array
      rhoDisplaced=> s % rhoDisplaced % array
      MontPot     => s % MontPot % array
      pressure    => s % pressure % array
      zMid        => s % zMid % array
      ssh         => s % ssh % array
      tracers     => s % tracers % array
      vertVelocityTop => s % vertVelocityTop % array
      BruntVaisalaFreqTop => s % BruntVaisalaFreqTop % array

      weightsOnEdge     => grid % weightsOnEdge % array
      kiteAreasOnVertex => grid % kiteAreasOnVertex % array
      cellsOnEdge       => grid % cellsOnEdge % array
      cellsOnVertex     => grid % cellsOnVertex % array
      verticesOnEdge    => grid % verticesOnEdge % array
      nEdgesOnCell      => grid % nEdgesOnCell % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnCell       => grid % edgesOnCell % array
      edgesOnEdge       => grid % edgesOnEdge % array
      edgesOnVertex     => grid % edgesOnVertex % array
      dcEdge            => grid % dcEdge % array
      dvEdge            => grid % dvEdge % array
      areaCell          => grid % areaCell % array
      areaTriangle      => grid % areaTriangle % array
      bottomDepth       => grid % bottomDepth % array
      fVertex           => grid % fVertex % array
      deriv_two         => grid % deriv_two % array
      maxLevelCell      => grid % maxLevelCell % array
      maxLevelEdgeTop   => grid % maxLevelEdgeTop % array
      maxLevelEdgeBot   => grid % maxLevelEdgeBot % array
      maxLevelVertexBot => grid % maxLevelVertexBot % array
      kiteIndexOnCell => grid % kiteIndexOnCell % array
      verticesOnCell => grid % verticesOnCell % array

      seaSurfacePressure => grid % seaSurfacePressure % array
                  
      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels
      vertexDegree = grid % vertexDegree

      boundaryCell => grid % boundaryCell % array

      edgeSignOnVertex => grid % edgeSignOnVertex % array
      edgeSignOnCell => grid % edgeSignOnCell % array

      !
      ! Compute height on cell edges at velocity locations
      !   Namelist options control the order of accuracy of the reconstructed h_edge value
      !

      ! initialize h_edge to avoid divide by zero and NaN problems.
      h_edge = -1.0e34
      coef_3rd_order = config_coef_3rd_order

      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,maxLevelEdgeTop(iEdge)
            h_edge(k,iEdge) = 0.5 * (h(k,cell1) + h(k,cell2))
         end do
      end do

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
      vertVelocityTop(:,:)=0.0
      ke(:,:) = 0.0
      v(:,:) = 0.0
      do iVertex = 1, nVertices
         invAreaTri1 = 1.0 / areaTriangle(iVertex)
         do i = 1, vertexDegree
            iEdge = edgesOnVertex(i, iVertex)
            do k = 1, maxLevelVertexBot(iVertex)
              r_tmp = dcEdge(iEdge) * u(k, iEdge)

              circulation(k, iVertex) = circulation(k, iVertex) + edgeSignOnVertex(i, iVertex) * r_tmp 
              vorticity(k, iVertex) = vorticity(k, iVertex) + edgeSignOnVertex(i, iVertex) * r_tmp * invAreaTri1
            end do
         end do
      end do

      allocate(div_hu(nVertLevels))
      do iCell = 1, nCells
         div_hu(:) = 0.0
         invAreaCell1 = 1.0 / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelCell(iCell)
               r_tmp = dvEdge(iEdge) * u(k, iEdge) * invAreaCell1

               divergence(k, iCell) = divergence(k, iCell) - edgeSignOnCell(i, iCell) * r_tmp
               div_hu(k)    = div_hu(k) - h_edge(k, iEdge) * edgeSignOnCell(i, iCell) * r_tmp 
               ke(k, iCell) = ke(k, iCell) + 0.25 * r_tmp * dcEdge(iEdge) * u(k,iEdge)
            end do
         end do
         ! Vertical velocity at bottom (maxLevelCell(iCell)+1) is zero, initialized above.
         do k=maxLevelCell(iCell),1,-1
            vertVelocityTop(k,iCell) = vertVelocityTop(k+1,iCell) - div_hu(k)
         end do         
      end do
      deallocate(div_hu)

      do iEdge=1,nEdges
         ! Compute v (tangential) velocities
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            do k = 1,maxLevelEdgeTop(iEdge) 
               v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
            end do
         end do
      end do

      !
      ! Compute kinetic energy in each vertex
      !
      kev(:,:) = 0.0; kevc(:,:) = 0.0
      do iVertex = 1, nVertices*ke_vertex_flag
        do i = 1, vertexDegree
          iEdge = edgesOnVertex(i, iVertex)
          r_tmp = dcEdge(iEdge) * dvEdge(iEdge) * 0.25 / areaTriangle(iVertex)
          do k = 1, nVertLevels
            kev(k, iVertex) = kev(k, iVertex) + r_tmp * u(k, iEdge)**2
          end do
        end do
      end do

      do iCell = 1, nCells*ke_vertex_flag
        invAreaCell1 = 1.0 / areaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          j = kiteIndexOnCell(i, iCell)
          iVertex = verticesOnCell(i, iCell)
          do k = 1, nVertLevels
            kevc(k, iCell) = kevc(k, iCell) + kiteAreasOnVertex(j, iVertex) * kev(k, iVertex) * invAreaCell1
          end do
        end do
      end do

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
      do iEdge = 1, nEdges
        vertex1 = verticesOnEdge(1, iEdge)
        vertex2 = verticesOnEdge(2, iEdge)
        do k = 1, maxLevelEdgeBot(iEdge)
          Vor_edge(k, iEdge) = 0.5 * (Vor_vertex(k, vertex1) + Vor_vertex(k, vertex2))
        end do
      end do

      do iCell = 1, nCells
        invAreaCell1 = 1.0 / areaCell(iCell)

        do i = 1, nEdgesOnCell(iCell)
          j = kiteIndexOnCell(i, iCell)
          iVertex = verticesOnCell(i, iCell)
          do k = 1, maxLevelCell(iCell)
            Vor_cell(k, iCell) = Vor_cell(k, iCell) + kiteAreasOnVertex(j, iVertex) * Vor_vertex(k, iVertex) * invAreaCell1
          end do
        end do
      end do

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
      if (config_vert_coord_movement.ne.'isopycnal') then
         call mpas_timer_start("equation of state", .false., diagEOSTimer)

         ! compute in-place density
         call ocn_equation_of_state_rho(s, grid, 0, 'relative', err)

         ! compute rhoDisplaced, the potential density referenced to the top layer
         call ocn_equation_of_state_rho(s, grid, 1, 'relative', err)

         call mpas_timer_stop("equation of state", diagEOSTimer)
      endif

      !
      ! Pressure
      ! This section must be after computing rho
      !
      if (config_pressure_gradient_type.eq.'MontgomeryPotential') then

        ! For Isopycnal model.
        ! Compute pressure at top of each layer, and then
        ! Montgomery Potential.
        allocate(pTop(nVertLevels))
        do iCell=1,nCells

           ! assume atmospheric pressure at the surface is zero for now.
           pTop(1) = 0.0
           ! For isopycnal mode, p is the Montgomery Potential.
           ! At top layer it is g*SSH, where SSH may be off by a 
           ! constant (ie, bottomDepth can be relative to top or bottom)
           MontPot(1,iCell) = gravity &
              * (bottomDepth(iCell) + sum(h(1:nVertLevels,iCell)))

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
           ! Pressure for generalized coordinates.
           ! Pressure at top surface may be due to atmospheric pressure
           ! or an ice-shelf depression. 
           pressure(1,iCell) = seaSurfacePressure(iCell) + rho(1,iCell)*gravity &
              * 0.5*h(1,iCell)

           do k=2,maxLevelCell(iCell)
              pressure(k,iCell) = pressure(k-1,iCell)  &
                + 0.5*gravity*(  rho(k-1,iCell)*h(k-1,iCell) &
                               + rho(k  ,iCell)*h(k  ,iCell))
           end do

           ! Compute zMid, the z-coordinate of the middle of the layer.
           ! This is used for the rho g grad z momentum term.
           ! Note the negative sign, since bottomDepth is positive
           ! and z-coordinates are negative below the surface.
           k = maxLevelCell(iCell)
           zMid(k:nVertLevels,iCell) = -bottomDepth(iCell) + 0.5*h(k,iCell)

           do k=maxLevelCell(iCell)-1, 1, -1
              zMid(k,iCell) = zMid(k+1,iCell)  &
                + 0.5*(  h(k+1,iCell) &
                       + h(k  ,iCell))
           end do

        end do

      endif

      !
      ! Brunt-Vaisala frequency
      !
      coef = -gravity/config_rho0
      do iCell=1,nCells
         BruntVaisalaFreqTop(1,iCell) = 0.0
         do k=2,maxLevelCell(iCell)
            BruntVaisalaFreqTop(k,iCell) = coef * (rhoDisplaced(k-1,iCell) - rhoDisplaced(k,iCell)) & 
              / (zMid(k-1,iCell) - zMid(k,iCell))
          end do
      end do

      !
      ! Sea Surface Height
      !
      do iCell=1,nCells
         ! Start at the bottom where we know the depth, and go up.
         ! The bottom depth for this cell is bottomDepth(iCell).
         ! Note the negative sign, since bottomDepth is positive
         ! and z-coordinates are negative below the surface.

         ssh(iCell) = - bottomDepth(iCell) + sum(h(1:maxLevelCell(iCell),iCell))

      end do

      !
      ! Apply the GM closure as a bolus velocity
      !
      if (config_h_kappa .GE. epsilon(0D0)) then
         call ocn_gm_compute_uBolus(s,grid)
      else
         uBolusGM = 0.0
      end if

   end subroutine ocn_diagnostic_solve!}}}

!***********************************************************************
!
!  routine ocn_wtop
!
!> \brief   Computes vertical transport
!> \author  Mark Petersen
!> \date    23 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the vertical transport through the top of each 
!>  cell.
!
!-----------------------------------------------------------------------
   subroutine ocn_wtop(grid,h,h_edge,u,wTop, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         h    !< Input: thickness

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         h_edge     !< Input: h interpolated to an edge

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         u     !< Input: transport

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         wTop     !< Output: vertical transport at top of cell

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j, cov
      real (kind=RKIND) :: flux, vorticity_abs, h_vertex, workpv, rho0Inv, hSum, invAreaCell

      integer :: nCells, nEdges, nVertices, nVertLevels, vertexDegree


      real (kind=RKIND), dimension(:), pointer :: &
        dvEdge, areaCell, vertCoordMovementWeights
      real (kind=RKIND), dimension(:), allocatable:: div_hu, h_tend_col
      real (kind=RKIND) :: div_hu_btr

      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, &
        verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex, &
        boundaryEdge, boundaryCell, edgeSignOnCell
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, &
        maxLevelCell, maxLevelEdgeTop, maxLevelEdgeBot, &
        maxLevelVertexBot,  maxLevelVertexTop

      err = 0

      nEdgesOnCell      => grid % nEdgesOnCell % array
      areaCell          => grid % areaCell % array
      cellsOnEdge       => grid % cellsOnEdge % array
      edgesOnCell       => grid % edgesOnCell % array
      edgeSignOnCell    => grid % edgeSignOnCell % array
      maxLevelCell      => grid % maxLevelCell % array
      maxLevelEdgeBot   => grid % maxLevelEdgeBot % array
      dvEdge            => grid % dvEdge % array
      vertCoordMovementWeights => grid % vertCoordMovementWeights % array

      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertLevels = grid % nVertLevels


      if (config_vert_coord_movement.eq.'isopycnal') then
        ! set vertical transport to zero in isopycnal case
        wTop=0.0
        return
      end if

      allocate(div_hu(nVertLevels), h_tend_col(nVertLevels))

      !
      ! Compute div(h^{edge} u) for each cell
      ! See Ringler et al. (2010) jcp paper, eqn 19, 21, and fig. 3.
      !

      do iCell=1,nCells
        div_hu(:) = 0.0
        div_hu_btr = 0.0
        hSum = 0.0
        invAreaCell = 1.0 / areaCell(iCell)

        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)

          do k = 1, maxLevelEdgeBot(iEdge)
            flux = h_edge(k, iEdge) * u(k, iEdge) * dvEdge(iEdge) * edgeSignOnCell(i, iCell) * invAreaCell
            div_hu(k) = div_hu(k) - flux
            div_hu_btr = div_hu_btr - flux
          end do
        end do

        do k = 1, maxLevelCell(iCell)
           h_tend_col(k) = - vertCoordMovementWeights(k) * h(k, iCell) * div_hu_btr
           hSum = hSum + vertCoordMovementWeights(k) * h(k, iCell)
        end do

        if(hSum > 0.0) then
           h_tend_col = h_tend_col / hSum
        end if

        ! Vertical transport through layer interface at top and bottom is zero.
        wTop(1,iCell) = 0.0
        wTop(maxLevelCell(iCell)+1,iCell) = 0.0
        do k=maxLevelCell(iCell),2,-1
           wTop(k,iCell) = wTop(k+1,iCell) - div_hu(k) - h_tend_col(k)
        end do
      end do

      deallocate(div_hu, h_tend_col)

   end subroutine ocn_wtop!}}}

!***********************************************************************
!
!  routine ocn_fuperp
!
!> \brief   Computes f u_perp
!> \author  Mark Petersen
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
!  routine ocn_diagnostics_init
!
!> \brief   Initializes flags used within diagnostics routines.
!> \author  Mark Petersen
!> \date    4 November 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes flags related to quantities computed within
!>  other diagnostics routines.
!
!-----------------------------------------------------------------------
    subroutine ocn_diagnostics_init(err)!{{{
        integer, intent(out) :: err !< Output: Error flag

        err = 0

        if(config_include_KE_vertex) then
            ke_vertex_flag = 1
            ke_cell_flag = 0
        else
            ke_vertex_flag = 0
            ke_cell_flag = 1
        endif

        if (trim(config_time_integrator) == 'RK4') then
            ! For RK4, PV includes f: PV = (eta+f)/h.
            fCoef = 1
        elseif (trim(config_time_integrator) == 'split_explicit' &
          .or.trim(config_time_integrator) == 'unsplit_explicit') then
            ! For split explicit, PV is eta/h because the Coriolis term 
            ! is added separately to the momentum tendencies.
            fCoef = 0
        end if

    end subroutine ocn_diagnostics_init!}}}

!***********************************************************************

end module ocn_diagnostics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
