module time_integration

   use grid_types
   use configure
   use constants
   use dmpar
   use vector_reconstruction
   use spline_interpolation

   contains

   subroutine timestep(domain, dt, timeStamp)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Advance model state forward in time by the specified time step
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
   !                 plus grid meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(in) :: timeStamp

      type (dm_info) :: dminfo
      type (block_type), pointer :: block

      if (trim(config_time_integration) == 'RK4') then
         call rk4(domain, dt)
      else
         write(0,*) 'Abort: Unknown time integration option '&
           //trim(config_time_integration)
         write(0,*) 'Currently, only ''RK4'' is supported.'
         call dmpar_abort(dminfo)
      end if

      block => domain % blocklist
      do while (associated(block))
         block % state % time_levs(2) % state % xtime % scalar = timeStamp

         if (isNaN(sum(block % state % time_levs(2) % state % u % array))) then
            write(0,*) 'Abort: NaN detected'
            call dmpar_abort(dminfo)
         endif

         block => block % next
      end do

   end subroutine timestep


   subroutine rk4(domain, dt)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Advance model state forward in time by the specified time step using 
   !   4th order Runge-Kutta
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
   !                 plus grid meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt

      integer :: iCell, k, i
      type (block_type), pointer :: block
      type (state_type) :: provis

      integer :: rk_step, iEdge

      real (kind=RKIND), dimension(4) :: rk_weights, rk_substep_weights

      block => domain % blocklist
      call allocate_state(provis, &
                          block % mesh % nCells, block % mesh % nEdges, block % mesh % maxEdges, block % mesh % maxEdges2, &
                          block % mesh % nVertices, block % mesh % vertexDegree, block % mesh % nVertLevels )

      !
      ! Initialize time_levs(2) with state at current time
      ! Initialize first RK state
      ! Couple tracers time_levs(2) with h in time-levels
      ! Initialize RK weights
      !
      block => domain % blocklist
      do while (associated(block))

         block % state % time_levs(2) % state % u % array(:,:) = block % state % time_levs(1) % state % u % array(:,:)
         block % state % time_levs(2) % state % h % array(:,:) = block % state % time_levs(1) % state % h % array(:,:)
         do iCell=1,block % mesh % nCells  ! couple tracers to h
           do k=1,block % mesh % maxLevelCell % array(iCell)
             block % state % time_levs(2) % state % tracers % array(:,k,iCell) = block % state % time_levs(1) % state % tracers % array(:,k,iCell) &
                                                                       * block % state % time_levs(1) % state % h % array(k,iCell)
            end do
         end do

         call copy_state(provis, block % state % time_levs(1) % state)

         block => block % next
      end do

      rk_weights(1) = dt/6.
      rk_weights(2) = dt/3.
      rk_weights(3) = dt/3.
      rk_weights(4) = dt/6.

      rk_substep_weights(1) = dt/2.
      rk_substep_weights(2) = dt/2.
      rk_substep_weights(3) = dt
      rk_substep_weights(4) = 0.


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN RK loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do rk_step = 1, 4
! ---  update halos for diagnostic variables

        block => domain % blocklist
        do while (associated(block))
           call dmpar_exch_halo_field2dReal(domain % dminfo, provis % pv_edge % array(:,:), &
                                            block % mesh % nVertLevels, block % mesh % nEdges, &
                                            block % parinfo % edgesToSend, block % parinfo % edgesToRecv)

           if (config_h_mom_eddy_visc4 > 0.0) then
              call dmpar_exch_halo_field2dReal(domain % dminfo, block % state % time_levs(2) % state % divergence % array(:,:), &
                                               block % mesh % nVertLevels, block % mesh % nCells, &
                                               block % parinfo % cellsToSend, block % parinfo % cellsToRecv)
              call dmpar_exch_halo_field2dReal(domain % dminfo, block % state % time_levs(2) % state % vorticity % array(:,:), &
                                               block % mesh % nVertLevels, block % mesh % nVertices, &
                                               block % parinfo % verticesToSend, block % parinfo % verticesToRecv)
           end if

           block => block % next
        end do

! ---  compute tendencies

        block => domain % blocklist
        do while (associated(block))
           call compute_tend(block % tend, provis, block % mesh)
           call compute_scalar_tend(block % tend, provis, block % mesh)
           call enforce_boundaryEdge(block % tend, block % mesh)
           block => block % next
        end do

! ---  update halos for prognostic variables

        block => domain % blocklist
        do while (associated(block))
           call dmpar_exch_halo_field2dReal(domain % dminfo, block % tend % u % array(:,:), &
                                            block % mesh % nVertLevels, block % mesh % nEdges, &
                                            block % parinfo % edgesToSend, block % parinfo % edgesToRecv)
           call dmpar_exch_halo_field2dReal(domain % dminfo, block % tend % h % array(:,:), &
                                            block % mesh % nVertLevels, block % mesh % nCells, &
                                            block % parinfo % cellsToSend, block % parinfo % cellsToRecv)
           call dmpar_exch_halo_field3dReal(domain % dminfo, block % tend % tracers % array(:,:,:), &
                                            block % tend % num_tracers, block % mesh % nVertLevels, block % mesh % nCells, &
                                            block % parinfo % cellsToSend, block % parinfo % cellsToRecv)
           block => block % next
        end do

! ---  compute next substep state

        if (rk_step < 4) then
           block => domain % blocklist
           do while (associated(block))

              provis % u % array(:,:)       = block % state % time_levs(1) % state % u % array(:,:)  &
                                         + rk_substep_weights(rk_step) * block % tend % u % array(:,:)
              provis % h % array(:,:)       = block % state % time_levs(1) % state % h % array(:,:)  &
                                         + rk_substep_weights(rk_step) * block % tend % h % array(:,:)
              do iCell=1,block % mesh % nCells
                 do k=1,block % mesh % maxLevelCell % array(iCell)
                    provis % tracers % array(:,k,iCell) = ( &
                                                                      block % state % time_levs(1) % state % h % array(k,iCell) * &
                                                                      block % state % time_levs(1) % state % tracers % array(:,k,iCell)  &
                                      + rk_substep_weights(rk_step) * block % tend % tracers % array(:,k,iCell) &
                                                                     ) / provis % h % array(k,iCell)
                 end do

              end do
              if (config_test_case == 1) then    ! For case 1, wind field should be fixed
                 provis % u % array(:,:) = block % state % time_levs(1) % state % u % array(:,:)
              end if
              call compute_solve_diagnostics(dt, provis, block % mesh)
              block => block % next
           end do
        end if



!--- accumulate update (for RK4)

        block => domain % blocklist
        do while (associated(block))
           block % state % time_levs(2) % state % u % array(:,:) = block % state % time_levs(2) % state % u % array(:,:) &
                                   + rk_weights(rk_step) * block % tend % u % array(:,:) 
           block % state % time_levs(2) % state % h % array(:,:) = block % state % time_levs(2) % state % h % array(:,:) &
                                   + rk_weights(rk_step) * block % tend % h % array(:,:) 

           do iCell=1,block % mesh % nCells
              do k=1,block % mesh % maxLevelCell % array(iCell)
                 block % state % time_levs(2) % state % tracers % array(:,k,iCell) =  &
                                                                       block % state % time_levs(2) % state % tracers % array(:,k,iCell) &
                                               + rk_weights(rk_step) * block % tend % tracers % array(:,k,iCell)
              end do
           end do

           block => block % next
        end do

      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END RK loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !
      !  A little clean up at the end: decouple new scalar fields and compute diagnostics for new state
      !
      block => domain % blocklist
      do while (associated(block))
         do iCell=1,block % mesh % nCells
            do k=1,block % mesh % maxLevelCell % array(iCell)
               block % state % time_levs(2) % state % tracers % array(:,k,iCell) = &
                                                                     block % state % time_levs(2) % state % tracers % array(:,k,iCell)  &
                                                                   / block % state % time_levs(2) % state % h % array(k,iCell)
            end do
         end do

         if (config_test_case == 1) then    ! For case 1, wind field should be fixed
            block % state % time_levs(2) % state % u % array(:,:) = block % state % time_levs(1) % state % u % array(:,:)
         end if

         call compute_solve_diagnostics(dt, block % state % time_levs(2) % state, block % mesh)

         call reconstruct(block % state % time_levs(2) % state, block % mesh)

         block => block % next
      end do

      call deallocate_state(provis)

   end subroutine rk4


   subroutine compute_tend(tend, s, grid)
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
      type (mesh_type), intent(in) :: grid

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, &
        vertex1, vertex2, eoe, i, j

      integer :: nCells, nEdges, nVertices, nVertLevels, nEdgesSolve
      real (kind=RKIND) :: flux, vorticity_abs, h_vertex, workpv, q, &
        upstream_bias, wTopEdge, rho0Inv, r
      real (kind=RKIND), dimension(:), pointer :: &
        h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle, &
        zMidZLevel, zTopZLevel 
      real (kind=RKIND), dimension(:,:), pointer :: &
        weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, pressure, &
        tend_h, tend_u, circulation, vorticity, ke, ke_edge, pv_edge, &
        MontPot, wTop, divergence
      type (dm_info) :: dminfo

      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, &
        maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot
      integer, dimension(:,:), pointer :: &
        cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, &
        edgesOnEdge, edgesOnVertex
      real (kind=RKIND) :: u_diffusion
      real (kind=RKIND), dimension(:), allocatable:: fluxVertTop,w_dudzTopEdge, vertViscTop

      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_divergence
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_u
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_circulation, delsq_vorticity


      real (kind=RKIND), dimension(:,:), pointer :: u_src
      real (kind=RKIND), parameter :: rho_ref = 1000.0

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
      tend_u      => tend % u % array
                  
      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nEdgesSolve = grid % nEdgesSolve
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels

      u_src => grid % u_src % array

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

      if (config_vert_grid_type.eq.'isopycnal') then
 
         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1,nVertLevels
               flux = u(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge)
               tend_h(k,cell1) = tend_h(k,cell1) - flux
               tend_h(k,cell2) = tend_h(k,cell2) + flux
            end do
         end do
         do iCell=1,nCells
            do k=1,nVertLevels
               tend_h(k,iCell) = tend_h(k,iCell) / areaCell(iCell)
            end do
         end do

      elseif (config_vert_grid_type.eq.'zlevel') then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1,min(1,maxLevelEdgeTop(iEdge))
               flux = u(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge)
               tend_h(k,cell1) = tend_h(k,cell1) - flux
               tend_h(k,cell2) = tend_h(k,cell2) + flux
            end do
         end do
         do iCell=1,nCells
            tend_h(1,iCell) = tend_h(1,iCell) / areaCell(iCell)
         end do

      endif ! config_vert_grid_type

      !
      ! height tendency: vertical advection term -d/dz(hw)
      !
      ! Vertical advection computed for top layer of a z grid only.
      if (config_vert_grid_type.eq.'zlevel') then
        do iCell=1,nCells
           tend_h(1,iCell) =   tend_h(1,iCell) + wTop(2,iCell)
        end do
      endif ! coordinate type

      !
      ! velocity tendency: start accumulating tendency terms
      !
      tend_u(:,:) = 0.0

      !
      ! velocity tendency: vertical advection term -w du/dz
      !
      if (config_vert_grid_type.eq.'zlevel') then
        allocate(w_dudzTopEdge(nVertLevels+1))
        w_dudzTopEdge(1) = 0.0
        do iEdge=1,grid % nEdgesSolve
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)

          do k=2,maxLevelEdgeTop(iEdge)
            ! Average w from cell center to edge
            wTopEdge = 0.5*(wTop(k,cell1)+wTop(k,cell2))

            ! compute dudz at vertical interface with first order derivative.
            w_dudzTopEdge(k) = wTopEdge * (u(k-1,iEdge)-u(k,iEdge)) &
                         / (zMidZLevel(k-1) - zMidZLevel(k))
          end do
          w_dudzTopEdge(maxLevelEdgeTop(iEdge)+1) = 0.0

          ! Average w*du/dz from vertical interface to vertical middle of cell
          do k=1,maxLevelEdgeTop(iEdge)
             tend_u(k,iEdge) = - 0.5 * (w_dudzTopEdge(k) + w_dudzTopEdge(k+1))
          enddo
        enddo
        deallocate(w_dudzTopEdge)
      endif

      !
      ! velocity tendency: pressure gradient
      !
      rho0Inv = 1.0/config_rho0
      if (config_vert_grid_type.eq.'isopycnal') then
        do iEdge=1,grid % nEdgesSolve
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)
          do k=1,maxLevelEdgeTop(iEdge)
             tend_u(k,iEdge) = tend_u(k,iEdge)     &
               - (MontPot(k,cell2) - MontPot(k,cell1))/dcEdge(iEdge)
           end do
        enddo
      elseif (config_vert_grid_type.eq.'zlevel') then
        do iEdge=1,grid % nEdgesSolve
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)
          do k=1,maxLevelEdgeTop(iEdge)
             tend_u(k,iEdge) = tend_u(k,iEdge)     &
               - rho0Inv*(  pressure(k,cell2) &
                          - pressure(k,cell1) )/dcEdge(iEdge)
          end do
        enddo
      endif

      !
      ! velocity tendency: del2 dissipation, \nu_2 \nabla^2 u
      !   computed as \nu( \nabla divergence + k \times \nabla vorticity )
      !   strictly only valid for config_h_mom_eddy_visc2 == constant
      !
      if ( config_h_mom_eddy_visc2 > 0.0 ) then
         do iEdge=1,grid % nEdgesSolve
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)

            do k=1,maxLevelEdgeTop(iEdge)

               ! Here -( vorticity(k,vertex2) - vorticity(k,vertex1) ) / dvEdge(iEdge)
               ! is - \nabla vorticity pointing from vertex 2 to vertex 1, or equivalently 
               !    + k \times \nabla vorticity pointing from cell1 to cell2.

               u_diffusion = ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
                            -( vorticity(k,vertex2) - vorticity(k,vertex1) ) / dvEdge(iEdge)
               u_diffusion = config_h_mom_eddy_visc2 * u_diffusion

               tend_u(k,iEdge) = tend_u(k,iEdge) + u_diffusion

            end do
         end do
      end if

      !
      ! velocity tendency: del4 dissipation, -\nu_4 \nabla^4 u
      !   computed as \nabla^2 u = \nabla divergence + k \times \nabla vorticity
      !   applied recursively.
      !   strictly only valid for config_h_mom_eddy_visc4 == constant
      !
      if ( config_h_mom_eddy_visc4 > 0.0 ) then

         allocate(delsq_divergence(nVertLevels, nCells+1))
         allocate(delsq_u(nVertLevels, nEdges+1))
         allocate(delsq_circulation(nVertLevels, nVertices+1))
         allocate(delsq_vorticity(nVertLevels, nVertices+1))

         delsq_u(:,:) = 0.0

         ! Compute \nabla^2 u = \nabla divergence + k \times \nabla vorticity
         do iEdge=1,grid % nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)

            do k=1,maxLevelEdgeTop(iEdge)

               delsq_u(k,iEdge) = & 
                  ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
                 -( vorticity(k,vertex2) - vorticity(k,vertex1)) / dvEdge(iEdge)

            end do
         end do

         ! vorticity using \nabla^2 u
         delsq_circulation(:,:) = 0.0
         do iEdge=1,nEdges
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)
            do k=1,maxLevelEdgeTop(iEdge)
               delsq_circulation(k,vertex1) = delsq_circulation(k,vertex1) &
                  - dcEdge(iEdge) * delsq_u(k,iEdge)
               delsq_circulation(k,vertex2) = delsq_circulation(k,vertex2) &
                  + dcEdge(iEdge) * delsq_u(k,iEdge)
            end do
         end do
         do iVertex=1,nVertices
            r = 1.0 / areaTriangle(iVertex)
            do k=1,maxLevelVertexBot(iVertex)
               delsq_vorticity(k,iVertex) = delsq_circulation(k,iVertex) * r
            end do
         end do

         ! Divergence using \nabla^2 u
         delsq_divergence(:,:) = 0.0
         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1,maxLevelEdgeTop(iEdge)
              delsq_divergence(k,cell1) = delsq_divergence(k,cell1) &
                + delsq_u(k,iEdge)*dvEdge(iEdge)
              delsq_divergence(k,cell2) = delsq_divergence(k,cell2) &
                - delsq_u(k,iEdge)*dvEdge(iEdge)
            end do
         end do
         do iCell = 1,nCells
            r = 1.0 / areaCell(iCell)
            do k = 1,maxLevelCell(iCell)
               delsq_divergence(k,iCell) = delsq_divergence(k,iCell) * r
            end do
         end do

         ! Compute - \kappa \nabla^4 u 
         ! as  \nabla div(\nabla^2 u) + k \times \nabla ( k \cross curl(\nabla^2 u) )
         do iEdge=1,grid % nEdgesSolve
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)

            do k=1,maxLevelEdgeTop(iEdge)

               u_diffusion = (  delsq_divergence(k,cell2) &
                              - delsq_divergence(k,cell1) ) / dcEdge(iEdge)  &
                            -(  delsq_vorticity(k,vertex2) &
                              - delsq_vorticity(k,vertex1) ) / dvEdge(iEdge)
 
               tend_u(k,iEdge) = tend_u(k,iEdge) - config_h_mom_eddy_visc4 * u_diffusion
            end do
         end do

         deallocate(delsq_divergence)
         deallocate(delsq_u)
         deallocate(delsq_circulation)
         deallocate(delsq_vorticity)

      end if

      !
      ! velocity tendency: nonlinear Coriolis term and grad of kinetic energy
      !
      do iEdge=1,grid % nEdgesSolve
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k=1,maxLevelEdgeTop(iEdge)

            q = 0.0
            do j = 1,nEdgesOnEdge(iEdge)
               eoe = edgesOnEdge(j,iEdge)
               workpv = 0.5 * (pv_edge(k,iEdge) + pv_edge(k,eoe))
               q = q + weightsOnEdge(j,iEdge) * u(k,eoe) * workpv * h_edge(k,eoe) 
            end do
            tend_u(k,iEdge) = tend_u(k,iEdge)     &
                   + q     &
                   - (   ke(k,cell2) - ke(k,cell1) ) / dcEdge(iEdge)

         end do
      end do

      !
      ! velocity tendency: forcing and bottom drag
      !
      ! mrp 101115 note: in order to include flux boundary conditions, we will need to 
      ! know the bottom edge with nonzero velocity and place the drag there.

      do iEdge=1,grid % nEdgesSolve

        k = maxLevelEdgeTop(iEdge)

        ! efficiency note: it would be nice to avoid this
        ! if within a do.  This could be done with
        ! k =  max(maxLevelEdgeTop(iEdge),1)
        ! and then tend_u(1,iEdge) is just not used for land cells.

        if (k>0) then

           ! forcing in top layer only
           tend_u(1,iEdge) =  tend_u(1,iEdge) &
              + u_src(1,iEdge)/rho_ref/h_edge(1,iEdge)

           ! bottom drag is the same as POP:
           ! -c |u| u  where c is unitless and 1.0e-3.
           ! see POP Reference guide, section 3.4.4.

           tend_u(k,iEdge) = tend_u(k,iEdge)  &
             - 1.0e-3*u(k,iEdge) &
               *sqrt(2.0*ke_edge(k,iEdge))/h_edge(k,iEdge)

        endif

      enddo

      !
      ! velocity tendency: vertical mixing d/dz( nu_v du/dz))
      !
      allocate(vertViscTop(nVertLevels+1))
      if (config_vert_visc_type.eq.'const') then
        vertViscTop = config_vert_viscosity
      elseif (config_vert_visc_type.eq.'tanh') then
        if (config_vert_grid_type.ne.'zlevel') then
          write(0,*) 'Abort: config_vert_visc_type.eq.tanh may only', &
            ' use config_vert_grid_type of zlevel at this time'
          call dmpar_abort(dminfo)
        endif
  
        do k=1,nVertLevels+1
          vertViscTop(k) = -(config_vmixTanhViscMax-config_vmixTanhViscMin)/2.0 &
            *tanh(-(zTopZLevel(k)-config_vmixTanhZMid) &
                  /config_vmixTanhZWidth) &
            + (config_vmixTanhViscMax+config_vmixTanhViscMin)/2
        enddo
      else
        write(0,*) 'Abort: unrecognized config_vert_visc_type'
        call dmpar_abort(dminfo)
      endif

      allocate(fluxVertTop(nVertLevels+1))
      fluxVertTop(1) = 0.0
      do iEdge=1,grid % nEdgesSolve
         
         do k=2,maxLevelEdgeTop(iEdge)
           fluxVertTop(k) = vertViscTop(k) &
              * ( u(k-1,iEdge) - u(k,iEdge) ) &
              * 2 / (h_edge(k-1,iEdge) + h_edge(k,iEdge))
         enddo
         fluxVertTop(maxLevelEdgeTop(iEdge)+1) = 0.0

         do k=1,maxLevelEdgeTop(iEdge)
           tend_u(k,iEdge) = tend_u(k,iEdge) &
             + (fluxVertTop(k) - fluxVertTop(k+1)) &
             / h_edge(k,iEdge)
         enddo

      end do
      deallocate(fluxVertTop, vertViscTop)

   end subroutine compute_tend


   subroutine compute_scalar_tend(tend, s, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed scalar tendencies
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (tend_type), intent(inout) :: tend
      type (state_type), intent(in) :: s
      type (mesh_type), intent(in) :: grid

      integer :: i, k, iCell, iEdge, iTracer, cell1, cell2, upwindCell,&
        nEdges, nCells, nCellsSolve, nVertLevels, num_tracers
      real (kind=RKIND) :: invAreaCell1, invAreaCell2, tracer_turb_flux
      real (kind=RKIND) :: flux, tracer_edge, r
      real (kind=RKIND), dimension(:), pointer :: &
        h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: &
        u,h,wTop, h_edge
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        tracers, tend_tr
      integer, dimension(:,:), pointer :: boundaryEdge
      type (dm_info) :: dminfo

      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge, &
        maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot
      integer, dimension(:,:), pointer :: cellsOnEdge, boundaryCell
      real (kind=RKIND), dimension(:), pointer :: zTopZLevel,zMidZLevel, &
         hRatioZLevelK, hRatioZLevelKm1
      real (kind=RKIND), dimension(:), allocatable:: vertDiffTop, tracer2ndDer, tracersIn, tracersOut, posZMidZLevel, &
            posZTopZLevel
      real (kind=RKIND), dimension(:,:), allocatable:: fluxVertTop, boundaryMask
      real (kind=RKIND), dimension(:,:,:), allocatable::tr_flux, tr_div, delsq_tracer, tracerTop


      real (kind=RKIND) :: d2fdx2_cell1, d2fdx2_cell2
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND) :: coef_3rd_order, flux3Coef, cSignWTop

      u           => s % u % array
      h           => s % h % array
      boundaryCell=> grid % boundaryCell % array
      wTop        => s % wTop % array
      tracers     => s % tracers % array
      h_edge      => s % h_edge % array

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

      coef_3rd_order = 0.
      if (config_tracer_adv_order == 3) coef_3rd_order = 1.0
      if (config_tracer_adv_order == 3 .and. config_monotonic) coef_3rd_order = 0.25

      if (config_tracer_adv_order == 2) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1,maxLevelEdgeTop(iEdge)
               do iTracer=1,num_tracers
                  tracer_edge = 0.5 * (tracers(iTracer,k,cell1) + tracers(iTracer,k,cell2))
                  flux = u(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge) * tracer_edge
                  tend_tr(iTracer,k,cell1) = tend_tr(iTracer,k,cell1) - flux/areaCell(cell1)
                  tend_tr(iTracer,k,cell2) = tend_tr(iTracer,k,cell2) + flux/areaCell(cell2)
               end do
            end do
         end do

      else if (config_tracer_adv_order == 3) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            do k=1,maxLevelEdgeTop(iEdge)

               d2fdx2_cell1 = 0.0
               d2fdx2_cell2 = 0.0

               do iTracer=1,num_tracers

                  !-- if not a boundary cell
                  if(boundaryCell(k,cell1).eq.0.and.boundaryCell(k,cell2).eq.0) then

                     d2fdx2_cell1 = deriv_two(1,1,iEdge) * tracers(iTracer,k,cell1)
                     d2fdx2_cell2 = deriv_two(1,2,iEdge) * tracers(iTracer,k,cell2)

                     !-- all edges of cell 1
                     do i=1, grid % nEdgesOnCell % array (cell1)
                        d2fdx2_cell1 = d2fdx2_cell1 + &
                        deriv_two(i+1,1,iEdge) * tracers(iTracer,k,grid % CellsOnCell % array (i,cell1))
                     end do

                     !-- all edges of cell 2
                     do i=1, grid % nEdgesOnCell % array (cell2)
                        d2fdx2_cell2 = d2fdx2_cell2 + &
                        deriv_two(i+1,2,iEdge) * tracers(iTracer,k,grid % CellsOnCell % array (i,cell2))
                     end do

                  endif

                  !-- if u > 0:
                  if (u(k,iEdge) > 0) then
                     flux = dvEdge(iEdge) * u(k,iEdge) * h_edge(k,iEdge) * (          &
                          0.5*(tracers(iTracer,k,cell1) + tracers(iTracer,k,cell2))      &
                          -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                          -(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12. )
                  !-- else u <= 0:
                  else
                     flux = dvEdge(iEdge) *  u(k,iEdge) * h_edge(k,iEdge) * (          &
                          0.5*(tracers(iTracer,k,cell1) + tracers(iTracer,k,cell2))      &
                          -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                          +(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12. )
                  end if

                  !-- update tendency
                  tend_tr(iTracer,k,cell1) = tend_tr(iTracer,k,cell1) - flux/areaCell(cell1)
                  tend_tr(iTracer,k,cell2) = tend_tr(iTracer,k,cell2) + flux/areaCell(cell2)
               enddo
            end do
         end do

      else if (config_tracer_adv_order == 4) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            do k=1,maxLevelEdgeTop(iEdge)

               d2fdx2_cell1 = 0.0
               d2fdx2_cell2 = 0.0

               do iTracer=1,num_tracers

                  !-- if not a boundary cell
                  if(boundaryCell(k,cell1).eq.0.and.boundaryCell(k,cell2).eq.0) then

                     d2fdx2_cell1 = deriv_two(1,1,iEdge) * tracers(iTracer,k,cell1)
                     d2fdx2_cell2 = deriv_two(1,2,iEdge) * tracers(iTracer,k,cell2)

                     !-- all edges of cell 1
                     do i=1, grid % nEdgesOnCell % array (cell1)
                        d2fdx2_cell1 = d2fdx2_cell1 + &
                        deriv_two(i+1,1,iEdge) * tracers(iTracer,k,grid % CellsOnCell % array (i,cell1))
                     end do

                     !-- all edges of cell 2
                     do i=1, grid % nEdgesOnCell % array (cell2)
                         d2fdx2_cell2 = d2fdx2_cell2 + &
                         deriv_two(i+1,2,iEdge) * tracers(iTracer,k,grid % CellsOnCell % array (i,cell2))
                     end do

                  endif

                  flux = dvEdge(iEdge) *  u(k,iEdge) * h_edge(k,iEdge) * (          &
                       0.5*(tracers(iTracer,k,cell1) + tracers(iTracer,k,cell2))      &
                          -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12. )

                  !-- update tendency
                  tend_tr(iTracer,k,cell1) = tend_tr(iTracer,k,cell1) - flux/areaCell(cell1)
                  tend_tr(iTracer,k,cell2) = tend_tr(iTracer,k,cell2) + flux/areaCell(cell2)
               enddo
            end do
         end do

      endif   ! if (config_tracer_adv_order == 2 )


      !
      ! tracer tendency: vertical advection term -d/dz( h \phi w)
      !

      if (config_vert_grid_type.eq.'zlevel') then

         allocate(tracerTop(num_tracers,nVertLevels+1,nCells))

         ! Tracers at the top and bottom boundary are assigned nearest 
         ! cell-centered value, regardless of tracer interpolation method.
         ! wTop=0 at top and bottom sets the boundary condition.
         do iCell=1,nCellsSolve 
            do iTracer=1,num_tracers
               tracerTop(iTracer,1,iCell) = tracers(iTracer,1,iCell)
               tracerTop(iTracer,maxLevelCell(iCell)+1,iCell) = &
               tracers(iTracer,maxLevelCell(iCell),iCell)
            end do
         end do

         if (config_vert_tracer_adv.eq.'stencil'.and. &
             config_vert_tracer_adv_order.eq.2) then

            ! Compute tracerTop using centered stencil, a simple average.

            do iCell=1,nCellsSolve 
               do k=2,maxLevelCell(iCell)
                  do iTracer=1,num_tracers
                     tracerTop(iTracer,k,iCell) = &
                        ( tracers(iTracer,k-1,iCell) &
                         +tracers(iTracer,k  ,iCell))/2.0
                  end do
               end do
            end do
         
         elseif (config_vert_tracer_adv.eq.'stencil'.and. &
             config_vert_tracer_adv_order.eq.3) then

            ! Compute tracerTop using 3rd order stencil.  This is the same
            ! as 4th order, but includes upwinding.

            ! Hardwire flux3Coeff at 1.0 for now.  Could add this to the 
            ! namelist, if desired.
            flux3Coef = 1.0
            do iCell=1,nCellsSolve 
               k=2
               do iTracer=1,num_tracers
                 tracerTop(iTracer,k,iCell) = &
                      hRatioZLevelK(k)  *tracers(iTracer,k-1,iCell) &
                    + hRatioZLevelKm1(k)*tracers(iTracer,k  ,iCell)
               end do
               do k=3,maxLevelCell(iCell)-1
                  cSignWTop = sign(flux3Coef,wTop(k,iCell))
                  do iTracer=1,num_tracers
                     tracerTop(iTracer,k,iCell) = &
                        ( (-1.+   cSignWTop)*tracers(iTracer,k-2,iCell) &
                         +( 7.-3.*cSignWTop)*tracers(iTracer,k-1,iCell) &
                         +( 7.+3.*cSignWTop)*tracers(iTracer,k  ,iCell) &
                         +(-1.-   cSignWTop)*tracers(iTracer,k+1,iCell) &
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

         elseif (config_vert_tracer_adv.eq.'stencil'.and. &
             config_vert_tracer_adv_order.eq.4) then

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

         elseif (config_vert_tracer_adv.eq.'spline'.and. &
             config_vert_tracer_adv_order.eq.2) then

            ! Compute tracerTop using linear interpolation.

            do iCell=1,nCellsSolve 
               do k=2,maxLevelCell(iCell)
                  do iTracer=1,num_tracers
                     ! Note hRatio on the k side is multiplied by tracer at k-1
                     ! and hRatio on the Km1 (k-1) side is mult. by tracer at k.
                     tracerTop(iTracer,k,iCell) = &
                          hRatioZLevelK(k)  *tracers(iTracer,k-1,iCell) &
                        + hRatioZLevelKm1(k)*tracers(iTracer,k  ,iCell)
                  end do
               end do
            end do
         
         elseif (config_vert_tracer_adv.eq.'spline'.and. &
             config_vert_tracer_adv_order.eq.3) then

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

                  call CubicSplineCoefficients(posZMidZLevel, &
                     tracersIn, maxLevelCell(iCell), tracer2ndDer)

                  call InterpolateCubicSpline( &
                     posZMidZLevel, tracersIn, tracer2ndDer, maxLevelCell(iCell), &
                     posZTopZLevel, tracersOut, maxLevelCell(iCell)-1 )

                  tracerTop(iTracer,2:maxLevelCell(iCell),iCell) = tracersOut(1:maxLevelCell(iCell)-1)

               end do
            end do

            deallocate(tracer2ndDer)
            deallocate(tracersIn,tracersOut, posZMidZLevel, posZTopZLevel)

        else

            print *, 'Abort: Incorrect combination of ', &
              'config_vert_tracer_adv and config_vert_tracer_adv_order.'
            print *, 'Use:'
            print *, 'config_vert_tracer_adv=''stencil'' and config_vert_tracer_adv_order=2,3,4 or'
            print *, 'config_vert_tracer_adv=''spline''  and config_vert_tracer_adv_order=2,3'
            call dmpar_abort(dminfo)

         endif ! vertical tracer advection method

         do iCell=1,nCellsSolve 
            do k=1,maxLevelCell(iCell)  
               do iTracer=1,num_tracers
                  tend_tr(iTracer,k,iCell) = tend_tr(iTracer,k,iCell) &
                     - (   wTop(k  ,iCell)*tracerTop(iTracer,k  ,iCell) &
                         - wTop(k+1,iCell)*tracerTop(iTracer,k+1,iCell))
               end do
            end do
         end do

         deallocate(tracerTop)

      endif ! ZLevel

      !
      ! tracer tendency: del2 horizontal tracer diffusion, div(h \kappa_2 \nabla \phi)
      !
      if ( config_h_tracer_eddy_diff2 > 0.0 ) then

         !
         ! compute a boundary mask to enforce insulating boundary conditions in the horizontal
         !
         allocate(boundaryMask(nVertLevels, nEdges+1))
         boundaryMask = 1.0
         where(boundaryEdge.eq.1) boundaryMask=0.0

         do iEdge=1,grid % nEdges
            cell1 = grid % cellsOnEdge % array(1,iEdge)
            cell2 = grid % cellsOnEdge % array(2,iEdge)
            invAreaCell1 = 1.0/areaCell(cell1)
            invAreaCell2 = 1.0/areaCell(cell2)

            do k=1,maxLevelEdgeTop(iEdge)
              do iTracer=1,num_tracers
                 ! \kappa_2 \nabla \phi on edge
                 tracer_turb_flux = config_h_tracer_eddy_diff2 &
                    *(  tracers(iTracer,k,cell2) &
                      - tracers(iTracer,k,cell1))/dcEdge(iEdge)

                 ! div(h \kappa_2 \nabla \phi) at cell center
                 flux = dvEdge (iEdge) * h_edge(k,iEdge) &
                    * tracer_turb_flux * boundaryMask(k, iEdge)
                 tend_tr(iTracer,k,cell1) = tend_tr(iTracer,k,cell1) + flux * invAreaCell1
                 tend_tr(iTracer,k,cell2) = tend_tr(iTracer,k,cell2) - flux * invAreaCell2
              end do
            end do

         end do

        deallocate(boundaryMask)

      end if

      !
      ! tracer tendency: del4 horizontal tracer diffusion, &
      !    div(h \kappa_4 \nabla [div(h \nabla \phi)])
      !
      if ( config_h_tracer_eddy_diff4 > 0.0 ) then

         !
         ! compute a boundary mask to enforce insulating boundary conditions in the horizontal
         !
         allocate(boundaryMask(nVertLevels, nEdges+1))
         boundaryMask = 1.0
         where(boundaryEdge.eq.1) boundaryMask=0.0

         allocate(delsq_tracer(num_tracers,nVertLevels, nCells+1))

         delsq_tracer(:,:,:) = 0.

         ! first del2: div(h \nabla \phi) at cell center
         do iEdge=1,grid % nEdges
            cell1 = grid % cellsOnEdge % array(1,iEdge)
            cell2 = grid % cellsOnEdge % array(2,iEdge)

            do k=1,maxLevelEdgeTop(iEdge)
              do iTracer=1,num_tracers
                 delsq_tracer(iTracer,k,cell1) = delsq_tracer(iTracer,k,cell1) &
                    + dvEdge(iEdge)*h_edge(k,iEdge) &
                      *(tracers(iTracer,k,cell2) - tracers(iTracer,k,cell1)) &
                      /dcEdge(iEdge) * boundaryMask(k,iEdge)
                 delsq_tracer(iTracer,k,cell2) = delsq_tracer(iTracer,k,cell2) &
                    - dvEdge(iEdge)*h_edge(k,iEdge) &
                    *(tracers(iTracer,k,cell2) - tracers(iTracer,k,cell1)) &
                    /dcEdge(iEdge) * boundaryMask(k,iEdge)
              end do
            end do

         end do

         do iCell = 1,nCells
            r = 1.0 / areaCell(iCell)
            do k=1,maxLevelCell(iCell)
            do iTracer=1,num_tracers
               delsq_tracer(iTracer,k,iCell) = delsq_tracer(iTracer,k,iCell) * r
            end do
            end do
         end do

         ! second del2: div(h \nabla [delsq_tracer]) at cell center
         do iEdge=1,grid % nEdges
            cell1 = grid % cellsOnEdge % array(1,iEdge)
            cell2 = grid % cellsOnEdge % array(2,iEdge)
            invAreaCell1 = 1.0 / areaCell(cell1)
            invAreaCell2 = 1.0 / areaCell(cell2)

            do k=1,maxLevelEdgeTop(iEdge)
               do iTracer=1,num_tracers
                  tracer_turb_flux = config_h_tracer_eddy_diff4 &
                     *(  delsq_tracer(iTracer,k,cell2)  &
                       - delsq_tracer(iTracer,k,cell1))/dcEdge(iEdge)
                  flux = dvEdge (iEdge) * tracer_turb_flux

                  tend_tr(iTracer,k,cell1) = tend_tr(iTracer,k,cell1) & 
                     - flux * invAreaCell1 * boundaryMask(k,iEdge)
                  tend_tr(iTracer,k,cell2) = tend_tr(iTracer,k,cell2) &
                     + flux * invAreaCell2 * boundaryMask(k,iEdge)

               enddo
            enddo
         end do

         deallocate(delsq_tracer)

      end if

      !
      ! tracer tendency: vertical diffusion h d/dz( \kappa_v d\phi/dz)
      !
      allocate(vertDiffTop(nVertLevels+1))
      if (config_vert_diff_type.eq.'const') then
        vertDiffTop = config_vert_diffusion
      elseif (config_vert_diff_type.eq.'tanh') then
        if (config_vert_grid_type.ne.'zlevel') then
          write(0,*) 'Abort: config_vert_diff_type.eq.tanh may only', &
            ' use config_vert_grid_type of zlevel at this time'
          call dmpar_abort(dminfo)
        endif
  
        do k=1,nVertLevels+1
          vertDiffTop(k) = -(config_vmixTanhDiffMax-config_vmixTanhDiffMin)/2.0 &
            *tanh(-(zTopZLevel(k)-config_vmixTanhZMid) &
                  /config_vmixTanhZWidth) &
            + (config_vmixTanhDiffMax+config_vmixTanhDiffMin)/2
        enddo
      else
        write(0,*) 'Abort: unrecognized config_vert_diff_type'
        call dmpar_abort(dminfo)
      endif

      allocate(fluxVertTop(num_tracers,nVertLevels+1))
      fluxVertTop(:,1) = 0.0
      do iCell=1,nCellsSolve 

         do k=2,maxLevelCell(iCell)
           do iTracer=1,num_tracers
             ! compute \kappa_v d\phi/dz
             fluxVertTop(iTracer,k) = vertDiffTop(k) &
                * (tracers(iTracer,k-1,iCell) - tracers(iTracer,k,iCell) )&
                * 2 / (h(k-1,iCell) + h(k,iCell))
           enddo
         enddo
         fluxVertTop(:,maxLevelCell(iCell)+1) = 0.0

         do k=1,maxLevelCell(iCell)
           do iTracer=1,num_tracers
             ! This is h d/dz( fluxVertTop) but h and dz cancel, so 
             ! reduces to delta( fluxVertTop)
             tend_tr(iTracer,k,iCell) = tend_tr(iTracer,k,iCell) &
               + fluxVertTop(iTracer,k) - fluxVertTop(iTracer,k+1)
           enddo
         enddo

      enddo ! iCell loop
      deallocate(fluxVertTop, vertDiffTop)

          ! print some diagnostics - for debugging
!         print *, 'after vertical mixing',&
! 'iTracer,k, minval(tend_tr(itracer,k,:)), maxval(tend_tr(itracer,k,:))'
!         do iTracer=1,num_tracers
!         do k = 1,nVertLevels
!            print '(2i5,20es12.4)', iTracer,k, &
!              minval(tend_tr(itracer,k,:)), maxval(tend_tr(itracer,k,:))
!         enddo
!         enddo


   end subroutine compute_scalar_tend


   subroutine compute_solve_diagnostics(dt, s, grid)
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

      integer :: nCells, nEdges, nVertices, nVertLevels, vertexDegree


      real (kind=RKIND), dimension(:), pointer :: &
        h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle, &
        hZLevel
      real (kind=RKIND), dimension(:,:), pointer :: &
        weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, w, pressure,&
        circulation, vorticity, ke, ke_edge, MontPot, wTop, &
        pv_edge, pv_vertex, pv_cell, gradPVn, gradPVt, divergence, &
        rho, temperature, salinity
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

      coef_3rd_order = 0.
      if (config_thickness_adv_order == 3) coef_3rd_order = 1.0
      if (config_thickness_adv_order == 3 .and. config_monotonic) coef_3rd_order = 0.25

      if (config_thickness_adv_order == 2) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1,maxLevelEdgeTop(iEdge)
               h_edge(k,iEdge) = 0.5 * (h(k,cell1) + h(k,cell2))
            end do
         end do

      else if (config_thickness_adv_order == 3) then

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

      else  if (config_thickness_adv_order == 4) then

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

      endif   ! if(config_thickness_adv_order == 2)

      !
      ! set the velocity and height at dummy address
      !    used -1e34 so error clearly occurs if these values are used.
      !
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
      do iVertex = 1,nVertices
         do k=1,maxLevelVertexBot(iVertex)
            h_vertex = 0.0
            do i=1,vertexDegree
               h_vertex = h_vertex + h(k,cellsOnVertex(i,iVertex)) * kiteAreasOnVertex(i,iVertex)
            end do
            h_vertex = h_vertex / areaTriangle(iVertex)

            pv_vertex(k,iVertex) = (fVertex(iVertex) + vorticity(k,iVertex)) / h_vertex
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
      if (config_vert_grid_type.eq.'zlevel') then
         call equation_of_state(s,grid)
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
                 - div_u(k,iCell)/areaCell(iCell)*h(k,iCell)
           end do
        end do
        deallocate(div_u)

      endif

   end subroutine compute_solve_diagnostics


   subroutine enforce_boundaryEdge(tend, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Enforce any boundary conditions on the normal velocity at each edge
   !
   ! Input: grid - grid metadata
   !
   ! Output: tend_u set to zero at boundaryEdge == 1 locations
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      implicit none

      type (tend_type), intent(inout) :: tend
      type (mesh_type), intent(in) :: grid

      integer, dimension(:,:), pointer :: boundaryEdge
      real (kind=RKIND), dimension(:,:), pointer :: tend_u
      integer :: nCells, nEdges, nVertices, nVertLevels
      integer :: iEdge, k

      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels

      boundaryEdge         => grid % boundaryEdge % array
      tend_u      => tend % u % array

      if(maxval(boundaryEdge).le.0) return

      do iEdge = 1,nEdges
        do k = 1,nVertLevels

          if(boundaryEdge(k,iEdge).eq.1) then
             tend_u(k,iEdge) = 0.0
          endif

        enddo
       enddo

   end subroutine enforce_boundaryEdge


   subroutine equation_of_state(s, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This module contains routines necessary for computing the density
   !  from model temperature and salinity using an equation of state.
   !
   ! Input: grid - grid metadata
   !        s - state: tracers
   !
   ! Output: s - state: computed density
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none

      type (state_type), intent(inout) :: s
      type (mesh_type), intent(in) :: grid

      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:,:), pointer :: rho
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      integer :: nCells, iCell, k
      type (dm_info) :: dminfo

      rho         => s % rho % array
      tracers     => s % tracers % array

      maxLevelCell      => grid % maxLevelCell % array
      nCells      = grid % nCells

      if (config_eos_type.eq.'linear') then

         do iCell=1,nCells
            do k=1,maxLevelCell(iCell)
               ! Linear equation of state
               rho(k,iCell) = 1000.0*(  1.0 &
                  - 2.5e-4*tracers(s % index_temperature,k,iCell) &
                  + 7.6e-4*tracers(s % index_salinity,k,iCell))
            end do
         end do

      elseif (config_eos_type.eq.'jm') then

         call equation_of_state_jm(s, grid)  

      else 
         print *, ' Incorrect choice of config_eos_type:',&
            config_eos_type
         call dmpar_abort(dminfo)
      endif

   end subroutine equation_of_state


   subroutine equation_of_state_jm(s, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This module contains routines necessary for computing the density
   !  from model temperature and salinity using an equation of state.
   !
   !  The UNESCO equation of state computed using the
   !  potential-temperature-based bulk modulus from Jackett and
   !  McDougall, JTECH, Vol.12, pp 381-389, April, 1995.
   !
   ! Input: grid - grid metadata
   !        s - state: tracers
   !
   ! Output: s - state: computed density
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (state_type), intent(in) :: s
      type (mesh_type), intent(in) :: grid


      type (dm_info) :: dminfo
      integer :: iEdge, iCell, iVertex, k

      integer :: nCells, nEdges, nVertices, nVertLevels


      real (kind=RKIND), dimension(:), pointer :: &
        zMidZLevel, pRefEOS
      real (kind=RKIND), dimension(:,:), pointer :: &
        rho
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers

      integer, dimension(:), pointer :: maxLevelCell

   real (kind=RKIND) :: &
      TQ,SQ,             &! adjusted T,S
      BULK_MOD,          &! Bulk modulus
      RHO_S,             &! density at the surface
      DRDT0,             &! d(density)/d(temperature), for surface
      DRDS0,             &! d(density)/d(salinity   ), for surface
      DKDT,              &! d(bulk modulus)/d(pot. temp.)
      DKDS,              &! d(bulk modulus)/d(salinity  )
      SQR,DENOMK,        &! work arrays
      WORK1, WORK2, WORK3, WORK4, T2, depth

   real (kind=RKIND) :: & 
      tmin, tmax,        &! valid temperature range for level k
      smin, smax          ! valid salinity    range for level k

   real (kind=RKIND) :: p, p2 ! temporary pressure scalars

!-----------------------------------------------------------------------
!
!  UNESCO EOS constants and JMcD bulk modulus constants
!
!-----------------------------------------------------------------------

   !*** for density of fresh water (standard UNESCO)

   real (kind=RKIND), parameter ::              &
      unt0 =   999.842594,           &
      unt1 =  6.793952e-2,           &
      unt2 = -9.095290e-3,           &
      unt3 =  1.001685e-4,           &
      unt4 = -1.120083e-6,           &
      unt5 =  6.536332e-9

   !*** for dependence of surface density on salinity (UNESCO)

   real (kind=RKIND), parameter ::              &
      uns1t0 =  0.824493 ,           &
      uns1t1 = -4.0899e-3,           &
      uns1t2 =  7.6438e-5,           &
      uns1t3 = -8.2467e-7,           &
      uns1t4 =  5.3875e-9,           &
      unsqt0 = -5.72466e-3,          &
      unsqt1 =  1.0227e-4,           &
      unsqt2 = -1.6546e-6,           &
      uns2t0 =  4.8314e-4

   !*** from Table A1 of Jackett and McDougall

   real (kind=RKIND), parameter ::              &
      bup0s0t0 =  1.965933e+4,       &
      bup0s0t1 =  1.444304e+2,       &
      bup0s0t2 = -1.706103   ,       &
      bup0s0t3 =  9.648704e-3,       &
      bup0s0t4 = -4.190253e-5

   real (kind=RKIND), parameter ::              &
      bup0s1t0 =  5.284855e+1,       &
      bup0s1t1 = -3.101089e-1,       &
      bup0s1t2 =  6.283263e-3,       &
      bup0s1t3 = -5.084188e-5

   real (kind=RKIND), parameter ::              &
      bup0sqt0 =  3.886640e-1,       &
      bup0sqt1 =  9.085835e-3,       &
      bup0sqt2 = -4.619924e-4

   real (kind=RKIND), parameter ::              &
      bup1s0t0 =  3.186519   ,       &
      bup1s0t1 =  2.212276e-2,       &
      bup1s0t2 = -2.984642e-4,       &
      bup1s0t3 =  1.956415e-6 

   real (kind=RKIND), parameter ::              &
      bup1s1t0 =  6.704388e-3,       &
      bup1s1t1 = -1.847318e-4,       &
      bup1s1t2 =  2.059331e-7,       &
      bup1sqt0 =  1.480266e-4 

   real (kind=RKIND), parameter ::              &
      bup2s0t0 =  2.102898e-4,       &
      bup2s0t1 = -1.202016e-5,       &
      bup2s0t2 =  1.394680e-7,       &
      bup2s1t0 = -2.040237e-6,       &
      bup2s1t1 =  6.128773e-8,       &
      bup2s1t2 =  6.207323e-10

      rho         => s % rho % array
      tracers     => s % tracers % array

      nCells      = grid % nCells
      maxLevelCell      => grid % maxLevelCell % array
      nVertLevels = grid % nVertLevels
      zMidZLevel        => grid % zMidZLevel % array


!  Jackett and McDougall
      tmin = -2.0  ! valid pot. temp. range
      tmax = 40.0 
      smin =  0.0  ! valid salinity, in psu   
      smax = 42.0 

      ! This could be put in a startup routine.
      ! Note I am using zMidZlevel, so pressure on top level does
      ! not include SSH contribution.  I am not sure if that matters.

!  This function computes pressure in bars from depth in meters
!  using a mean density derived from depth-dependent global 
!  average temperatures and salinities from Levitus 1994, and 
!  integrating using hydrostatic balance.

      allocate(pRefEOS(nVertLevels))
      do k = 1,nVertLevels
        depth = -zMidZLevel(k)
        pRefEOS(k) = 0.059808*(exp(-0.025*depth) - 1.0) &
            + 0.100766*depth + 2.28405e-7*depth**2
      enddo

  do iCell=1,nCells
    do k=1,maxLevelCell(iCell)

      p   = pRefEOS(k)
      p2  = pRefEOS(k)*pRefEOS(k)

      SQ  = max(min(tracers(s%index_salinity,k,iCell),smax),smin)
      TQ  = max(min(tracers(s%index_temperature,k,iCell),tmax),tmin)

      SQR = sqrt(SQ)
      T2  = TQ*TQ

      !***
      !*** first calculate surface (p=0) values from UNESCO eqns.
      !***

      WORK1 = uns1t0 + uns1t1*TQ + & 
             (uns1t2 + uns1t3*TQ + uns1t4*T2)*T2
      WORK2 = SQR*(unsqt0 + unsqt1*TQ + unsqt2*T2)

      RHO_S = unt1*TQ + (unt2 + unt3*TQ + (unt4 + unt5*TQ)*T2)*T2 &
                      + (uns2t0*SQ + WORK1 + WORK2)*SQ

      !***
      !*** now calculate bulk modulus at pressure p from 
      !*** Jackett and McDougall formula
      !***

      WORK3 = bup0s1t0 + bup0s1t1*TQ +                    &
             (bup0s1t2 + bup0s1t3*TQ)*T2 +                &
              p *(bup1s1t0 + bup1s1t1*TQ + bup1s1t2*T2) + &
              p2*(bup2s1t0 + bup2s1t1*TQ + bup2s1t2*T2)
      WORK4 = SQR*(bup0sqt0 + bup0sqt1*TQ + bup0sqt2*T2 + &
                   bup1sqt0*p)

      BULK_MOD  = bup0s0t0 + bup0s0t1*TQ +                    &
                  (bup0s0t2 + bup0s0t3*TQ + bup0s0t4*T2)*T2 + &
                  p *(bup1s0t0 + bup1s0t1*TQ +                &
                     (bup1s0t2 + bup1s0t3*TQ)*T2) +           &
                  p2*(bup2s0t0 + bup2s0t1*TQ + bup2s0t2*T2) + &
                  SQ*(WORK3 + WORK4)

      DENOMK = 1.0/(BULK_MOD - p)

      rho(k,iCell) = (unt0 + RHO_S)*BULK_MOD*DENOMK

    end do
  end do

      deallocate(pRefEOS)

   end subroutine equation_of_state_jm

end module time_integration
