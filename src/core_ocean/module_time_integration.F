module time_integration

   use grid_types
   use configure
   use constants
   use dmpar
   use vector_reconstruction

   contains

   subroutine timestep(domain, dt)
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
         block % state % time_levs(2) % state % xtime % scalar &
           = block % state % time_levs(1) % state % xtime % scalar + dt

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

      integer :: rk_step

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
           do k=1,block % mesh % nVertLevels
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
                 do k=1,block % mesh % nVertLevels
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
              do k=1,block % mesh % nVertLevels
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
            do k=1,block % mesh % nVertLevels
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

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j

      integer :: nCells, nEdges, nVertices, nVertLevels
      real (kind=RKIND) :: h_mom_eddy_visc2, h_mom_eddy_visc4
      real (kind=RKIND) :: flux, vorticity_abs, h_vertex, workpv, q, &
        upstream_bias, wTopEdge, rho0Inv, r
      real (kind=RKIND), dimension(:), pointer :: &
        h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle, &
        zMidZLevel, zTopZLevel 
      real (kind=RKIND), dimension(:,:), pointer :: &
        weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, wTop, &
        tend_h, tend_u, circulation, vorticity, ke, ke_edge, pv_edge, &
        divergence, MontPot, pZLevel, zMidEdge, zTopEdge
      type (dm_info) :: dminfo

      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge
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
      ke_edge          => s % ke_edge % array
      pv_edge     => s % pv_edge % array
      MontPot     => s % MontPot % array
      pZLevel     => s % pZLevel % array
      zTopEdge    => s % zTopEdge % array
      zMidEdge    => s % zMidEdge % array

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

      tend_h      => tend % h % array
      tend_u      => tend % u % array
                  
      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels

      u_src => grid % u_src % array

      h_mom_eddy_visc2 = config_h_mom_eddy_visc2
      h_mom_eddy_visc4 = config_h_mom_eddy_visc4

      !
      ! height tendency: horizontal advection term -\nabla\cdot ( hu)
      !
      ! See Ringler et al. (2010) jcp paper, eqn 19, 21, and fig. 3. 
      ! for explanation of divergence operator.
      tend_h(:,:) = 0.0
      do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCells) then
               do k=1,nVertLevels
                  flux = u(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge)
                  tend_h(k,cell1) = tend_h(k,cell1) - flux
               end do 
            end if
            if (cell2 <= nCells) then
               do k=1,nVertLevels
                  flux = u(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge)
                  tend_h(k,cell2) = tend_h(k,cell2) + flux
               end do 
            end if
      end do 
      do iCell=1,nCells
         do k=1,nVertLevels
            tend_h(k,iCell) = tend_h(k,iCell) / areaCell(iCell)
         end do
      end do

      !
      ! height tendency: vertical advection term -d/dz(hw)
      !
      if (config_vert_grid_type.eq.'zlevel') then

        do iCell=1,nCells

           tend_h(1,iCell) =   tend_h(1,iCell) + wTop(2,iCell)

           ! This next loop is to verify that h for levels below the first
           ! remain constant.  At a later time this could be replaced
           ! by just tend_h(2:nVertLevels,:) = 0.0, and then there is 
           ! no need to compute the horizontal tend_h term for k=2:nVertLevels
           ! on a zlevel grid, above.
           do k=2,nVertLevels
              tend_h(k,iCell) =   tend_h(k,iCell) &
               - (wTop(k,iCell) - wTop(k+1,iCell))
            end do

        end do
      endif ! coordinate type

      !
      ! velocity tendency: start accumulating tendency terms
      !
      tend_u(:,:) = 0.0

      !
      ! velocity tendency: vertical advection term -w du/dz
      !
      allocate(w_dudzTopEdge(nVertLevels+1))
      w_dudzTopEdge(1) = 0.0
      w_dudzTopEdge(nVertLevels+1) = 0.0
      if (config_vert_grid_type.eq.'zlevel') then
       do iEdge=1,grid % nEdgesSolve
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k=2,nVertLevels
           ! Average w from cell center to edge
           wTopEdge = 0.5*(wTop(k,cell1)+wTop(k,cell2))

           ! compute dudz at vertical interface with first order derivative.
           w_dudzTopEdge(k) = wTopEdge * (u(k-1,iEdge)-u(k,iEdge)) &
                        / (zMidZLevel(k-1) - zMidZLevel(k))
         end do

         ! Average w*du/dz from vertical interface to vertical middle of cell
         do k=1,nVertLevels
            tend_u(k,iEdge) = - 0.5 * (w_dudzTopEdge(k) + w_dudzTopEdge(k+1))
         enddo
       enddo
      endif
      deallocate(w_dudzTopEdge)

      !
      ! velocity tendency: pressure gradient
      !
      rho0Inv = 1.0/config_rho0
      if (config_vert_grid_type.eq.'isopycnal') then
        do iEdge=1,grid % nEdgesSolve
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)
          do k=1,nVertLevels
             tend_u(k,iEdge) = tend_u(k,iEdge)     &
               - (MontPot(k,cell2) - MontPot(k,cell1))/dcEdge(iEdge)
           end do
        enddo
      elseif (config_vert_grid_type.eq.'zlevel') then
        do iEdge=1,grid % nEdgesSolve
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)
          do k=1,nVertLevels
             tend_u(k,iEdge) = tend_u(k,iEdge)     &
               - rho0Inv*(  pZLevel(k,cell2) &
                          - pZLevel(k,cell1) )/dcEdge(iEdge)
          end do
        enddo
      endif

      !
      ! velocity tendency: del2 dissipation, \nu_2 \nabla^2 u
      !   computed as \nu( \nabla divergence + k \times \nabla vorticity )
      !   strictly only valid for h_mom_eddy_visc2 == constant
      !
      if ( h_mom_eddy_visc2 > 0.0 ) then
         do iEdge=1,grid % nEdgesSolve
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)

            do k=1,nVertLevels

               ! Here -( vorticity(k,vertex2) - vorticity(k,vertex1) ) / dvEdge(iEdge)
               ! is - \nabla vorticity pointing from vertex 2 to vertex 1, or equivalently 
               !    + k \times \nabla vorticity pointing from cell1 to cell2.

               u_diffusion = ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
                            -( vorticity(k,vertex2) - vorticity(k,vertex1) ) / dvEdge(iEdge)
               u_diffusion = h_mom_eddy_visc2 * u_diffusion

               tend_u(k,iEdge) = tend_u(k,iEdge) + u_diffusion

            end do
         end do
      end if

      !
      ! velocity tendency: del4 dissipation, -\nu_4 \nabla^4 u
      !   computed as \nabla^2 u = \nabla divergence + k \times \nabla vorticity
      !   applied recursively.
      !   strictly only valid for h_mom_eddy_visc4 == constant
      !
      if ( h_mom_eddy_visc4 > 0.0 ) then

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

            do k=1,nVertLevels

               u_diffusion = ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
                            -( vorticity(k,vertex2) - vorticity(k,vertex1)) / dvEdge(iEdge)
 
               delsq_u(k,iEdge) = delsq_u(k,iEdge) + u_diffusion

            end do
         end do

         ! vorticity using \nabla^2 u
         delsq_circulation(:,:) = 0.0
         do iEdge=1,nEdges
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)
            do k=1,nVertLevels
               delsq_circulation(k,vertex1) = delsq_circulation(k,vertex1) &
                  - dcEdge(iEdge) * delsq_u(k,iEdge)
               delsq_circulation(k,vertex2) = delsq_circulation(k,vertex2) &
                  + dcEdge(iEdge) * delsq_u(k,iEdge)
            end do
         end do
         do iVertex=1,nVertices
            r = 1.0 / areaTriangle(iVertex)
            do k=1,nVertLevels
               delsq_vorticity(k,iVertex) = delsq_circulation(k,iVertex) * r
            end do
         end do

         ! Divergence using \nabla^2 u
         delsq_divergence(:,:) = 0.0
         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1,nVertLevels
              delsq_divergence(k,cell1) = delsq_divergence(k,cell1) &
                + delsq_u(k,iEdge)*dvEdge(iEdge)
              delsq_divergence(k,cell2) = delsq_divergence(k,cell2) &
                - delsq_u(k,iEdge)*dvEdge(iEdge)
            end do
         end do
         do iCell = 1,nCells
            r = 1.0 / areaCell(iCell)
            do k = 1,nVertLevels
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

            do k=1,nVertLevels

               u_diffusion = (  delsq_divergence(k,cell2) &
                              - delsq_divergence(k,cell1) ) / dcEdge(iEdge)  &
                            -(  delsq_vorticity(k,vertex2) &
                              - delsq_vorticity(k,vertex1) ) / dvEdge(iEdge)
 
               tend_u(k,iEdge) = tend_u(k,iEdge) - h_mom_eddy_visc4 * u_diffusion
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

         do k=1,nVertLevels

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
      do iEdge=1,grid % nEdgesSolve

         ! forcing in top layer only
         tend_u(1,iEdge) =  tend_u(1,iEdge) &
           + u_src(1,iEdge)/rho_ref/h_edge(1,iEdge)

         ! bottom drag is the same as POP:
         ! -c |u| u  where c is unitless and 1.0e-3.
         ! see POP Reference guide, section 3.4.4.
         tend_u(nVertLevels,iEdge) = tend_u(nVertLevels,iEdge)  &
           - 1.0e-3*u(nVertLevels,iEdge) &
             *sqrt(2.0*ke_edge(nVertLevels,iEdge))/h_edge(nVertLevels,iEdge)

         ! old bottom drag, just linear friction
         ! du/dt = u/tau where tau=100 days.
         !tend_u(nVertLevels,iEdge) = tend_u(nVertLevels,iEdge)  &
         !  - u(nVertLevels,iEdge)/(100.0*86400.0)

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

      allocate(fluxVertTop(1:nVertLevels+1))
      fluxVertTop(1) = 0.0
      fluxVertTop(nVertLevels+1) = 0.0
      do iEdge=1,grid % nEdgesSolve
         do k=2,nVertLevels
           fluxVertTop(k) = vertViscTop(k) &
              * ( u(k-1,iEdge) - u(k,iEdge) ) &
              / (zMidEdge(k-1,iEdge) - zMidEdge(k,iEdge))
         enddo
         do k=1,nVertLevels
           tend_u(k,iEdge) = tend_u(k,iEdge) &
             + (fluxVertTop(k) - fluxVertTop(k+1)) &
              /(zTopEdge(k,iEdge) - zTopEdge(k+1,iEdge))
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
        nEdges, nCells, nVertLevels, num_tracers
      real (kind=RKIND) :: invAreaCell1, invAreaCell2, tracer_turb_flux
      real (kind=RKIND) :: flux, tracer_edge, r, dist
      real (kind=RKIND), dimension(:), pointer :: &
        h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: &
        u,h,wTop, h_edge, zMid, zTop
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        tracers, tend_tr
      integer, dimension(:,:), pointer :: boundaryEdge
      type (dm_info) :: dminfo

      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge
      integer, dimension(:,:), pointer :: cellsOnEdge, boundaryCell
      real (kind=RKIND), dimension(:), pointer :: zTopZLevel
      real (kind=RKIND), dimension(:,:), allocatable:: fluxVertTop, tracerTop, boundaryMask
      real (kind=RKIND), dimension(:,:,:), allocatable::tr_flux, tr_div, delsq_tracer

      real (kind=RKIND), dimension(:), allocatable:: vertDiffTop

      real (kind=RKIND) :: d2fdx2_cell1, d2fdx2_cell2
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND) :: coef_3rd_order


      u           => s % u % array
      h           => s % h % array
      boundaryCell=> grid % boundaryCell % array
      wTop        => s % wTop % array
      tracers     => s % tracers % array
      h_edge      => s % h_edge % array
      zMid        => s % zMid % array
      zTop        => s % zTop % array

      tend_tr     => tend % tracers % array
                  
      areaCell          => grid % areaCell % array
      cellsOnEdge       => grid % cellsOnEdge % array
      dvEdge            => grid % dvEdge % array
      dcEdge            => grid % dcEdge % array
      zTopZLevel        => grid % zTopZLevel % array
      boundaryEdge      => grid % boundaryEdge % array

      nEdges      = grid % nEdges
      nCells      = grid % nCells
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
      coef_3rd_order = 0.
      if (config_tracer_adv_order == 3) coef_3rd_order = 1.0
      if (config_tracer_adv_order == 3 .and. config_monotonic) coef_3rd_order = 0.25

      if (config_tracer_adv_order == 2) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= grid%nCells .and. cell2 <= grid%nCells) then
               do k=1,nVertLevels
                  do iTracer=1,num_tracers
                     tracer_edge = 0.5 * (tracers(iTracer,k,cell1) + tracers(iTracer,k,cell2))
                     flux = u(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge) * tracer_edge
                     tend_tr(iTracer,k,cell1) = tend_tr(iTracer,k,cell1) - flux/areaCell(cell1)
                     tend_tr(iTracer,k,cell2) = tend_tr(iTracer,k,cell2) + flux/areaCell(cell2)
                  end do
               end do
            end if
         end do

      else if (config_tracer_adv_order == 3) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            !-- if a cell not on the most outside ring of the halo
            if (cell1 <= grid%nCells .and. cell2 <= grid%nCells) then

               do k=1,nVertLevels

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
            end if
         end do

      else if (config_tracer_adv_order == 4) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            !-- if an edge is not on the outer-most ring of the halo
            if (cell1 <= grid%nCells .and. cell2 <= grid%nCells) then

               do k=1,nVertLevels

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
            end if
         end do

      endif   ! if (config_tracer_adv_order == 2 )


      !
      ! tracer tendency: vertical advection term -d/dz( h \phi w)
      !
      allocate(tracerTop(num_tracers,nVertLevels+1))
      tracerTop(:,1)=0.0
      tracerTop(:,nVertLevels+1)=0.0
      do iCell=1,grid % nCellsSolve 

         if (config_vert_tracer_adv.eq.'centered') then
           do k=2,nVertLevels
             do iTracer=1,num_tracers
               tracerTop(iTracer,k) = ( tracers(iTracer,k-1,iCell) &
                                       +tracers(iTracer,k  ,iCell))/2.0
             end do
           end do
         
         elseif (config_vert_tracer_adv.eq.'upwind') then
           do k=2,nVertLevels
             if (wTop(k,iCell)>0.0) then
               upwindCell = k
             else
               upwindCell = k-1
             endif
             do iTracer=1,num_tracers
               tracerTop(iTracer,k) = tracers(iTracer,upwindCell,iCell)
             end do
           end do

         endif

         do k=1,nVertLevels  
            do iTracer=1,num_tracers
               tend_tr(iTracer,k,iCell) = tend_tr(iTracer,k,iCell) &
                  - (   wTop(k  ,iCell)*tracerTop(iTracer,k  ) &
                      - wTop(k+1,iCell)*tracerTop(iTracer,k+1))
            end do
         end do

      enddo ! iCell
      deallocate(tracerTop)

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

            do k=1,grid % nVertLevels
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

            do k=1,grid % nVertLevels
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

         do iCell = 1, nCells
            r = 1.0 / areaCell(iCell)
            do k=1,nVertLevels
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

            do k=1,grid % nVertLevels
            do iTracer=1,num_tracers
               tracer_turb_flux = config_h_tracer_eddy_diff4 &
                  *(  delsq_tracer(iTracer,k,cell2)  &
                    - delsq_tracer(iTracer,k,cell1))/dcEdge(iEdge)
               flux = dvEdge (iEdge) * tracer_turb_flux

               tend_tr(iTracer,k,cell1) = tend_tr(iTracer,k,cell1) & 
                  - flux * invAreaCell1 * boundaryMask(k,iEdge)
               tend_tr(iTracer,k,cell2) = tend_tr(iTracer,k,cell2) &
                  + flux * invAreaCell2 * boundaryMask(k,iEdge)

            end do
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
      fluxVertTop(:,nVertLevels+1) = 0.0
      do iCell=1,grid % nCellsSolve 
         do k=2,nVertLevels
           do iTracer=1,num_tracers
             ! compute \kappa_v d\phi/dz
             fluxVertTop(iTracer,k) = vertDiffTop(k) &
                * (tracers(iTracer,k-1,iCell) - tracers(iTracer,k,iCell) )&
                / (zMid(k-1,iCell) -zMid(k,iCell))
           enddo
         enddo

         do k=1,nVertLevels
           dist = zTop(k,iCell) - zTop(k+1,iCell)
           do iTracer=1,num_tracers
             tend_tr(iTracer,k,iCell) = tend_tr(iTracer,k,iCell) &
               + h(k,iCell)*(fluxVertTop(iTracer,k) - fluxVertTop(iTracer,k+1))/dist
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
      real (kind=RKIND) :: flux, vorticity_abs, h_vertex, workpv, delta_p

      integer :: nCells, nEdges, nVertices, nVertLevels


      real (kind=RKIND), dimension(:), pointer :: &
        h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle, &
        hZLevel, ssh
      real (kind=RKIND), dimension(:,:), pointer :: &
        weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, w, wTop, &
        circulation, vorticity, ke, ke_edge, &
        pv_edge, pv_vertex, pv_cell, gradPVn, gradPVt, divergence, &
        zMid, zTop, zMidEdge, zTopEdge, p, pTop, MontPot, rho, temperature, salinity, pZLevel
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      real (kind=RKIND), dimension(:,:), allocatable:: div_u
      character :: c1*6

      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex, boundaryEdge, boundaryCell
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge
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
      zMid        => s % zMid % array
      zTop        => s % zTop % array
      zMidEdge    => s % zMidEdge % array
      zTopEdge    => s % zTopEdge % array
      p           => s % p % array
      pZLevel     => s % pZLevel % array
      pTop        => s % pTop % array
      MontPot     => s % MontPot % array
      ssh         => s % ssh  % array

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
                  
      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels

      boundaryEdge => grid % boundaryEdge % array
      boundaryCell => grid % boundaryCell % array

      !
      ! Find those cells that have an edge on the boundary
      !
      boundaryCell(:,:) = 0
      do iEdge=1,nEdges
       do k=1,nVertLevels
         if(boundaryEdge(k,iEdge).eq.1) then
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           boundaryCell(k,cell1) = 1
           boundaryCell(k,cell2) = 1
         endif
       enddo
      enddo


      !
      ! Compute height on cell edges at velocity locations
      !   Namelist options control the order of accuracy of the reconstructed h_edge value
      !

      coef_3rd_order = 0.
      if (config_thickness_adv_order == 3) coef_3rd_order = 1.0
      if (config_thickness_adv_order == 3 .and. config_monotonic) coef_3rd_order = 0.25

      if (config_thickness_adv_order == 2) then

         do iEdge=1,grid % nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= grid%nCells .and. cell2 <= grid%nCells) then
               do k=1,grid % nVertLevels
                  h_edge(k,iEdge) = 0.5 * (h(k,cell1) + h(k,cell2))
               end do
            end if
         end do

      else if (config_thickness_adv_order == 3) then

         do iEdge=1,grid%nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            !-- if a cell not on the most outside ring of the halo
            if (cell1 <= grid%nCells .and. cell2 <= grid%nCells) then

               do k=1,grid % nVertLevels

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
            end if      ! if (cell1 <=
         end do         ! do iEdge

      else  if (config_thickness_adv_order == 4) then

         do iEdge=1,grid%nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            !-- if a cell not on the most outside ring of the halo
            if (cell1 <= grid%nCells .and. cell2 <= grid%nCells) then

               do k=1,grid % nVertLevels

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
            end if      ! if (cell1 <=
         end do         ! do iEdge

      endif   ! if(config_thickness_adv_order == 2)

      !
      ! set the velocity in the nEdges+1 slot to zero, this is a dummy address
      !    used to when reading for edges that do not exist
      !
      u(:,nEdges+1) = 0.0

      !
      ! Compute circulation and relative vorticity at each vertex
      !
      circulation(:,:) = 0.0
      do iEdge=1,nEdges
         if (verticesOnEdge(1,iEdge) <= nVertices) then
            do k=1,nVertLevels
               circulation(k,verticesOnEdge(1,iEdge)) = circulation(k,verticesOnEdge(1,iEdge)) - dcEdge(iEdge) * u(k,iEdge)
            end do
         end if
         if (verticesOnEdge(2,iEdge) <= nVertices) then
            do k=1,nVertLevels
               circulation(k,verticesOnEdge(2,iEdge)) = circulation(k,verticesOnEdge(2,iEdge)) + dcEdge(iEdge) * u(k,iEdge)
            end do
         end if
      end do
      do iVertex=1,nVertices
         do k=1,nVertLevels
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
         if (cell1 <= nCells) then
            do k=1,nVertLevels
              divergence(k,cell1) = divergence(k,cell1) + u(k,iEdge)*dvEdge(iEdge)
            enddo
         endif
         if(cell2 <= nCells) then
            do k=1,nVertLevels
              divergence(k,cell2) = divergence(k,cell2) - u(k,iEdge)*dvEdge(iEdge)
            enddo
         end if
      end do
      do iCell = 1,nCells
        r = 1.0 / areaCell(iCell)
        do k = 1,nVertLevels
           divergence(k,iCell) = divergence(k,iCell) * r
        enddo
      enddo

      !
      ! Compute kinetic energy in each cell
      !
      ke(:,:) = 0.0
      do iCell=1,nCells
         do i=1,nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            do k=1,nVertLevels
               ke(k,iCell) = ke(k,iCell) + 0.25 * dcEdge(iEdge) * dvEdge(iEdge) * u(k,iEdge)**2.0
            end do
         end do
         do k=1,nVertLevels
            ke(k,iCell) = ke(k,iCell) / areaCell(iCell)
         end do
      end do

      !
      !
      ! Compute v (tangential) velocities
      !
      v(:,:) = 0.0
      do iEdge = 1,nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            if (eoe <= nEdges) then
               do k = 1,nVertLevels
                 v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
              end do
            end if
         end do
      end do

      !
      ! Compute ke on cell edges at velocity locations for quadratic bottom drag. 
      !
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         if (cell1 <= nCells .and. cell2 <= nCells) then
            do k=1,nVertLevels
               ke_edge(k,iEdge) = 0.5 * (ke(k,cell1) + ke(k,cell2))
            end do
         else
            do k=1,nVertLevels
               ke_edge(k,iEdge) = 0.0
            end do
         end if
      end do

      !
      ! Compute height at vertices, pv at vertices, and average pv to edge locations
      !  ( this computes pv_vertex at all vertices bounding real cells and distance-1 ghost cells )
      !
      VTX_LOOP: do iVertex = 1,nVertices
         do i=1,grid % vertexDegree
            if (cellsOnVertex(i,iVertex) > nCells) cycle VTX_LOOP
         end do
         do k=1,nVertLevels
            h_vertex = 0.0
            do i=1,grid % vertexDegree
               h_vertex = h_vertex + h(k,cellsOnVertex(i,iVertex)) * kiteAreasOnVertex(i,iVertex)
            end do
            h_vertex = h_vertex / areaTriangle(iVertex)

            pv_vertex(k,iVertex) = (fVertex(iVertex) + vorticity(k,iVertex)) / h_vertex
         end do
      end do VTX_LOOP


      !
      ! Compute gradient of PV in the tangent direction
      !   ( this computes gradPVt at all edges bounding real cells and distance-1 ghost cells )
      !
      do iEdge = 1,nEdges
         do k = 1,nVertLevels
           gradPVt(k,iEdge) = (pv_vertex(k,verticesOnEdge(2,iEdge)) - pv_vertex(k,verticesOnEdge(1,iEdge))) / &
                              dvEdge(iEdge)
         enddo
      enddo

      !
      ! Compute pv at the edges
      !   ( this computes pv_edge at all edges bounding real cells and distance-1 ghost cells )
      !
      pv_edge(:,:) = 0.0
      do iVertex = 1,nVertices
        do i=1,grid % vertexDegree
          iEdge = edgesOnVertex(i,iVertex)
          if(iEdge <= nEdges) then
            do k=1,nVertLevels
              pv_edge(k,iEdge) =  pv_edge(k,iEdge)  + 0.5 * pv_vertex(k,iVertex)
            enddo
          endif
        end do
      end do

      !
      ! Modify PV edge with upstream bias. 
      !
      do iEdge = 1,nEdges
         do k = 1,nVertLevels
           pv_edge(k,iEdge) = pv_edge(k,iEdge) - config_apvm_upwinding * v(k,iEdge) * dt * gradPVt(k,iEdge)
         enddo
      enddo


      !
      ! Compute pv at cell centers
      !    ( this computes pv_cell for all real cells and distance-1 ghost cells )
      !
      pv_cell(:,:) = 0.0
      do iVertex = 1, nVertices
       do i=1,grid % vertexDegree
         iCell = cellsOnVertex(i,iVertex)
         if (iCell <= nCells) then
           do k = 1,nVertLevels
             pv_cell(k,iCell) = pv_cell(k,iCell) + kiteAreasOnVertex(i, iVertex) * pv_vertex(k, iVertex) / areaCell(iCell)
           enddo
         endif
       enddo
      enddo

      !
      ! Compute gradient of PV in normal direction
      !   ( this computes gradPVn for all edges bounding real cells )
      !
      gradPVn(:,:) = 0.0
      do iEdge = 1,nEdges
        if( cellsOnEdge(1,iEdge) <= nCells .and. cellsOnEdge(2,iEdge) <= nCells) then
          do k = 1,nVertLevels
            gradPVn(k,iEdge) = (pv_cell(k,cellsOnEdge(2,iEdge)) - pv_cell(k,cellsOnEdge(1,iEdge))) / &
                                 dcEdge(iEdge)
          enddo
        endif
      enddo


      ! Modify PV edge with upstream bias.
      !
      do iEdge = 1,nEdges
         do k = 1,nVertLevels
           pv_edge(k,iEdge) = pv_edge(k,iEdge) - config_apvm_upwinding * u(k,iEdge) * dt * gradPVn(k,iEdge)
         enddo
      enddo

      !
      ! Compute sea surface height
      !
      do iCell=1,nCells
        ssh(iCell) = h(1,iCell) - hZLevel(1)
      enddo

      !
      ! equation of state
      !
      ! For a isopycnal model, density should remain constant.
      if (config_vert_grid_type.eq.'zlevel') then
        do iCell=1,nCells
          do k=1,nVertLevels
            ! Linear equation of state, for the time being
            rho(k,iCell) = 1000.0*(  1.0 &
               - 2.5e-4*tracers(s % index_temperature,k,iCell) &
               + 7.6e-4*tracers(s % index_salinity,k,iCell))
          end do
        end do
      endif


      ! For Isopycnal model.
      ! Compute mid- and bottom-depth of each layer, from bottom up.
      ! Then compute mid- and bottom-pressure of each layer, and 
      ! Montgomery Potential, from top down
      !
      do iCell=1,nCells

         ! h_s is ocean topography: elevation above lowest point, 
         ! and is positive. z coordinates are pos upward, and z=0
         ! is at lowest ocean point.
         zTop(nVertLevels+1,iCell) = h_s(iCell) 
         do k=nVertLevels,1,-1
            zMid(k,iCell) = zTop(k+1,iCell) + 0.5*h(k,iCell)
            zTop(k,iCell) = zTop(k+1,iCell) +     h(k,iCell)
         end do

         ! assume atmospheric pressure at the surface is zero for now.
         pTop(1,iCell) = 0.0
         do k=1,nVertLevels
            delta_p = rho(k,iCell)*gravity* h(k,iCell)
            p(k  ,iCell) = pTop(k,iCell) + 0.5*delta_p
            pTop(k+1,iCell) = pTop(k,iCell) + delta_p
         end do

         MontPot(1,iCell) = gravity * zTop(1,iCell) 
         do k=2,nVertLevels
            ! from delta M = p delta / rho
            MontPot(k,iCell) = MontPot(k-1,iCell) &
               + pTop(k,iCell)*(1.0/rho(k,iCell) - 1.0/rho(k-1,iCell)) 
         end do

      end do

      do iEdge=1,nEdges
       cell1 = cellsOnEdge(1,iEdge)
       cell2 = cellsOnEdge(2,iEdge)
       if(cell1<=nCells .and. cell2<=nCells) then
         do k=1,nVertLevels
           zTopEdge(k,iEdge) = (zTop(k,cell1)+zTop(k,cell2))/2.0
           zMidEdge(k,iEdge) = (zMid(k,cell1)+zMid(k,cell2))/2.0
         enddo
         zTopEdge(nVertLevels+1,iEdge) = ( zTop(nVertLevels+1,cell1) &
                                         + zTop(nVertLevels+1,cell2))/2.0
        endif
      enddo


      ! For z-level model.
      ! Compute pressure at middle of each level.  
      ! pZLevel and p should only differ at k=1, where p is 
      ! pressure at middle of layer including SSH, and pZLevel is
      ! pressure at a depth of hZLevel(1)/2.
      !
      do iCell=1,nCells
         ! compute pressure for z-level coordinates
         ! assume atmospheric pressure at the surface is zero for now.
         pZLevel(1,iCell) = rho(1,iCell)*gravity &
            * (h(1,iCell)-0.5*hZLevel(1)) 
         do k=2,nVertLevels
            pZLevel(k,iCell) = pZLevel(k-1,iCell)  &
              + 0.5*gravity*(  rho(k-1,iCell)*hZLevel(k-1) &
                             + rho(k  ,iCell)*hZLevel(k  ))
         end do

      end do

      ! compute vertical velocity
      if (config_vert_grid_type.eq.'isopycnal') then
        ! set vertical velocity to zero in isopycnal case
        wTop=0.0  

      elseif (config_vert_grid_type.eq.'zlevel') then

        !
        ! Compute div(u) for each cell
        ! See Ringler et al. (2010) jcp paper, eqn 19, 21, and fig. 3.
        !
        allocate(div_u(nVertLevels,nCells))
        div_u(:,:) = 0.0
        do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCells) then
               do k=1,nVertLevels
                  flux = u(k,iEdge) * dvEdge(iEdge) 
                  div_u(k,cell1) = div_u(k,cell1) + flux
               end do 
            end if
            if (cell2 <= nCells) then
               do k=1,nVertLevels
                  flux = u(k,iEdge) * dvEdge(iEdge) 
                  div_u(k,cell2) = div_u(k,cell2) - flux
               end do 
            end if
        end do 

        do iCell=1,nCells
           do k=1,nVertLevels
              div_u(k,iCell) = div_u(k,iCell) / areaCell(iCell)
           end do

           ! Vertical velocity at bottom is zero.
           ! this next line can be set permanently somewhere upon startup.
           wTop(nVertLevels+1,iCell) = 0.0
           do k=nVertLevels,1,-1
              wTop(k,iCell) = wTop(k+1,iCell) - div_u(k,iCell)*h(k,iCell)
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

end module time_integration
