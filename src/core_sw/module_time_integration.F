module time_integration

   use vector_reconstruction
   use grid_types
   use configure
   use constants
   use dmpar


   contains


   subroutine timestep(domain, dt)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Advance model state forward in time by the specified time step
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
   !                 plus grid meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt

      type (block_type), pointer :: block

      if (trim(config_time_integration) == 'RK4') then
         call rk4(domain, dt)
      else
         write(0,*) 'Unknown time integration option '//trim(config_time_integration)
         write(0,*) 'Currently, only ''RK4'' is supported.'
         stop
      end if

      block => domain % blocklist
      do while (associated(block))
         block % state % time_levs(2) % state % xtime % scalar = block % state % time_levs(1) % state % xtime % scalar + dt
         block => block % next
      end do

   end subroutine timestep


   subroutine rk4(domain, dt)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Advance model state forward in time by the specified time step using 
   !   4th order Runge-Kutta
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
   !                 plus grid meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt

      integer :: iCell, k
      type (block_type), pointer :: block
      type (state_type) :: provis

      integer :: rk_step

      real (kind=RKIND), dimension(4) :: rk_weights, rk_substep_weights

      block => domain % blocklist
      call allocate_state(provis, &
                          block % mesh % nCells, block % mesh % nEdges, block % mesh % maxEdges, block % mesh % maxEdges2, &
                          block % mesh % nVertices, block % mesh % vertexDegree, block % mesh % nVertLevels, &
                          block % mesh % nTracers)

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


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! BEGIN RK loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      do rk_step = 1, 4

! ---  update halos for diagnostic variables

        block => domain % blocklist
        do while (associated(block))
           call dmpar_exch_halo_field2dReal(domain % dminfo, provis % pv_edge % array(:,:), &
                                            block % mesh % nVertLevels, block % mesh % nEdges, &
                                            block % parinfo % edgesToSend, block % parinfo % edgesToRecv)
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
                                            block % mesh % nTracers, block % mesh % nVertLevels, block % mesh % nCells, &
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
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! END RK loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


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
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Compute height and normal wind tendencies, as well as diagnostic variables
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed tendencies for prognostic variables
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (tend_type), intent(inout) :: tend
      type (state_type), intent(in) :: s
      type (mesh_type), intent(in) :: grid

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j
      real (kind=RKIND) :: flux, vorticity_abs, workpv, q, upstream_bias

      integer :: nCells, nEdges, nVertices, nVertLevels
      real (kind=RKIND), dimension(:), pointer :: h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: vh, weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, tend_h, tend_u, &
                                                    circulation, vorticity, ke, pv_edge, divergence, h_vertex
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge
      real (kind=RKIND) :: u_diffusion, visc


      visc = config_visc

      h           => s % h % array
      u           => s % u % array
      v           => s % v % array
      h_edge      => s % h_edge % array
      circulation => s % circulation % array
      vorticity   => s % vorticity % array
      divergence  => s % divergence % array
      ke          => s % ke % array
      pv_edge     => s % pv_edge % array
      vh          => s % vh % array

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

      tend_h      => tend % h % array
      tend_u      => tend % u % array
                  
      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels


      !
      ! Compute height tendency for each cell
      !
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
      do iCell=1,grid % nCellsSolve
         do k=1,nVertLevels
            tend_h(k,iCell) = tend_h(k,iCell) / areaCell(iCell)
         end do
      end do

#ifdef LANL_FORMULATION
      !
      ! Compute u (normal) velocity tendency for each edge (cell face)
      !
      tend_u(:,:) = 0.0
      do iEdge=1,grid % nEdgesSolve
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)
         
         do k=1,nVertLevels

         !
         ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
         !                    only valid for visc == constant
         !
            u_diffusion =   ( divergence(k,cell2)  -  divergence(k,cell1) ) / dcEdge(iEdge) &
                           -(vorticity(k,vertex2)  - vorticity(k,vertex1) ) / dvEdge(iEdge)
            u_diffusion = visc * u_diffusion

            q = 0.0
            do j = 1,nEdgesOnEdge(iEdge)
               eoe = edgesOnEdge(j,iEdge)
               workpv = 0.5 * (pv_edge(k,iEdge) + pv_edge(k,eoe))
               q = q + weightsOnEdge(j,iEdge) * u(k,eoe) * workpv * h_edge(k,eoe) 
            end do

            tend_u(k,iEdge) =       &
                              q     &
                              + u_diffusion &
                              - (   ke(k,cell2) - ke(k,cell1) + &
                                    gravity * (h(k,cell2) + h_s(cell2) - h(k,cell1) - h_s(cell1)) &
                                  ) / dcEdge(iEdge)
         end do
      end do
#endif

#ifdef NCAR_FORMULATION
      !
      ! Compute u (normal) velocity tendency for each edge (cell face)
      !
      tend_u(:,:) = 0.0
      do iEdge=1,grid % nEdgesSolve
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k=1,nVertLevels
            vorticity_abs = fEdge(iEdge) + (circulation(k,vertex1) + circulation(k,vertex2)) / &
                                           (areaTriangle(vertex1) + areaTriangle(vertex2))

            workpv = 2.0 * vorticity_abs / (h(k,cell1) + h(k,cell2))

            tend_u(k,iEdge) = workpv * vh(k,iEdge) - &
                              (ke(k,cell2) - ke(k,cell1) + &
                                 gravity * (h(k,cell2) + h_s(cell2) - h(k,cell1) - h_s(cell1)) &
                              ) / &
                              dcEdge(iEdge)
         end do
      end do
#endif

   end subroutine compute_tend


   subroutine compute_scalar_tend(tend, s, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed scalar tendencies
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (tend_type), intent(inout) :: tend
      type (state_type), intent(in) :: s
      type (mesh_type), intent(in) :: grid

      integer :: iCell, iEdge, k, iTracer, cell1, cell2
      real (kind=RKIND) :: flux, tracer_edge

      tend % tracers % array(:,:,:) = 0.0
      do iEdge=1,grid % nEdges
            cell1 = grid % cellsOnEdge % array(1,iEdge)
            cell2 = grid % cellsOnEdge % array(2,iEdge)
            if (cell1 <= grid%nCells .and. cell2 <= grid%nCells) then
               do k=1,grid % nVertLevels
                  do iTracer=1,grid % nTracers
                     tracer_edge = 0.5 * (s % tracers % array(iTracer,k,cell1) + s % tracers % array(iTracer,k,cell2))
                     flux = s % u % array(k,iEdge) * grid % dvEdge % array(iEdge) * s % h_edge % array(k,iEdge) * tracer_edge
                     tend % tracers % array(iTracer,k,cell1) = tend % tracers % array(iTracer,k,cell1) - flux
                     tend % tracers % array(iTracer,k,cell2) = tend % tracers % array(iTracer,k,cell2) + flux
                  end do 
               end do 
            end if
      end do 

      do iCell=1,grid % nCellsSolve
         do k=1,grid % nVertLevelsSolve
            do iTracer=1,grid % nTracers
               tend % tracers % array(iTracer,k,iCell) = tend % tracers % array(iTracer,k,iCell) / grid % areaCell % array(iCell)
            end do
         end do
      end do

   end subroutine compute_scalar_tend


   subroutine compute_solve_diagnostics(dt, s, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Compute diagnostic fields used in the tendency computations
   !
   ! Input: grid - grid metadata
   !
   ! Output: s - computed diagnostics
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      real (kind=RKIND), intent(in) :: dt
      type (state_type), intent(inout) :: s
      type (mesh_type), intent(in) :: grid


      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j, cov
      real (kind=RKIND) :: flux, vorticity_abs, workpv

      integer :: nCells, nEdges, nVertices, nVertLevels
      real (kind=RKIND), dimension(:), pointer :: h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: vh, weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, tend_h, tend_u, &
                                                    circulation, vorticity, ke, pv_edge, pv_vertex, pv_cell, gradPVn, gradPVt, divergence, &
                                                    h_vertex, vorticity_cell
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex, boundaryEdge
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge
      real (kind=RKIND) :: r, h1, h2


      h           => s % h % array
      u           => s % u % array
      v           => s % v % array
      vh          => s % vh % array
      h_edge      => s % h_edge % array
      h_vertex    => s % h_vertex % array
      tend_h      => s % h % array
      tend_u      => s % u % array
      circulation => s % circulation % array
      vorticity   => s % vorticity % array
      divergence  => s % divergence % array
      ke          => s % ke % array
      pv_edge     => s % pv_edge % array
      pv_vertex   => s % pv_vertex % array
      pv_cell     => s % pv_cell % array
      vorticity_cell => s % vorticity_cell % array
      gradPVn     => s % gradPVn % array
      gradPVt     => s % gradPVt % array

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
                  
      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels

      boundaryEdge => grid % boundaryEdge % array

      !
      ! Compute height on cell edges at velocity locations
      !
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         if (cell1 <= nCells .and. cell2 <= nCells) then
            do k=1,nVertLevels
               h_edge(k,iEdge) = 0.5 * (h(k,cell1) + h(k,cell2))
            end do
         end if
      end do

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

#ifdef NCAR_FORMULATION
      !
      ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
      !
      vh(:,:) = 0.0
      do iEdge=1,grid % nEdgesSolve
         do j=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(j,iEdge)
            do k=1,nVertLevels
               vh(k,iEdge) = vh(k,iEdge) + weightsOnEdge(j,iEdge) * u(k,eoe) * h_edge(k,eoe)
            end do
         end do
      end do
#endif


      !
      ! Compute height at vertices, pv at vertices, and average pv to edge locations
      !  ( this computes pv_vertex at all vertices bounding real cells and distance-1 ghost cells )
      !
      VTX_LOOP: do iVertex = 1,nVertices
         do i=1,grid % vertexDegree
            if (cellsOnVertex(i,iVertex) > nCells) cycle VTX_LOOP
         end do
         do k=1,nVertLevels
            h_vertex(k,iVertex) = 0.0
            do i=1,grid % vertexDegree
               h_vertex(k,iVertex) = h_vertex(k,iVertex) + h(k,cellsOnVertex(i,iVertex)) * kiteAreasOnVertex(i,iVertex)
            end do
            h_vertex(k,iVertex) = h_vertex(k,iVertex) / areaTriangle(iVertex)

            pv_vertex(k,iVertex) = (fVertex(iVertex) + vorticity(k,iVertex)) / h_vertex(k,iVertex)
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
           pv_edge(k,iEdge) = pv_edge(k,iEdge) - 0.5 * v(k,iEdge) * dt * gradPVt(k,iEdge)
         enddo
      enddo


      !
      ! Compute pv at cell centers
      !    ( this computes pv_cell for all real cells and distance-1 ghost cells )
      !
      pv_cell(:,:) = 0.0
      vorticity_cell(:,:) = 0.0
      do iVertex = 1, nVertices
       do i=1,grid % vertexDegree
         iCell = cellsOnVertex(i,iVertex)
         if (iCell <= nCells) then
           do k = 1,nVertLevels
             pv_cell(k,iCell) = pv_cell(k,iCell) + kiteAreasOnVertex(i, iVertex) * pv_vertex(k, iVertex) / areaCell(iCell)
             vorticity_cell(k,iCell) = vorticity_cell(k,iCell) + kiteAreasOnVertex(i, iVertex) * vorticity(k, iVertex) / areaCell(iCell)
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
           pv_edge(k,iEdge) = pv_edge(k,iEdge) - 0.5 * u(k,iEdge) * dt * gradPVn(k,iEdge)
         enddo
      enddo

      !
      ! set pv_edge = fEdge / h_edge at boundary points
      !
   !  if (maxval(boundaryEdge).ge.0) then
   !  do iEdge = 1,nEdges
   !     cell1 = cellsOnEdge(1,iEdge)
   !     cell2 = cellsOnEdge(2,iEdge)
   !     do k = 1,nVertLevels
   !       if(boundaryEdge(k,iEdge).eq.1) then
   !         v(k,iEdge) = 0.0
   !         if(cell1.gt.0) then
   !            h1 = h(k,cell1)
   !            pv_edge(k,iEdge) = fEdge(iEdge) / h1
   !            h_edge(k,iEdge) = h1
   !         else
   !            h2 = h(k,cell2)
   !            pv_edge(k,iEdge) = fEdge(iEdge) / h2
   !            h_edge(k,iEdge) = h2
   !         endif
   !       endif
   !     enddo
   !  enddo
   !  endif


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
      tend_u               => tend % u % array

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
