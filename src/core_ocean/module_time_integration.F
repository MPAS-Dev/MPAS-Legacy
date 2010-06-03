module time_integration

   use grid_types
   use configure
   use constants
   use dmpar
   ! xsad 10-02-05:
   use vector_reconstruction
   ! xsad 10-02-05 end


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
      integer errorcode,ierr

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
        block % time_levs(2) % state % xtime % scalar = block % time_levs(1) % state % xtime % scalar + dt
   !     ! mrp 100310  I added this to avoid running with NaNs
   !     if (isNaN(sum(block % time_levs(2) % state % u % array))) then
   !        print *, 'Stopping: NaN detected'
   !        call MPI_abort(MPI_COMM_WORLD,errorcode,ierr)
   !     endif
   !     ! mrp 100310 end
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

      integer :: iCell, k
      type (block_type), pointer :: block

      integer, parameter :: PROVIS = 1
      integer, parameter :: TEND   = 2
      integer :: rk_step

      real (kind=RKIND), dimension(4) :: rk_weights, rk_substep_weights


      !
      ! Initialize time_levs(2) with state at current time
      ! Initialize first RK state
      ! Couple tracers time_levs(2) with h in time-levels
      ! Initialize RK weights
      !
      block => domain % blocklist
      do while (associated(block))

         block % time_levs(2) % state % u % array(:,:) = block % time_levs(1) % state % u % array(:,:)
         block % time_levs(2) % state % h % array(:,:) = block % time_levs(1) % state % h % array(:,:)
         do iCell=1,block % mesh % nCells  ! couple tracers to h
           do k=1,block % mesh % nVertLevels
             block % time_levs(2) % state % tracers % array(:,k,iCell) = block % time_levs(1) % state % tracers % array(:,k,iCell) &
                                                                       * block % time_levs(1) % state % h % array(k,iCell)
            end do
         end do

         call copy_state(block % time_levs(1) % state, block % intermediate_step(PROVIS))

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
           call dmpar_exch_halo_field2dReal(domain % dminfo, block % intermediate_step(PROVIS) % pv_edge % array(:,:), &
                                            block % mesh % nVertLevels, block % mesh % nEdges, &
                                            block % parinfo % edgesToSend, block % parinfo % edgesToRecv)
           block => block % next
        end do

! ---  compute tendencies

        block => domain % blocklist
        do while (associated(block))

           call compute_tend(block % intermediate_step(TEND), block % intermediate_step(PROVIS), block % mesh)
           call compute_scalar_tend(block % intermediate_step(TEND), block % intermediate_step(PROVIS), block % mesh)
           call enforce_boundaryEdge(block % intermediate_step(TEND), block % mesh)
           block => block % next
        end do

! ---  update halos for prognostic variables

        block => domain % blocklist
        do while (associated(block))
           call dmpar_exch_halo_field2dReal(domain % dminfo, block % intermediate_step(TEND) % u % array(:,:), &
                                            block % mesh % nVertLevels, block % mesh % nEdges, &
                                            block % parinfo % edgesToSend, block % parinfo % edgesToRecv)
           call dmpar_exch_halo_field2dReal(domain % dminfo, block % intermediate_step(TEND) % h % array(:,:), &
                                            block % mesh % nVertLevels, block % mesh % nCells, &
                                            block % parinfo % cellsToSend, block % parinfo % cellsToRecv)
           call dmpar_exch_halo_field3dReal(domain % dminfo, block % intermediate_step(TEND) % tracers % array(:,:,:), &
                                            block % mesh % nTracers, block % mesh % nVertLevels, block % mesh % nCells, &
                                            block % parinfo % cellsToSend, block % parinfo % cellsToRecv)
           block => block % next
        end do

! ---  compute next substep state

        if (rk_step < 4) then
           block => domain % blocklist
           do while (associated(block))
              block % intermediate_step(PROVIS) % u % array(:,:)       = block % time_levs(1) % state % u % array(:,:)  &
                                         + rk_substep_weights(rk_step) * block % intermediate_step(TEND) % u % array(:,:)
              block % intermediate_step(PROVIS) % h % array(:,:)       = block % time_levs(1) % state % h % array(:,:)  &
                                         + rk_substep_weights(rk_step) * block % intermediate_step(TEND) % h % array(:,:)
              do iCell=1,block % mesh % nCells
                 do k=1,block % mesh % nVertLevels
                    block % intermediate_step(PROVIS) % tracers % array(:,k,iCell) = ( &
                                                                      block % time_levs(1) % state % h % array(k,iCell) * &
                                                                      block % time_levs(1) % state % tracers % array(:,k,iCell)  &
                                      + rk_substep_weights(rk_step) * block % intermediate_step(TEND) % tracers % array(:,k,iCell) &
                                                                     ) / block % intermediate_step(PROVIS) % h % array(k,iCell)
                 end do
              end do
              if (config_test_case == 1) then    ! For case 1, wind field should be fixed
                 block % intermediate_step(PROVIS) % u % array(:,:) = block % time_levs(1) % state % u % array(:,:)
              end if
              call compute_solve_diagnostics(dt, block % intermediate_step(PROVIS), block % mesh)
              block => block % next
           end do
        end if

!--- accumulate update (for RK4)

        block => domain % blocklist
        do while (associated(block))
           block % time_levs(2) % state % u % array(:,:) = block % time_levs(2) % state % u % array(:,:) &
                                   + rk_weights(rk_step) * block % intermediate_step(TEND) % u % array(:,:) 
           block % time_levs(2) % state % h % array(:,:) = block % time_levs(2) % state % h % array(:,:) &
                                   + rk_weights(rk_step) * block % intermediate_step(TEND) % h % array(:,:) 
           do iCell=1,block % mesh % nCells
              do k=1,block % mesh % nVertLevels
                 block % time_levs(2) % state % tracers % array(:,k,iCell) =  &
                                                                       block % time_levs(2) % state % tracers % array(:,k,iCell) &
                                               + rk_weights(rk_step) * block % intermediate_step(TEND) % tracers % array(:,k,iCell)
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
               block % time_levs(2) % state % tracers % array(:,k,iCell) = &
                                                                     block % time_levs(2) % state % tracers % array(:,k,iCell)  &
                                                                   / block % time_levs(2) % state % h % array(k,iCell)
            end do
         end do

         if (config_test_case == 1) then    ! For case 1, wind field should be fixed
            block % time_levs(2) % state % u % array(:,:) = block % time_levs(1) % state % u % array(:,:)
         end if

         call compute_solve_diagnostics(dt, block % time_levs(2) % state, block % mesh)

         call reconstruct(block % time_levs(2) % state, block % mesh)

         block => block % next
      end do

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

      type (grid_state), intent(inout) :: tend
      type (grid_state), intent(in) :: s
      type (grid_meta), intent(in) :: grid

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j
      real (kind=RKIND) :: flux, vorticity_abs, h_vertex, workpv, q, upstream_bias, dist
      real (kind=RKIND), allocatable, dimension(:) :: fluxVert

      integer :: nCells, nEdges, nVertices, nVertLevels
      real (kind=RKIND), dimension(:), pointer :: h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle, zSurfaceEdge
      real (kind=RKIND), dimension(:,:), pointer :: vh, weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, tend_h, tend_u, &
                                                    circulation, vorticity, ke, pv_edge, divergence, zBotEdge, zMidEdge
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge
      real (kind=RKIND) :: u_diffusion, visc

      !mrp 100112:
      real (kind=RKIND), dimension(:,:), pointer :: MontPot
      !mrp 100112 end

      real (kind=RKIND), dimension(:,:), pointer :: u_src
      real (kind=RKIND), parameter :: rho_ref = 1000.0

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
      MontPot     => s % MontPot % array
      zBotEdge    => s % zBotEdge % array
      zSurfaceEdge=> s % zSurfaceEdge % array
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

      tend_h      => tend % h % array
      tend_u      => tend % u % array
                  
      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels

      u_src => grid % u_src % array


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
                   - (   ke(k,cell2) - ke(k,cell1)  &
                       + MontPot(k,cell2) - MontPot(k,cell1) &
                         ) / dcEdge(iEdge)
         end do
      end do

     do iEdge=1,grid % nEdgesSolve
      ! surface forcing
        tend_u(1,iEdge) =  tend_u(1,iEdge) + u_src(1,iEdge)/rho_ref/h_edge(1,iEdge)

      ! bottom drag
        tend_u(nVertLevels,iEdge) = tend_u(nVertLevels,iEdge) - u(nVertLevels,iEdge)/(100.0*86400.0)

     enddo

! vertical mixing
      allocate(fluxVert(0:nVertLevels))
      do iEdge=1,grid % nEdgesSolve
        fluxVert(0) = 0.0
        fluxVert(nVertLevels) = 0.0
        do k=nVertLevels-1,1,-1
          fluxVert(k) = ( u(k,iEdge) - u(k+1,iEdge) ) / (zMidEdge(k,iEdge) - zMidEdge(k+1,iEdge))
        enddo
        fluxVert = 1.0e-4 * fluxVert
        do k=1,nVertLevels
          if(k.eq.1) then
             dist = zSurfaceEdge(iEdge) - zBotEdge(k,iEdge)
          else
             dist = zBotEdge(k-1,iEdge) - zBotEdge(k,iEdge)
          endif
          tend_u(k,iEdge) = tend_u(k,iEdge) + (fluxVert(k-1) - fluxVert(k))/dist
        enddo
     enddo
     deallocate(fluxVert)

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
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed scalar tendencies
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_state), intent(inout) :: tend
      type (grid_state), intent(in) :: s
      type (grid_meta), intent(in) :: grid

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
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields used in the tendency computations
   !
   ! Input: grid - grid metadata
   !
   ! Output: s - computed diagnostics
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), intent(in) :: dt
      type (grid_state), intent(inout) :: s
      type (grid_meta), intent(in) :: grid


      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j, cov
      real (kind=RKIND) :: flux, vorticity_abs, h_vertex, workpv

      integer :: nCells, nEdges, nVertices, nVertLevels
      real (kind=RKIND), dimension(:), pointer :: h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: vh, weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, tend_h, tend_u, &
                                                    circulation, vorticity, ke, pv_edge, pv_vertex, pv_cell, gradPVn, gradPVt, divergence
      real (kind=RKIND), dimension(:,:), pointer :: &
        zMid, zBot, pmid, pbot, MontPot, rho, zBotEdge, zMidEdge
      real (kind=RKIND), dimension(:), pointer :: zSurface, zSurfaceEdge
      real (kind=RKIND) :: delta_p
      character :: c1*6

      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex, boundaryEdge
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge
      real (kind=RKIND) :: r, h1, h2


      h           => s % h % array
      u           => s % u % array
      v           => s % v % array
      vh          => s % vh % array
      h_edge      => s % h_edge % array
      tend_h      => s % h % array
      tend_u      => s % u % array
      circulation => s % circulation % array
      vorticity   => s % vorticity % array
      divergence  => s % divergence % array
      ke          => s % ke % array
      pv_edge     => s % pv_edge % array
      pv_vertex   => s % pv_vertex % array
      pv_cell     => s % pv_cell % array
      gradPVn     => s % gradPVn % array
      gradPVt     => s % gradPVt % array
      !mrp 100112:
      rho         => s % rho % array
      zMid        => s % zMid % array
      zBot        => s % zBot % array
      zMidEdge    => s % zMidEdge % array
      zBotEdge    => s % zBotEdge % array
      zSurfaceEdge=> s % zSurfaceEdge % array
      pmid        => s % pmid % array
      pbot        => s % pbot % array
      MontPot     => s % MontPot % array
      zSurface    => s % zSurface % array

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
         elseif(cell1 <= nCells) then
            do k=1,nVertLevels
               h_edge(k,iEdge) = h(k,cell1)
            end do
         elseif(cell2 <= nCells) then
            do k=1,nVertLevels
               h_edge(k,iEdge) = h(k,cell2)
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
            if (cellsOnVertex(i,iVertex) > nVertices) cycle VTX_LOOP
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
           pv_edge(k,iEdge) = pv_edge(k,iEdge) - 0.5 * v(k,iEdge) * dt * gradPVt(k,iEdge)
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
           pv_edge(k,iEdge) = pv_edge(k,iEdge) - 0.5 * u(k,iEdge) * dt * gradPVn(k,iEdge)
         enddo
      enddo

      !
      ! Compute mid- and bottom-depth of each layer, from bottom up.
      ! Then compute mid- and bottom-pressure of each layer, and 
      ! Montgomery Potential, from top down
      !
      do iCell=1,nCells

         ! h_s is ocean topography: elevation above lowest point, 
         ! and is positive. z coordinates are pos upward, and z=0
         ! is at lowest ocean point.
         zBot(nVertLevels,iCell) = h_s(iCell) 
         zMid(nVertLevels,iCell) = zBot(nVertLevels,iCell) + 0.5*h(nVertLevels,iCell)
         do k=nVertLevels-1,1,-1
            zBot(k,iCell) = zBot(k+1,iCell) + h(k+1,iCell)
            zMid(k,iCell) = zBot(k,iCell) + 0.5*h(k,iCell)
         end do
         ! rather than using zBot(0,iCell), I am adding another 2D variable.
         zSurface(iCell) = zBot(1,iCell) + h(1,iCell)

         ! assume pressure at the surface is zero for now.
         pmid(1,iCell) = 0.5*rho(1,iCell)*gravity* h(1,iCell) ! + psurf(iCell)
         pbot(1,iCell) =     rho(1,iCell)*gravity* h(1,iCell) ! + psurf(iCell)
         MontPot(1,iCell) = gravity * zSurface(iCell) 
         do k=2,nVertLevels
            delta_p = rho(k,iCell)*gravity* h(k,iCell)
            pmid(k,iCell) = pbot(k-1,iCell) + 0.5*delta_p
            pbot(k,iCell) = pbot(k-1,iCell) + delta_p

            ! from delta M = p delta / rho
            MontPot(k,iCell) = MontPot(k-1,iCell) &
               + pbot(k-1,iCell)*(1.0/rho(k,iCell) - 1.0/rho(k-1,iCell)) 
         end do
      end do

      do iEdge=1,nEdges
       cell1 = cellsOnEdge(1,iEdge)
       cell2 = cellsOnEdge(2,iEdge)
       if(cell1<=nCells .and. cell2<=nCells) then
        do k=1,nVertLevels
          zBotEdge(k,iEdge) = (zBot(k,cell1)+zBot(k,cell2))/2.0
          zMidEdge(k,iEdge) = (zMid(k,cell1)+zMid(k,cell2))/2.0
        enddo
        zSurfaceEdge(iEdge) = (zSurface(cell1) + zSurface(cell2))/2.0
        endif
      enddo


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

      type (grid_state), intent(inout) :: tend
      type (grid_meta), intent(in) :: grid

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
