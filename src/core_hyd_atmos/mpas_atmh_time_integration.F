module atmh_time_integration

   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_dmpar
   use mpas_vector_reconstruction


   contains


   subroutine atmh_timestep(domain, dt, timeStamp)
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
      character(len=*), intent(in) :: timeStamp

      type (block_type), pointer :: block

      if (trim(config_time_integration) == 'SRK3') then
         call atmh_srk3(domain, dt)
      else
         write(0,*) 'Unknown time integration option '//trim(config_time_integration)
         write(0,*) 'Currently, only ''SRK3'' is supported.'
         stop
      end if

      block => domain % blocklist
      do while (associated(block))
         block % state % time_levs(2) % state % xtime % scalar = timeStamp
         block => block % next
      end do

   end subroutine atmh_timestep


   subroutine atmh_srk3(domain, dt)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Advance model state forward in time by the specified time step using 
   !   time-split RK3 scheme
   !
   ! Hydrostatic (primitive eqns.) solver
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

      integer, parameter :: TEND   = 1
      integer :: rk_step, number_of_sub_steps

      real (kind=RKIND), dimension(3) :: rk_timestep, rk_sub_timestep
      integer, dimension(3) :: number_sub_steps
      integer :: small_step
      logical, parameter :: debug = .false.
      logical, parameter :: debug_mass_conservation = .true.

      real (kind=RKIND) :: domain_mass, scalar_mass, scalar_min, scalar_max
      real (kind=RKIND) :: global_domain_mass, global_scalar_mass, global_scalar_min, global_scalar_max

      !
      ! Initialize time_levs(2) with state at current time
      ! Initialize RK weights
      !

      number_of_sub_steps = config_number_of_sub_steps

      rk_timestep(1) = dt/3.
      rk_timestep(2) = dt/2.
      rk_timestep(3) = dt

      rk_sub_timestep(1) = dt/3.
      rk_sub_timestep(2) = dt/real(number_of_sub_steps)
      rk_sub_timestep(3) = dt/real(number_of_sub_steps)

      number_sub_steps(1) = 1
      number_sub_steps(2) = number_of_sub_steps/2
      number_sub_steps(3) = number_of_sub_steps

      if(debug) write(0,*) ' copy step in rk solver '

      block => domain % blocklist
      do while (associated(block))
         call mpas_copy_state( block % state % time_levs(2) % state, block % state % time_levs(1) % state )
         block => block % next
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! BEGIN RK loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


      do rk_step = 1, 3

        if(debug) write(0,*) ' rk substep ', rk_step

        call mpas_dmpar_exch_halo_field(domain % blocklist % mesh % qtot)
        call mpas_dmpar_exch_halo_field(domain % blocklist % mesh % cqu)
        call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % h)
        call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % pressure)
        call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % geopotential)
        call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % alpha)
        call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % pv_edge)
        if (config_h_mom_eddy_visc4 > 0.0) then
           call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % divergence)
           call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % vorticity)
        end if

        if(debug) write(0,*) ' rk substep ', rk_step

        block => domain % blocklist
        do while (associated(block))
           call atmh_compute_dyn_tend( block % tend, block % state % time_levs(2) % state, block % mesh )
           block => block % next
        end do

        if(debug) write(0,*) ' returned from dyn_tend '

        !
        ! ---  update halos for tendencies
        !
        call mpas_dmpar_exch_halo_field(domain % blocklist % tend % u)
        call mpas_dmpar_exch_halo_field(domain % blocklist % tend % theta)


        ! ---  advance over sub_steps

        block => domain % blocklist
        do while (associated(block))
           !
           ! Note: Scalars in new time level shouldn't be overwritten, since their provisional values 
           !    from the previous RK step are needed to compute new scalar tendencies in advance_scalars. 
           !    A cleaner way of preserving scalars should be added in future.
           !
           block % mesh % scalars_old % array(:,:,:) = block % state % time_levs(2) % state % scalars % array(:,:,:)
           call mpas_copy_state( block % state % time_levs(2) % state, block % state % time_levs(1) % state )
           block % state % time_levs(2) % state % scalars % array(:,:,:) = block % mesh % scalars_old % array(:,:,:)
           block => block % next
        end do

        if(debug) write(0,*) ' returned from copy_state '

        do small_step = 1, number_sub_steps(rk_step)

           if(debug) write(0,*) ' small step ',small_step
      
           block => domain % blocklist
           do while (associated(block))
              call atmh_advance_dynamics( block % tend, block % state % time_levs(2) % state,  &
                                     block % mesh,                                                           &
                                     small_step, number_sub_steps(rk_step), rk_sub_timestep(rk_step) )
              block => block % next
           end do

          if(debug) write(0,*) ' dynamics advance complete '
  
           !  will need communications here?
           !
           ! ---  update halos for prognostic variables
           !
!!              call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % u)
           call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % h_edge)
!!              call mpas_dmpar_exch_halo_field(domain % blocklist % mesh % uhAvg)
!!              call mpas_dmpar_exch_halo_field(domain % blocklist % mesh % wwAvg)
           call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % theta)
!!              call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % h)
!!              call mpas_dmpar_exch_halo_field(domain % blocklist % tend % h)
           call mpas_dmpar_exch_halo_field(domain % blocklist % mesh % dpsdt)
           call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % surface_pressure)
           call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % alpha)
           call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % ww)
           call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % pressure)
!!              call mpas_dmpar_exch_halo_field(domain % blocklist % mesh % pressure_old)
           call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % geopotential)

        end do

        if(debug) write(0,*) ' advance scalars '


        ! ---  advance scalars with time integrated mass fluxes

        block => domain % blocklist
        do while (associated(block))
           !
           ! Note: The advance_scalars_mono routine can be used without limiting, and thus, encompasses 
           !       the functionality of the advance_scalars routine; however, it is noticeably slower, 
           !       so we keep the advance_scalars routine as well
           !
           if (rk_step < 3 .or. (.not. config_monotonic .and. .not. config_positive_definite)) then
              call atmh_advance_scalars( block % tend,                                                               &
                                    block % state % time_levs(1) % state, block % state % time_levs(2) % state, &
                                    block % mesh, rk_timestep(rk_step) )
           else
               call atmh_advance_scalars_mono(block % tend, &
                                              block % state % time_levs(1) % state, block % state % time_levs(2) % state, &
                                              block % mesh, rk_timestep(rk_step), rk_step, 3)
           end if
           block => block % next
        end do

        call mpas_dmpar_exch_halo_field(domain % blocklist % tend % scalars)
        call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % scalars)

        
        if(debug) write(0,*) ' advance scalars complete '

        ! --- compute some diagnostic quantities for the next timestep

        block => domain % blocklist
        do while (associated(block))
           call atmh_compute_solver_constants( block % state % time_levs(2) % state, block % mesh )
           call atmh_compute_state_diagnostics( block % state % time_levs(2) % state, block % mesh )
           call atmh_compute_solve_diagnostics( dt, block % state % time_levs(2) % state, block % mesh )
           block => block % next
        end do

        if(debug) write(0,*) ' diagnostics complete '
      

        !  might need communications here *****************************

      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! END RK loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      !
      ! Compute full velocity vectors at cell centers, and compute vertical velocity diagnostic
      !
      block => domain % blocklist
      do while (associated(block))
         call mpas_reconstruct(block % mesh, block % state % time_levs(2) % state % u % array, &
                          block % diag % uReconstructX % array,                           &
                          block % diag % uReconstructY % array,                           &
                          block % diag % uReconstructZ % array,                           &
                          block % diag % uReconstructZonal % array,                       &
                          block % diag % uReconstructMeridional % array                   &
                         )

         call atmh_compute_w(block % state % time_levs(2) % state, block % state % time_levs(1) % state, block % mesh, dt)
         block => block % next
      end do

      if(debug) write(0,*) ' rk step complete - mass diagnostics '

      if(debug .or. debug_mass_conservation) then
         domain_mass = 0.
         scalar_mass = 0.
         block => domain % blocklist
         scalar_min = block % state % time_levs(2) % state % scalars % array (2,1,1)
         scalar_max = block % state % time_levs(2) % state % scalars % array (2,1,1)
         do while(associated(block))
           do iCell = 1, block % mesh % nCellsSolve
             domain_mass = domain_mass + block % state % time_levs(2) % state % surface_pressure % array (iCell) * &
                                         block % mesh % areaCell % array (iCell) &
                                       - block % state % time_levs(2) % state % pressure % array (block % mesh % nVertLevels + 1, 1) * &
                                         block % mesh % areaCell % array (iCell)
             do k=1, block % mesh % nVertLevels   ! Could be nVertLevelsSolve?
               scalar_mass = scalar_mass - block % state % time_levs(2) % state % scalars % array (2,k,iCell) * &
                                           block % state % time_levs(2) % state % h % array (k,iCell) * &
                                           block % mesh % dnw % array (k) * &
                                           block % mesh % areaCell % array (iCell)
               scalar_min = min(scalar_min,block % state % time_levs(2) % state % scalars % array (2,k,iCell))
               scalar_max = max(scalar_max,block % state % time_levs(2) % state % scalars % array (2,k,iCell))
             end do
           end do
           block => block % next
         end do
         call mpas_dmpar_sum_real(domain % dminfo, domain_mass, global_domain_mass)
         call mpas_dmpar_sum_real(domain % dminfo, scalar_mass, global_scalar_mass)
         call mpas_dmpar_min_real(domain % dminfo, scalar_min, global_scalar_min)
         call mpas_dmpar_max_real(domain % dminfo, scalar_max, global_scalar_max)
         write(0,*) ' mass in the domain = ',global_domain_mass
         write(0,*) ' scalar mass in the domain = ',global_scalar_mass
         write(0,*) ' scalar_min, scalar_max ',global_scalar_min, global_scalar_max
      end if


   end subroutine atmh_srk3

!------------------------------------------------------------------------------------------------------------------

   subroutine atmh_compute_solver_constants(s, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Compute height and normal wind tendencies, as well as diagnostic variables
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed diagnostics (parallel velocities, v; mass fluxes, vh; 
   !                circulation; vorticity; and kinetic energy, ke) and the 
   !                tendencies for height (h) and u (u)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (state_type), intent(in) :: s
      type (mesh_type), intent(inout) :: grid

      integer :: iEdge, iCell, k, cell1, cell2, iq

      integer :: nCells, nEdges, nVertLevels

      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertLevels = grid % nVertLevels

      grid % qtot % array = 0.
      grid % cqu % array = 1.

      if (s % num_scalars > 0) then

        do iCell = 1, nCells
          do k = 1, nVertLevels
            do iq = s % moist_start, s % moist_end
              grid % qtot % array(k,iCell) = grid % qtot % array(k,iCell) + s % scalars % array (iq, k, iCell)
            end do
          end do
        end do

        do iEdge = 1, nEdges
          do k = 1, nVertLevels
            cell1 = grid % cellsOnEdge % array(1,iEdge)
            cell2 = grid % cellsOnEdge % array(2,iEdge)
            grid % cqu % array(k,iEdge) = 1./( 1. + 0.5*(grid % qtot % array(k,cell1)+grid % qtot % array(k,cell2)) )
          end do
        end do

      end if

      end subroutine atmh_compute_solver_constants

!------------------------------------------------------------------------------------------------------------------

   subroutine atmh_compute_state_diagnostics(s, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Compute diagnostic variables
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed diagnostics (parallel velocities, v; mass fluxes, vh; 
   !                circulation; vorticity; and kinetic energy, ke) and the 
   !                tendencies for height (h) and u (u)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (state_type), intent(inout) :: s
      type (mesh_type), intent(inout) :: grid

      real (kind=RKIND), dimension(:,:,:), pointer :: scalar
      real (kind=RKIND), dimension(:,:), pointer :: h, pressure, qtot, alpha, geopotential, theta
      real (kind=RKIND), dimension(:,:), pointer :: theta_old, ww_old, u_old, u, ww, h_edge_old, h_edge, h_old
      real (kind=RKIND), dimension(:), pointer :: surface_pressure, dbn, dnu, dnw

      integer :: iEdge, iCell, k, cell1, cell2, iq
      integer :: nCells, nEdges, nVertLevels

      real (kind=RKIND) :: p0,tm,ptop,ptmp

      h                => s % h % array
      theta            => s % theta % array
      pressure         => s % pressure % array
      qtot             => grid % qtot % array
      surface_pressure => s % surface_pressure % array
      alpha            => s % alpha % array
      geopotential     => s % geopotential % array
      scalar           => s % scalars % array
      theta_old        => grid % theta_old % array
      u_old            => grid % u_old % array
      ww_old           => grid % ww_old % array
      h_old            => grid % h_old % array
      h_edge_old       => grid % h_edge_old % array
      h_edge           => s % h_edge % array
      u                => s % u % array
      ww               => s % ww % array

      dbn              => grid % dbn % array
      dnu              => grid % dnu % array
      dnw              => grid % dnw % array

      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertLevels = grid % nVertLevels



!      ptop        = grid % ptop
!      p0          = grid % p0
!       ptop = 219.4067
       p0 = 1e+05
       ptop = pressure(nVertLevels+1,1)

!       write(0,*) ' ptop in compute_state_diagnostics ',ptop

!*****************************

      do iCell = 1, nCells
        do k=1,nVertLevels
          h(k,iCell) = (1.-dbn(k))*(p0-ptop)+dbn(k)*(surface_pressure(iCell)-ptop)
        end do

        do k = nVertLevels, 1, -1
          pressure(k,iCell) = pressure(k+1,iCell) - dnw(k)*(1.+qtot(k,iCell))*h(k,iCell)
        end do

        do k=1, nVertLevels
          ! note that theta is not coupled here
          tm = (1.+1.61*scalar(1,k,iCell))*theta(k,iCell)  !  assume scalar(1) is qv here?
          alpha(k,iCell) = (rgas/p0)*tm*(0.5*(pressure(k+1,iCell)+pressure(k,iCell))/p0)**cvpm
          geopotential(k+1,iCell) = geopotential(k,iCell) - dnw(k)*h(k,iCell)*alpha(k,iCell)
        end do
      end do

      theta_old(:,:) = theta(:,:)
      ww_old(:,:) = ww(:,:)
      u_old(:,:) = u(:,:)
      h_edge_old(:,:) = h_edge(:,:)
      h_old(:,:) = h(:,:)

      end subroutine atmh_compute_state_diagnostics

!------------------------------------------------------------------------------------------

   subroutine atmh_compute_dyn_tend(tend, s, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Compute height and normal wind tendencies, as well as diagnostic variables
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed diagnostics (parallel velocities, v; mass fluxes, vh; 
   !                circulation; vorticity; and kinetic energy, ke) and the 
   !                tendencies for height (h) and u (u)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (tend_type), intent(inout) :: tend
      type (state_type), intent(in) :: s
      type (mesh_type), intent(in) :: grid

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j
      real (kind=RKIND) :: flux, vorticity_abs, h_vertex, workpv, q, upstream_bias

      integer :: nCells, nEdges, nVertices, nVertLevels
      real (kind=RKIND) :: h_mom_eddy_visc2,   v_mom_eddy_visc2,   h_mom_eddy_visc4
      real (kind=RKIND) :: h_theta_eddy_visc2, v_theta_eddy_visc2, h_theta_eddy_visc4
      real (kind=RKIND) :: u_diffusion
      real (kind=RKIND), dimension(:), pointer :: h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle, p_s
      real (kind=RKIND), dimension(:,:), pointer :: vh, weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, tend_h, tend_u, &
                                                    circulation, divergence, vorticity, ke, pv_edge, geopotential, theta, ww, &
                                                    h_diabatic, tend_theta
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge

      real (kind=RKIND), dimension( grid % nVertLevels + 1 ) :: wdtn, wdun
      real (kind=RKIND) :: theta_edge, theta_turb_flux, z1, z2, z3, z4, zm, z0, zp, r
      real (kind=RKIND) :: d2fdx2_cell1, d2fdx2_cell2

      real (kind=RKIND), dimension(:), pointer :: rdnu, rdnw, fnm, fnp

      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_theta, delsq_divergence
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_u
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_circulation, delsq_vorticity

      h            => s % h % array
      u            => s % u % array
      h_edge       => s % h_edge % array
      circulation  => s % circulation % array
      divergence   => s % divergence % array
      vorticity    => s % vorticity % array
      ke           => s % ke % array
      pv_edge      => s % pv_edge % array
      geopotential => s % geopotential % array
      theta        => s % theta % array

      weightsOnEdge     => grid % weightsOnEdge % array
      cellsOnEdge       => grid % cellsOnEdge % array  
      verticesOnEdge    => grid % verticesOnEdge % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array
      dcEdge            => grid % dcEdge % array
      dvEdge            => grid % dvEdge % array
      areaCell          => grid % areaCell % array
      areaTriangle      => grid % areaTriangle % array
      fEdge             => grid % fEdge % array
      deriv_two         => grid % deriv_two % array

      vh          => tend % vh % array
      tend_u      => tend % u % array
      tend_theta  => tend % theta % array
      h_diabatic  => grid % h_diabatic % array

      ww          => s % ww % array
      rdnu        => grid % rdnu % array
      rdnw        => grid % rdnw % array
      fnm         => grid % fnm % array
      fnp         => grid % fnp % array
                  
      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertLevels = grid % nVertLevels
      nVertices   = grid % nVertices

      h_mom_eddy_visc2 = config_h_mom_eddy_visc2
      h_mom_eddy_visc4 = config_h_mom_eddy_visc4
      v_mom_eddy_visc2 = config_v_mom_eddy_visc2
      h_theta_eddy_visc2 = config_h_theta_eddy_visc2
      h_theta_eddy_visc4 = config_h_theta_eddy_visc4
      v_theta_eddy_visc2 = config_v_theta_eddy_visc2


      !
      ! Compute u (normal) velocity tendency for each edge (cell face)
      !

      tend_u(:,:) = 0.0

#ifdef LANL_FORMULATION
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
            tend_u(k,iEdge) = q - (ke(k,cell2) - ke(k,cell1)) / dcEdge(iEdge)
         end do
      end do

#endif

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

      do iEdge=1,grid % nEdgesSolve
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k=1,nVertLevels
            vorticity_abs = fEdge(iEdge) + (circulation(k,vertex1) + circulation(k,vertex2)) / &
                                           (areaTriangle(vertex1) + areaTriangle(vertex2))

            workpv = 2.0 * vorticity_abs / (h(k,cell1) + h(k,cell2))

            tend_u(k,iEdge) = workpv * vh(k,iEdge) - (ke(k,cell2) - ke(k,cell1)) / dcEdge(iEdge)
         end do
      end do
#endif


      !
      !  horizontal mixing for u
      !
      if ( h_mom_eddy_visc2 > 0.0 ) then
         do iEdge=1,grid % nEdgesSolve
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)

            do k=1,nVertLevels

               !
               ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
               !                    only valid for h_mom_eddy_visc2 == constant
               !
               u_diffusion =   ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
                              -( vorticity(k,vertex2) - vorticity(k,vertex1) ) / dvEdge(iEdge)
               u_diffusion = h_mom_eddy_visc2 * u_diffusion
 
               tend_u(k,iEdge) = tend_u(k,iEdge) + u_diffusion
            end do
         end do
      end if

      if ( h_mom_eddy_visc4 > 0.0 ) then

         allocate(delsq_divergence(nVertLevels, nCells+1))
         allocate(delsq_u(nVertLevels, nEdges+1))
         allocate(delsq_circulation(nVertLevels, nVertices+1))
         allocate(delsq_vorticity(nVertLevels, nVertices+1))

         delsq_u(:,:) = 0.0

         do iEdge=1,grid % nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)

            do k=1,nVertLevels
   
               !
               ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
               !                    only valid for h_mom_eddy_visc4 == constant
               !
               u_diffusion =   ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
                              -( vorticity(k,vertex2) - vorticity(k,vertex1) ) / dvEdge(iEdge)
 
               delsq_u(k,iEdge) = delsq_u(k,iEdge) + u_diffusion
            end do
         end do

         delsq_circulation(:,:) = 0.0
         do iEdge=1,nEdges
            do k=1,nVertLevels
               delsq_circulation(k,verticesOnEdge(1,iEdge)) = delsq_circulation(k,verticesOnEdge(1,iEdge)) - dcEdge(iEdge) * delsq_u(k,iEdge)
               delsq_circulation(k,verticesOnEdge(2,iEdge)) = delsq_circulation(k,verticesOnEdge(2,iEdge)) + dcEdge(iEdge) * delsq_u(k,iEdge)
            end do
         end do
         do iVertex=1,nVertices
            r = 1.0 / areaTriangle(iVertex)
            do k=1,nVertLevels
               delsq_vorticity(k,iVertex) = delsq_circulation(k,iVertex) * r
            end do
         end do

         delsq_divergence(:,:) = 0.0
         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1,nVertLevels
              delsq_divergence(k,cell1) = delsq_divergence(k,cell1) + delsq_u(k,iEdge)*dvEdge(iEdge)
              delsq_divergence(k,cell2) = delsq_divergence(k,cell2) - delsq_u(k,iEdge)*dvEdge(iEdge)
            end do
         end do
         do iCell = 1,nCells
            r = 1.0 / areaCell(iCell)
            do k = 1,nVertLevels
               delsq_divergence(k,iCell) = delsq_divergence(k,iCell) * r
            end do
         end do

         do iEdge=1,grid % nEdgesSolve
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)

            do k=1,nVertLevels

               !
               ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
               !                    only valid for h_mom_eddy_visc4 == constant
               !
               u_diffusion =   ( delsq_divergence(k,cell2)  - delsq_divergence(k,cell1) ) / dcEdge(iEdge)  &
                              -( delsq_vorticity(k,vertex2) - delsq_vorticity(k,vertex1) ) / dvEdge(iEdge)
 
               tend_u(k,iEdge) = tend_u(k,iEdge) - h_mom_eddy_visc4 * u_diffusion
            end do
         end do

         deallocate(delsq_divergence)
         deallocate(delsq_u)
         deallocate(delsq_circulation)
         deallocate(delsq_vorticity)

      end if


      !
      !  vertical advection for u
      !
      do iEdge=1,grid % nEdgesSolve
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         wdun(1) = 0.
         do k=2,nVertLevels
            wdun(k) =                                                                                  &
                     (ww(k,cell1)/(h(k,cell1)+h(k-1,cell1)) +ww(k,cell2)/(h(k,cell2)+h(k-1,cell2)))*   &
                      rdnu(k)*(u(k,iEdge)-u(k-1,iEdge))
         end do
         wdun(nVertLevels+1) = 0.

         do k=1,nVertLevels
            tend_u(k,iEdge) = tend_u(k,iEdge) - 0.5*(wdun(k+1)+wdun(k))
         end do
      end do


      !
      !  vertical mixing for u - 2nd order 
      !
      if ( v_mom_eddy_visc2 > 0.0 ) then

         do iEdge=1,grid % nEdgesSolve
   
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
    
            do k=2,nVertLevels-1
    
               z1 = 0.5*(geopotential(k-1,cell1)+geopotential(k-1,cell2))/gravity
               z2 = 0.5*(geopotential(k  ,cell1)+geopotential(k  ,cell2))/gravity
               z3 = 0.5*(geopotential(k+1,cell1)+geopotential(k+1,cell2))/gravity
               z4 = 0.5*(geopotential(k+2,cell1)+geopotential(k+2,cell2))/gravity
     
               zm = 0.5*(z1+z2)
               z0 = 0.5*(z2+z3)
               zp = 0.5*(z3+z4)
     
               tend_u(k,iEdge) = tend_u(k,iEdge) + v_mom_eddy_visc2*(                 &
                                  (u(k+1,iEdge)-u(k  ,iEdge))/(zp-z0)                 &
                                 -(u(k  ,iEdge)-u(k-1,iEdge))/(z0-zm) )/(0.5*(zp-zm))
            end do
         end do

      end if


!----------- rhs for theta

      tend_theta(:,:) = 0.


      !
      !  horizontal mixing for theta - we could combine this with advection directly (i.e. as a turbulent flux),
      !  but here we can also code in hyperdiffusion if we wish (2nd order at present)
      !
      if ( h_theta_eddy_visc2 > 0.0 ) then

         do iEdge=1,grid % nEdges
            cell1 = grid % cellsOnEdge % array(1,iEdge)
            cell2 = grid % cellsOnEdge % array(2,iEdge)

            do k=1,grid % nVertLevels
               theta_turb_flux = h_theta_eddy_visc2*prandtl*(theta(k,cell2) - theta(k,cell1))/dcEdge(iEdge)
               flux = dvEdge (iEdge) * h_edge(k,iEdge) * theta_turb_flux
               tend_theta(k,cell1) = tend_theta(k,cell1) + flux
               tend_theta(k,cell2) = tend_theta(k,cell2) - flux
            end do 

         end do 

      end if 

      if ( h_theta_eddy_visc4 > 0.0 ) then

         allocate(delsq_theta(nVertLevels, nCells+1))

         delsq_theta(:,:) = 0.

         do iEdge=1,grid % nEdges
            cell1 = grid % cellsOnEdge % array(1,iEdge)
            cell2 = grid % cellsOnEdge % array(2,iEdge)

            do k=1,grid % nVertLevels
               delsq_theta(k,cell1) = delsq_theta(k,cell1) + dvEdge(iEdge)*h_edge(k,iEdge)*(theta(k,cell2) - theta(k,cell1))/dcEdge(iEdge)
               delsq_theta(k,cell2) = delsq_theta(k,cell2) - dvEdge(iEdge)*h_edge(k,iEdge)*(theta(k,cell2) - theta(k,cell1))/dcEdge(iEdge)
            end do 

         end do 

         do iCell = 1, nCells
            r = 1.0 / areaCell(iCell)
            do k=1,nVertLevels
               delsq_theta(k,iCell) = delsq_theta(k,iCell) * r
            end do
         end do

         do iEdge=1,grid % nEdges
            cell1 = grid % cellsOnEdge % array(1,iEdge)
            cell2 = grid % cellsOnEdge % array(2,iEdge)

            do k=1,grid % nVertLevels
               theta_turb_flux = h_theta_eddy_visc4*prandtl*(delsq_theta(k,cell2) - delsq_theta(k,cell1))/dcEdge(iEdge)
               flux = dvEdge (iEdge) * theta_turb_flux

               tend_theta(k,cell1) = tend_theta(k,cell1) - flux
               tend_theta(k,cell2) = tend_theta(k,cell2) + flux
            end do 

         end do 

         deallocate(delsq_theta)

      end if 


      !
      !  horizontal advection for theta
      !

      if (config_theta_adv_order == 2) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1,grid % nVertLevels
               flux = dvEdge(iEdge) *  h_edge(k,iEdge) * u(k,iEdge) * (      &
                                      0.5*(theta(k,cell1) + theta(k,cell2)) )
               tend_theta(k,cell1) = tend_theta(k,cell1) - flux
               tend_theta(k,cell2) = tend_theta(k,cell2) + flux
            end do 
         end do 

      else if (config_theta_adv_order == 3) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
  
            do k=1,grid % nVertLevels
   
               d2fdx2_cell1 = deriv_two(1,1,iEdge) * theta(k,cell1)
               d2fdx2_cell2 = deriv_two(1,2,iEdge) * theta(k,cell2)
               do i=1, grid % nEdgesOnCell % array (cell1)
                  d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * theta(k,grid % CellsOnCell % array (i,cell1))
               end do
               do i=1, grid % nEdgesOnCell % array (cell2)
                  d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * theta(k,grid % CellsOnCell % array (i,cell2))
               end do
 
!  3rd order stencil
               if( u(k,iEdge) > 0) then
                  flux = dvEdge(iEdge) *  h_edge(k,iEdge) * u(k,iEdge) * (          &
                                         0.5*(theta(k,cell1) + theta(k,cell2))      &
                                         -(dcEdge(iEdge) **2) * (d2fdx2_cell1) / 6. )
               else
                  flux = dvEdge(iEdge) *  h_edge(k,iEdge) * u(k,iEdge) * (          &
                                         0.5*(theta(k,cell1) + theta(k,cell2))      &
                                         -(dcEdge(iEdge) **2) * (d2fdx2_cell2) / 6. )
               end if
   
               tend_theta(k,cell1) = tend_theta(k,cell1) - flux
               tend_theta(k,cell2) = tend_theta(k,cell2) + flux
 
            end do 
         end do 

      else  if (config_theta_adv_order == 4) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            do k=1,grid % nVertLevels
   
               d2fdx2_cell1 = deriv_two(1,1,iEdge) * theta(k,cell1)
               d2fdx2_cell2 = deriv_two(1,2,iEdge) * theta(k,cell2)
               do i=1, grid % nEdgesOnCell % array (cell1)
                  d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * theta(k,grid % CellsOnCell % array (i,cell1))
               end do
               do i=1, grid % nEdgesOnCell % array (cell2)
                  d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * theta(k,grid % CellsOnCell % array (i,cell2))
               end do
   
               flux = dvEdge(iEdge) *  h_edge(k,iEdge) * u(k,iEdge) * (          &
                                      0.5*(theta(k,cell1) + theta(k,cell2))      &
                                       -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12. )
  
               tend_theta(k,cell1) = tend_theta(k,cell1) - flux
               tend_theta(k,cell2) = tend_theta(k,cell2) + flux
            end do 
 
         end do
      end if


      !
      !  vertical advection plus diabatic term
      !  Note: we are also dividing through by the cell area after the horizontal flux divergence
      !
      do iCell = 1, nCells
         wdtn(1) = 0.
         do k=2,nVertLevels
            wdtn(k) =  ww(k,icell)*(fnm(k)*theta(k,iCell)+fnp(k)*theta(k-1,iCell))
         end do
         wdtn(nVertLevels+1) = 0.
         do k=1,nVertLevels
            tend_theta(k,iCell) = tend_theta(k,iCell)/areaCell(iCell) -rdnw(k)*(wdtn(k+1)-wdtn(k))
!!           tend_theta(k,iCell) = tend_theta(k) + h(k,iCell)*h_diabatic(k,iCell)
         end do
      end do


      !
      !  vertical mixing for theta - 2nd order 
      !
      if ( v_theta_eddy_visc2 > 0.0 ) then

         do iCell = 1, grid % nCellsSolve
            do k=2,nVertLevels-1
               z1 = geopotential(k-1,iCell)/gravity
               z2 = geopotential(k  ,iCell)/gravity
               z3 = geopotential(k+1,iCell)/gravity
               z4 = geopotential(k+2,iCell)/gravity
     
               zm = 0.5*(z1+z2)
               z0 = 0.5*(z2+z3)
               zp = 0.5*(z3+z4)
     
               tend_theta(k,iCell) = tend_theta(k,iCell) + v_theta_eddy_visc2*prandtl*h(k,iCell)*(  &
                                        (theta(k+1,iCell)-theta(k  ,iCell))/(zp-z0)                 &
                                       -(theta(k  ,iCell)-theta(k-1,iCell))/(z0-zm) )/(0.5*(zp-zm))
            end do
         end do

      end if

   end subroutine atmh_compute_dyn_tend

!---------------------------------------------------------------------------------------------------------

   subroutine atmh_advance_dynamics(tend, s, grid, small_step, number_small_steps, dt)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Advance the dry dynamics a small timestep (forward-backward integration)
   !
   ! Input: s - current model state
   !        tend - large-timestep tendency (d*/dt)
   !        grid - grid metadata
   !        dt   - timestep
   !
   ! Output: s - model state advanced a timestep dt
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (tend_type), intent(in) :: tend
      type (state_type), intent(inout) :: s
      type (mesh_type), intent(in) :: grid
      real (kind=RKIND), intent(in) :: dt
      integer, intent(in) :: small_step, number_small_steps

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j
      real (kind=RKIND) :: flux, vorticity_abs, h_vertex, workpv, q, upstream_bias

      integer :: nEdges, nVertices, nVertLevels
      real (kind=RKIND), dimension(:), pointer :: h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle, dpsdt, &
                                                  surface_pressure
      real (kind=RKIND), dimension(:,:), pointer :: vh, weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, tend_h, tend_u, &
                                                    circulation, vorticity, ke, pv_edge, geopotential, alpha, theta,       &
                                                    pressure, pressure_old, tend_theta, uhAvg, wwAvg, ww, u_old,           &
                                                    theta_old, h_edge_old, qtot, ww_old, cqu, h_old
      real (kind=RKIND), dimension(:,:,:), pointer :: scalar

!      real (kind=RKIND), pointer :: smext, p0, ptop
      real (kind=RKIND) :: smext, smdiv, p0, ptop
      real (kind=RKIND) :: tm, ptmp, he_old

      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge

      real (kind=RKIND), dimension( grid % nVertLevels + 1) :: wdtn

      real (kind=RKIND), dimension(:), pointer :: dnw, dbn, rdnw, dnu, fnm, fnp
      real (kind=RKIND) :: maxpdt,minpdt, maxww, minww
      integer :: maxpt,minpt

      h            => s % h % array
      u            => s % u % array
      h_edge       => s % h_edge % array
      theta        => s % theta % array

!      u_old        => s_old % u % array
!      h_edge_old   => s_old % h_edge % array
!      theta_old    => s_old % theta % array
!      ww_old      => s_old % ww % array
!      h_old       => s_old % h % array
      u_old        => grid % u_old % array
      h_edge_old   => grid % h_edge_old % array
      theta_old    => grid % theta_old % array
      ww_old      => grid % ww_old % array
      h_old       => grid % h_old % array

      geopotential => s % geopotential % array
      alpha        => s % alpha % array
      surface_pressure     => s % surface_pressure % array
      pressure     => s % pressure % array
      pressure_old => grid % pressure_old % array

      cellsOnEdge       => grid % cellsOnEdge % array
      dcEdge            => grid % dcEdge % array
      dvEdge            => grid % dvEdge % array
      areaCell          => grid % areaCell % array
      tend_h      => tend % h % array
      tend_u      => tend % u % array
      tend_theta      => tend % theta % array
                  

      uhAvg       => grid % uhAvg % array
      wwAvg       => grid % wwAvg % array
      dpsdt       => grid % dpsdt % array
      qtot        => grid % qtot % array
      cqu         => grid % cqu % array
      ww          => s % ww % array
      scalar      => s % scalars % array

      dnw         => grid % dnw % array
      dbn         => grid % dbn % array
      dnu         => grid % dnu % array
      rdnw        => grid % rdnw % array
      fnm         => grid % fnm % array
      fnp         => grid % fnp % array

!      p0          => grid % p0
!      ptop        => grid % ptop
!      smext       => grid % smext

      nVertLevels = grid % nVertLevels
      nEdges = grid % nEdges

      p0 = 1.e+05
      ptop = pressure(nVertLevels+1,1)
      smext = 0.1
      smdiv = 0.1

!       write(0,*) ' ptop in advance_dynamics ',ptop

!---  begin computations

!  we assume that the pressure, alpha, geopotential are already properly set
!  in first small step of a set, couple theta

      if(small_step == 1) then

        do iCell=1,grid % nCells
           do k=1,nVertLevels
              theta(k,iCell) = theta(k,iCell)*h(k,iCell)
           end do
        end do

        uhAvg = 0.
        wwAvg = 0.
        pressure_old(:,:) = pressure(:,:)
        dpsdt(:) = 0.

      end if

      !
      !  update horizontal momentum
      !

      do iEdge=1,grid % nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,nVertLevels
            u(k,iEdge) = u(k,iEdge) + dt*tend_u(k,iEdge)                 &
                               -(0.5*dt/dcEdge(iEdge))*(                 &
                 (geopotential(k+1,cell2)-geopotential(k+1,cell1))       &
                +(geopotential(k  ,cell2)-geopotential(k  ,cell1))       &
                +cqu(k,iEdge)*(alpha(k,cell2)+alpha(k,cell1))*           &
                       0.5*(pressure(k+1,cell2)-pressure(k+1,cell1)      &
                           +pressure(k  ,cell2)-pressure(k  ,cell1)))    &
                      -smext*dcEdge(iEdge)*(dpsdt(cell2)-dpsdt(cell1))/h_edge(k,iEdge)
         end do
      end do

      !
      !  calculate omega, update theta
      !

      tend_h = 0.
      do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1,nVertLevels
               flux = u(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge)
               tend_h(k,cell1) = tend_h(k,cell1) + flux
               tend_h(k,cell2) = tend_h(k,cell2) - flux
            end do 
            do k=1,nVertLevels
               uhAvg(k,iEdge) = uhAvg(k,iEdge) + u(k,iEdge) * h_edge(k,iEdge)
            end do 
      end do 

      do iCell=1, grid % nCells
        dpsdt(iCell) = 0.
        do k=1,nVertLevels
          dpsdt(iCell) = dpsdt(iCell) + dnw(k)*tend_h(k,iCell)
        end do
        dpsdt(iCell) = dpsdt(iCell)/areaCell(iCell)

        surface_pressure(iCell) = surface_pressure(iCell) + dt*dpsdt(iCell)

        do k=1,nVertLevels
          h(k,iCell) = (1.-dbn(k))*(p0-ptop)+dbn(k)*(surface_pressure(iCell)-ptop)
        end do

        ! omega calculation

        ww(1,iCell) = 0.
        do k=2, nVertLevels
          ww(k,iCell) = ww(k-1,iCell)-dnw(k-1)*(dbn(k-1)*dpsdt(iCell)+tend_h(k-1,iCell)/areaCell(iCell))
          wwAvg(k,iCell) = wwAvg(k,iCell) + ww(k,iCell)
        end do
        ww(nVertLevels+1,iCell) = 0.

        ! theta update  - theta should be coupled at this point...

        wdtn(1) = 0.
        do k = 2, nVertLevels
          wdtn(k) = (ww(k,iCell)-ww_old(k,iCell))*(fnm(k)*theta_old(k,iCell)+fnp(k)*theta_old(k-1,iCell))
        end do
        wdtn(nVertLevels+1) = 0.

        do k = 1, nVertLevels
          theta(k,iCell) = theta(k,iCell) + dt*tend_theta(k,iCell)
          theta(k,iCell) = theta(k,iCell) - dt*rdnw(k)*(wdtn(k+1)-wdtn(k))
        end do
      end do

      !
      ! add in perturbation horizontal flux divergence
      !

      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,nVertLevels
            h_edge(k,iEdge) = 0.5*(h(k,cell1)+h(k,cell2))  !  here is update of h_edge
            he_old = 0.5*(h_old(k,cell1)+h_old(k,cell2))
            flux = 0.5*(u(k,iEdge) * h_edge(k,iEdge) - u_old(k,iEdge) * he_old)* &
                        (theta_old(k,cell1)+theta_old(k,cell2))*dvEdge(iEdge)
            theta(k,cell1) = theta(k,cell1) - dt*flux/areaCell(cell1)
            theta(k,cell2) = theta(k,cell2) + dt*flux/areaCell(cell2)
         end do
      end do


      !  compute some diagnostics using the new state

      do iCell = 1, grid % nCells

        do k = nVertLevels,1,-1
          pressure(k,iCell) = pressure(k+1,iCell) - dnw(k)*(1.+qtot(k,iCell))*h(k,iCell) 
        end do

        do k=1, nVertLevels
          tm = (1.+1.61*scalar(1,k,iCell))*theta(k,iCell)/h(k,iCell)  !  assume scalar(1) is qv here?
          alpha(k,iCell) = (rgas/p0)*tm*  &
              (0.5*(pressure(k,iCell)+pressure(k+1,iCell))/p0)**cvpm
          geopotential(k+1,iCell) = geopotential(k,iCell) - dnw(k)*h(k,iCell)*alpha(k,iCell)
        end do

        if(small_step /= number_small_steps) then
          do k=1, nVertLevels+1
            ptmp = pressure(k,iCell)
            pressure(k,iCell) = pressure(k,iCell) + smdiv*(pressure(k,iCell)-pressure_old(k,iCell))
            pressure_old(k,iCell) = ptmp
          end do
        end if

      end do

!  if last small step of a set, decouple theta

      if(small_step == number_small_steps) then
        do iCell=1,grid % nCells
           do k=1,nVertLevels
              theta(k,iCell) = theta(k,iCell)/h(k,iCell)
           end do
        end do
        uhAvg = uhAvg/real(number_small_steps)
        wwAvg = wwAvg/real(number_small_steps)
      end if

   end subroutine atmh_advance_dynamics


   subroutine atmh_advance_scalars( tend, s_old, s_new, grid, dt)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed scalar tendencies
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (tend_type), intent(in) :: tend
      type (state_type), intent(in) :: s_old
      type (state_type), intent(inout) :: s_new
      type (mesh_type), intent(in) :: grid
      real (kind=RKIND) :: dt

      integer :: i, iCell, iEdge, k, iScalar, cell1, cell2, num_scalars
      real (kind=RKIND) :: flux, scalar_edge, d2fdx2_cell1, d2fdx2_cell2

      real (kind=RKIND), dimension(:,:,:), pointer :: scalar_old, scalar_new, scalar_tend
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND), dimension(:,:), pointer :: uhAvg, h_old, h_new, wwAvg
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, areaCell
      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension( s_old % num_scalars, grid % nVertLevels + 1 ) :: wdtn
      integer :: nVertLevels

      real (kind=RKIND), dimension(:), pointer :: fnm, fnp, rdnw
      real (kind=RKIND) :: coef_3rd_order

      num_scalars = s_old % num_scalars

      coef_3rd_order = 0.
      if (config_scalar_adv_order == 3) coef_3rd_order = 1.0
      if (config_scalar_adv_order == 3 .and. config_monotonic) coef_3rd_order = 0.25

      scalar_old  => s_old % scalars % array
      scalar_new  => s_new % scalars % array
      deriv_two   => grid % deriv_two % array
      uhAvg       => grid % uhAvg % array
      dvEdge      => grid % dvEdge % array
      dcEdge      => grid % dcEdge % array
      cellsOnEdge => grid % cellsOnEdge % array
      scalar_tend => tend % scalars % array
      h_old       => s_old % h % array
      h_new       => s_new % h % array
      wwAvg       => grid % wwAvg % array
      areaCell    => grid % areaCell % array

      fnm         => grid % fnm % array
      fnp         => grid % fnp % array
      rdnw        => grid % rdnw % array

      nVertLevels = grid % nVertLevels

      scalar_tend = 0.  !  testing purposes - we have no sources or sinks

      !
      ! Runge Kutta integration, so we compute fluxes from scalar_new values, update starts form scalar_old
      !
      !
      !  horizontal flux divergence, accumulate in scalar_tend

      if (config_scalar_adv_order == 2) then

         do iEdge=1,grid%nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1,grid % nVertLevels
               do iScalar=1,num_scalars
                  scalar_edge = 0.5 * (scalar_new(iScalar,k,cell1) + scalar_new(iScalar,k,cell2))
                  flux = uhAvg(k,iEdge) * dvEdge(iEdge)  * scalar_edge
                  scalar_tend(iScalar,k,cell1) = scalar_tend(iScalar,k,cell1) - flux/areaCell(cell1)
                  scalar_tend(iScalar,k,cell2) = scalar_tend(iScalar,k,cell2) + flux/areaCell(cell2)
               end do 
            end do 
         end do 

      else if (config_scalar_adv_order == 3) then

         do iEdge=1,grid%nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
  
            do k=1,grid % nVertLevels
   
               do iScalar=1,num_scalars
                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * scalar_new(iScalar,k,cell1)
                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * scalar_new(iScalar,k,cell2)
                  do i=1, grid % nEdgesOnCell % array (cell1)
                     d2fdx2_cell1 = d2fdx2_cell1 + &
                                    deriv_two(i+1,1,iEdge) * scalar_new(iScalar,k,grid % CellsOnCell % array (i,cell1))
                  end do
                  do i=1, grid % nEdgesOnCell % array (cell2)
                     d2fdx2_cell2 = d2fdx2_cell2 + &
                                    deriv_two(i+1,2,iEdge) * scalar_new(iScalar,k,grid % CellsOnCell % array (i,cell2))
                  end do

                  if (uhAvg(k,iEdge) > 0) then
                     flux = dvEdge(iEdge) * uhAvg(k,iEdge) * (          &
                                            0.5*(scalar_new(iScalar,k,cell1) + scalar_new(iScalar,k,cell2))      &
                                             -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                                             -(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12. )
                  else
                     flux = dvEdge(iEdge) *  uhAvg(k,iEdge) * (          &
                                            0.5*(scalar_new(iScalar,k,cell1) + scalar_new(iScalar,k,cell2))      &
                                             -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                                             +(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12. )
                  end if

! old version of the above code, with coef_3rd_order assumed to be 1.0
!                  if (uhAvg(k,iEdge) > 0) then
!                     flux = dvEdge(iEdge) * uhAvg(k,iEdge) * (          &
!                                            0.5*(scalar_new(iScalar,k,cell1) + scalar_new(iScalar,k,cell2))      &
!                                            -(dcEdge(iEdge) **2) * (d2fdx2_cell1) / 6. )
!                  else
!                     flux = dvEdge(iEdge) *  uhAvg(k,iEdge) * (          &
!                                            0.5*(scalar_new(iScalar,k,cell1) + scalar_new(iScalar,k,cell2))      &
!                                            -(dcEdge(iEdge) **2) * (d2fdx2_cell2) / 6. )
!                  end if
   
                  scalar_tend(iScalar,k,cell1) = scalar_tend(iScalar,k,cell1) - flux/areaCell(cell1)
                  scalar_tend(iScalar,k,cell2) = scalar_tend(iScalar,k,cell2) + flux/areaCell(cell2)
 
               end do 
            end do 
         end do 

      else  if (config_scalar_adv_order == 4) then

         do iEdge=1,grid%nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            do k=1,grid % nVertLevels
   
               do iScalar=1,num_scalars
                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * scalar_new(iScalar,k,cell1)
                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * scalar_new(iScalar,k,cell2)
                  do i=1, grid % nEdgesOnCell % array (cell1)
                        d2fdx2_cell1 = d2fdx2_cell1 + &
                                       deriv_two(i+1,1,iEdge) * scalar_new(iScalar,k,grid % CellsOnCell % array (i,cell1))
                  end do
                  do i=1, grid % nEdgesOnCell % array (cell2)
                     d2fdx2_cell2 = d2fdx2_cell2 + &
                                    deriv_two(i+1,2,iEdge) * scalar_new(iScalar,k,grid % CellsOnCell % array (i,cell2))
                  end do
      
                  flux = dvEdge(iEdge) *  uhAvg(k,iEdge) * (          &
                                         0.5*(scalar_new(iScalar,k,cell1) + scalar_new(iScalar,k,cell2))      &
                                          -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12. )
     
                  scalar_tend(iScalar,k,cell1) = scalar_tend(iScalar,k,cell1) - flux/areaCell(cell1)
                  scalar_tend(iScalar,k,cell2) = scalar_tend(iScalar,k,cell2) + flux/areaCell(cell2)
               end do 
            end do 
 
         end do
      end if


      !
      !  vertical flux divergence
      !

      do iCell=1,grid % nCells

        wdtn(:,1) = 0.
        do k = 2, nVertLevels
          do iScalar=1,num_scalars
            wdtn(iScalar,k) = wwAvg(k,iCell)*(fnm(k)*scalar_new(iScalar,k,iCell)+fnp(k)*scalar_new(iScalar,k-1,iCell))
          end do
        end do
        wdtn(:,nVertLevels+1) = 0.

         do k=1,grid % nVertLevels   ! Could be nVertLevelsSolve?
            do iScalar=1,num_scalars
              scalar_new(iScalar,k,iCell) = (   scalar_old(iScalar,k,iCell)*h_old(k,iCell) &
                    + dt*( scalar_tend(iScalar,k,iCell) -rdnw(k)*(wdtn(iScalar,k+1)-wdtn(iScalar,k)) ) )/h_new(k,iCell)
                                                                                        
            end do
         end do
      end do

   end subroutine atmh_advance_scalars


   subroutine atmh_advance_scalars_mono(tend, s_old, s_new, grid, dt, rk_step, rk_order)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed scalar tendencies
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (tend_type), intent(in) :: tend
      type (state_type), intent(in) :: s_old
      type (state_type), intent(inout) :: s_new
      type (mesh_type), intent(in) :: grid
      real (kind=RKIND), intent(in) :: dt
      integer, intent(in) :: rk_step, rk_order

      type (block_type), pointer :: block

      integer :: i, iCell, iEdge, k, iScalar, cell_upwind, cell1, cell2, num_scalars
      real (kind=RKIND) :: flux, scalar_edge, d2fdx2_cell1, d2fdx2_cell2
      real (kind=RKIND) :: fdir, flux_upwind, h_flux_upwind, s_upwind

      real (kind=RKIND), dimension(:,:,:), pointer :: scalar_old, scalar_new, scalar_tend
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND), dimension(:,:), pointer :: uhAvg, h_old, h_new, wwAvg
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, areaCell
      integer, dimension(:,:), pointer :: cellsOnEdge

      type (field3DReal), pointer :: tempField
      type (field3DReal), target :: tempFieldTarget

      real (kind=RKIND), dimension( s_old % num_scalars, grid % nEdges+1) :: h_flux
      real (kind=RKIND), dimension(2, s_old % num_scalars, grid % nCells+1) :: v_flux, v_flux_upwind, s_update
      real (kind=RKIND), dimension(2, s_old % num_scalars, grid % nCells+1), target :: scale_out, scale_in
      real (kind=RKIND), dimension( s_old % num_scalars ) :: s_max, s_min, s_max_update, s_min_update

      integer :: nVertLevels, km0, km1, ktmp, kcp1, kcm1

      real (kind=RKIND), dimension(:), pointer :: fnm, fnp, rdnw
      real (kind=RKIND), parameter :: eps=1.e-20
      real (kind=RKIND) :: coef_3rd_order

      block => grid % block

      num_scalars = s_old % num_scalars

      scalar_old  => s_old % scalars % array
      scalar_new  => s_new % scalars % array
      deriv_two   => grid % deriv_two % array
      uhAvg       => grid % uhAvg % array
      dvEdge      => grid % dvEdge % array
      dcEdge      => grid % dcEdge % array
      cellsOnEdge => grid % cellsOnEdge % array
      scalar_tend => tend % scalars % array
      h_old       => s_old % h % array
      h_new       => s_new % h % array
      wwAvg       => grid % wwAvg % array
      areaCell    => grid % areaCell % array

      fnm         => grid % fnm % array
      fnp         => grid % fnp % array
      rdnw        => grid % rdnw % array

      nVertLevels = grid % nVertLevels

      scalar_tend = 0.  !  testing purposes - we have no sources or sinks

      !
      ! Runge Kutta integration, so we compute fluxes from scalar_new values, update starts from scalar_old
      !

      km1 = 1
      km0 = 2
      v_flux(km1,:,:) = 0.
      v_flux_upwind(km1,:,:) = 0.
      scale_out(:,:,:) = 1.
      scale_in(:,:,:) = 1.

      coef_3rd_order = 0.
      if (config_scalar_adv_order == 3) coef_3rd_order = 1.0
      if (config_scalar_adv_order == 3 .and. config_monotonic) coef_3rd_order = 0.25

      do k = 1, grid % nVertLevels
         kcp1 = min(k+1,grid % nVertLevels)
         kcm1 = max(k-1,1)

!  vertical flux

         do iCell=1,grid % nCells

            if (k < grid % nVertLevels) then
               cell_upwind = k
               if (wwAvg(k+1,iCell) >= 0) cell_upwind = k+1
               do iScalar=1,num_scalars
                  v_flux(km0,iScalar,iCell) = dt * wwAvg(k+1,iCell) *   &
                       (fnm(k+1) * scalar_new(iScalar,k+1,iCell) + fnp(k+1) * scalar_new(iScalar,k,iCell))
                  v_flux_upwind(km0,iScalar,iCell) = dt * wwAvg(k+1,iCell) * scalar_old(iScalar,cell_upwind,iCell)
                  v_flux(km0,iScalar,iCell) = v_flux(km0, iScalar,iCell) - v_flux_upwind(km0,iScalar,iCell)
!                  v_flux(iScalar,iCell,km0) = 0.  ! use only upwind - for testing
                  s_update(km0,iScalar,iCell) = scalar_old(iScalar,k,iCell) * h_old(k,iCell)  &
                            - rdnw(k) * (v_flux_upwind(km0,iScalar,iCell) - v_flux_upwind(km1,iScalar,iCell))
               end do
            else
               do iScalar=1,num_scalars
                  v_flux(km0,iScalar,iCell) = 0.
                  v_flux_upwind(km0,iScalar,iCell) = 0.
                  s_update(km0,iScalar,iCell) = scalar_old(iScalar,k,iCell) * h_old(k,iCell)  &
                            - rdnw(k) * (v_flux_upwind(km0,iScalar,iCell) - v_flux_upwind(km1,iScalar,iCell))
               end do
            end if

         end do

! horizontal flux

         if (config_scalar_adv_order == 2) then

            do iEdge=1,grid%nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               cell_upwind = cell2
               if (uhAvg(k,iEdge) >= 0) cell_upwind = cell1
               do iScalar=1,num_scalars
                  scalar_edge = 0.5 * (scalar_new(iScalar,k,cell1) + scalar_new(iScalar,k,cell2))
                  h_flux(iScalar,iEdge) = dt * uhAvg(k,iEdge) * dvEdge(iEdge) * scalar_edge
                  h_flux_upwind = dt * uhAvg(k,iEdge) * dvEdge(iEdge) * scalar_old(iScalar,k,cell_upwind)
                  h_flux(iScalar,iEdge) = h_flux(iScalar,iEdge) - h_flux_upwind
!                  h_flux(iScalar,iEdge) = 0.  ! use only upwind - for testing
                  s_update(km0,iScalar,cell1) = s_update(km0,iScalar,cell1) - h_flux_upwind / grid % areaCell % array(cell1)
                  s_update(km0,iScalar,cell2) = s_update(km0,iScalar,cell2) + h_flux_upwind / grid % areaCell % array(cell2)
               end do 
            end do 

         else if (config_scalar_adv_order >= 3) then

            do iEdge=1,grid%nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               cell_upwind = cell2
               if (uhAvg(k,iEdge) >= 0) cell_upwind = cell1
               do iScalar=1,num_scalars

                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * scalar_new(iScalar,k,cell1)
                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * scalar_new(iScalar,k,cell2)
                  do i=1, grid % nEdgesOnCell % array (cell1)
                     d2fdx2_cell1 = d2fdx2_cell1 + &
                                    deriv_two(i+1,1,iEdge) * scalar_new(iScalar,k,grid % CellsOnCell % array (i,cell1))
                  end do
                  do i=1, grid % nEdgesOnCell % array (cell2)
                     d2fdx2_cell2 = d2fdx2_cell2 + &
                                    deriv_two(i+1,2,iEdge) * scalar_new(iScalar,k,grid % CellsOnCell % array (i,cell2))
                  end do
   
                  if (uhAvg(k,iEdge) > 0) then
                     flux = dvEdge(iEdge) * uhAvg(k,iEdge) * (          &
                                            0.5*(scalar_new(iScalar,k,cell1) + scalar_new(iScalar,k,cell2))      &
                                             -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                                             -(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12. )
                  else
                     flux = dvEdge(iEdge) *  uhAvg(k,iEdge) * (          &
                                            0.5*(scalar_new(iScalar,k,cell1) + scalar_new(iScalar,k,cell2))      &
                                             -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                                             +(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12. )
                  end if
  
                  h_flux(iScalar,iEdge) = dt * flux
                  h_flux_upwind = dt * uhAvg(k,iEdge) * dvEdge(iEdge) * scalar_old(iScalar,k,cell_upwind)
                  h_flux(iScalar,iEdge) = h_flux(iScalar,iEdge) - h_flux_upwind
!                  h_flux(iScalar,iEdge) = 0.  ! use only upwind - for testing
                  s_update(km0,iScalar,cell1) = s_update(km0,iScalar,cell1) - h_flux_upwind / grid % areaCell % array(cell1)
                  s_update(km0,iScalar,cell2) = s_update(km0,iScalar,cell2) + h_flux_upwind / grid % areaCell % array(cell2)
               end do 
            end do 

         end if


         if ( (rk_step == rk_order) .and. (config_monotonic .or. config_positive_definite) ) then   

!*************************************************************************************************************
!---  limiter - we limit horizontal and vertical fluxes on level k 
!---  (these are h fluxes contributing to level k scalars, and v flux contributing to level k, k-1 scalars)

            do iCell=1,grid % nCells
  
               do iScalar=1,num_scalars
   
                  s_max(iScalar) = max(scalar_old(iScalar,k,iCell), scalar_old(iScalar,kcp1,iCell), scalar_old(iScalar,kcm1,iCell))
                  s_min(iScalar) = min(scalar_old(iScalar,k,iCell), scalar_old(iScalar,kcp1,iCell), scalar_old(iScalar,kcm1,iCell))
                  s_max_update(iScalar) = s_update(km0,iScalar,iCell)
                  s_min_update(iScalar) = s_update(km0,iScalar,iCell)
    
                  ! add in vertical flux to get max and min estimate
                  s_max_update(iScalar) = s_max_update(iScalar)  &
                     - rdnw(k) * (max(0.0_RKIND,v_flux(km0,iScalar,iCell)) - min(0.0_RKIND,v_flux(km1,iScalar,iCell)))
                  s_min_update(iScalar) = s_min_update(iScalar)  &
                     - rdnw(k) * (min(0.0_RKIND,v_flux(km0,iScalar,iCell)) - max(0.0_RKIND,v_flux(km1,iScalar,iCell)))
    
               end do
   
               do i = 1, grid % nEdgesOnCell % array(iCell)  ! go around the edges of each cell
                  do iScalar=1,num_scalars
    
                     s_max(iScalar)  = max(scalar_old(iScalar,k,grid % cellsOnCell % array(i,iCell)), s_max(iScalar))
                     s_min(iScalar)  = min(scalar_old(iScalar,k,grid % cellsOnCell % array(i,iCell)), s_min(iScalar))
     
                     iEdge = grid % EdgesOnCell % array (i,iCell)
                     if (iCell == cellsOnEdge(1,iEdge)) then
                        fdir = 1.0
                     else
                        fdir = -1.0
                     end if
                     flux = -fdir * h_flux(iScalar,iEdge)/grid % areaCell % array(iCell)
                     s_max_update(iScalar) = s_max_update(iScalar) + max(0.0_RKIND,flux)
                     s_min_update(iScalar) = s_min_update(iScalar) + min(0.0_RKIND,flux)
  
                  end do
   
               end do
   
               if( config_positive_definite ) s_min(:) = 0.
   
               do iScalar=1,num_scalars
                  scale_out (km0,iScalar,iCell) = 1.
                  scale_in (km0,iScalar,iCell) = 1.
                  s_max_update (iScalar) =  s_max_update (iScalar) / h_new (k,iCell)
                  s_min_update (iScalar) =  s_min_update (iScalar) / h_new (k,iCell)
                  s_upwind = s_update(km0,iScalar,iCell) / h_new(k,iCell)
                  if ( s_max_update(iScalar) > s_max(iScalar) .and. config_monotonic)   &
                     scale_in (km0,iScalar,iCell) = max(0.0_RKIND,(s_max(iScalar)-s_upwind)/(s_max_update(iScalar)-s_upwind+eps))
                  if ( s_min_update(iScalar) < s_min(iScalar) )   &
                     scale_out (km0,iScalar,iCell) = max(0.0_RKIND,(s_upwind-s_min(iScalar))/(s_upwind-s_min_update(iScalar)+eps))
                end do
  
            end do ! end loop over cells to compute scale factor


       tempField => tempFieldTarget
       tempField % block => block
       tempField % dimSizes(1) = 2
       tempField % dimSizes(2) = num_scalars
       tempField % dimSizes(3) = grid % nCells
       tempField % sendList => block % parinfo % cellsToSend
       tempField % recvList => block % parinfo % cellsToRecv
       tempField % copyList => block % parinfo % cellsToCopy
       tempField % prev => null()
       tempField % next => null()

       tempField % array => scale_in
       call mpas_dmpar_exch_halo_field(tempField)

       tempField % array => scale_out
       call mpas_dmpar_exch_halo_field(tempField)
     
       ! rescale the horizontal fluxes
 
            do iEdge = 1, grid % nEdges
               cell1 = grid % cellsOnEdge % array(1,iEdge)
               cell2 = grid % cellsOnEdge % array(2,iEdge)
               do iScalar=1,num_scalars
                  flux = h_flux(iScalar,iEdge)
                  if (flux > 0) then
                     flux = flux * min(scale_out(km0,iScalar,cell1), scale_in(km0,iScalar,cell2))
                  else
                     flux = flux * min(scale_in(km0,iScalar,cell1), scale_out(km0,iScalar,cell2))
                  end if
                  h_flux(iScalar,iEdge) = flux
               end do
            end do
 
       ! rescale the vertical flux
 
            do iCell=1,grid % nCells
               do iScalar=1,num_scalars
                  flux =  v_flux(km1,iScalar,iCell)
                  if (flux > 0) then
                     flux = flux * min(scale_out(km0,iScalar,iCell), scale_in(km1,iScalar,iCell))
                  else
                     flux = flux * min(scale_in(km0,iScalar,iCell), scale_out(km1,iScalar,iCell))
                  end if
                  v_flux(km1,iScalar,iCell) = flux
               end do
            end do

!  end of limiter
!*******************************************************************************************************************

         end if

!---  update

         do iCell=1,grid % nCells
            !  add in upper vertical flux that was just renormalized
            do iScalar=1,num_scalars
               s_update(km0,iScalar,iCell) = s_update(km0,iScalar,iCell) + rdnw(k) * v_flux(km1,iScalar,iCell)
               if (k > 1) s_update(km1,iScalar,iCell) = s_update(km1,iScalar,iCell) - rdnw(k-1)*v_flux(km1,iScalar,iCell)
            end do
         end do
 
         do iEdge=1,grid%nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do iScalar=1,num_scalars
               s_update(km0,iScalar,cell1) = s_update(km0,iScalar,cell1) - &
                   h_flux(iScalar,iEdge) / grid % areaCell % array(cell1)
               s_update(km0,iScalar,cell2) = s_update(km0,iScalar,cell2) + &
                   h_flux(iScalar,iEdge) / grid % areaCell % array(cell2)
            end do 
         end do 
 
         ! decouple from mass
         if (k > 1) then
            do iCell=1,grid % nCells
               do iScalar=1,num_scalars
                  s_update(km1,iScalar,iCell) = s_update(km1,iScalar,iCell) / h_new(k-1,iCell)
               end do
            end do
 
            do iCell=1,grid % nCells
               do iScalar=1,num_scalars
                  scalar_new(iScalar,k-1,iCell) = s_update(km1,iScalar,iCell) 
               end do
            end do
         end if
 
         ktmp = km1
         km1 = km0
         km0 = ktmp

      end do

      do iCell=1,grid % nCells
         do iScalar=1,num_scalars
            scalar_new(iScalar,grid % nVertLevels,iCell) = s_update(km1,iScalar,iCell) / h_new(grid%nVertLevels,iCell)
         end do
      end do

   end subroutine atmh_advance_scalars_mono


   subroutine atmh_compute_solve_diagnostics(dt, s, grid)
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
      real (kind=RKIND) :: flux, vorticity_abs, h_vertex, workpv, r

      integer :: nCells, nEdges, nVertices, nVertLevels
      real (kind=RKIND), dimension(:), pointer :: h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: vh, weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, tend_h, tend_u, &
                                                    circulation, vorticity, ke, pv_edge, pv_vertex, pv_cell, gradPVn, gradPVt, &
                                                    divergence
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge


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

      !
      ! Compute height on cell edges at velocity locations
      !
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,nVertLevels
            h_edge(k,iEdge) = 0.5 * (h(k,cell1) + h(k,cell2))
         end do
      end do

      !
      ! Compute circulation and relative vorticity at each vertex
      !
      circulation(:,:) = 0.0
      do iEdge=1,nEdges
         do k=1,nVertLevels
            circulation(k,verticesOnEdge(1,iEdge)) = circulation(k,verticesOnEdge(1,iEdge)) - dcEdge(iEdge) * u(k,iEdge)
            circulation(k,verticesOnEdge(2,iEdge)) = circulation(k,verticesOnEdge(2,iEdge)) + dcEdge(iEdge) * u(k,iEdge)
         end do
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
         do k=1,nVertLevels
           divergence(k,cell1) = divergence(k,cell1) + u(k,iEdge)*dvEdge(iEdge)
           divergence(k,cell2) = divergence(k,cell2) - u(k,iEdge)*dvEdge(iEdge)
         end do
      end do
      do iCell = 1,nCells
        r = 1.0 / areaCell(iCell)
        do k = 1,nVertLevels
           divergence(k,iCell) = divergence(k,iCell) * r
        end do
      end do


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
            do k = 1,nVertLevels
               v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
            end do
         end do
      end do


      ! tdr
      !
      ! Compute height at vertices, pv at vertices, and average pv to edge locations
      !  ( this computes pv_vertex at all vertices bounding real cells )
      !
      do iVertex = 1,nVertices
         do k=1,nVertLevels
            h_vertex = 0.0
            do i=1,grid % vertexDegree
               h_vertex = h_vertex + h(k,cellsOnVertex(i,iVertex)) * kiteAreasOnVertex(i,iVertex)
            end do
            h_vertex = h_vertex / areaTriangle(iVertex)

            pv_vertex(k,iVertex) = (fVertex(iVertex) + vorticity(k,iVertex)) / h_vertex
         end do
      end do
      ! tdr


      ! tdr
      !
      ! Compute gradient of PV in the tangent direction
      !   ( this computes gradPVt at all edges bounding real cells )
      !
      do iEdge = 1,nEdges
         do k = 1,nVertLevels
           gradPVt(k,iEdge) = (pv_vertex(k,verticesOnEdge(2,iEdge)) - pv_vertex(k,verticesOnEdge(1,iEdge))) / &
                              dvEdge(iEdge)
         end do
      end do

      ! tdr
      !
      ! Compute pv at the edges
      !   ( this computes pv_edge at all edges bounding real cells )
      !
      pv_edge(:,:) = 0.0
      do iVertex = 1,nVertices
        do i=1,grid % vertexDegree
           iEdge = edgesOnVertex(i,iVertex)
           do k=1,nVertLevels
              pv_edge(k,iEdge) =  pv_edge(k,iEdge)  + 0.5 * pv_vertex(k,iVertex)
           end do
        end do
      end do
      ! tdr

      ! tdr
      !
      ! Modify PV edge with upstream bias. 
      !
      do iEdge = 1,nEdges
         do k = 1,nVertLevels
           pv_edge(k,iEdge) = pv_edge(k,iEdge) - config_apvm_upwinding * v(k,iEdge) * dt * gradPVt(k,iEdge)
         end do
      end do


      ! tdr
      !
      ! Compute pv at cell centers
      !    ( this computes pv_cell for all real cells )
      !
      pv_cell(:,:) = 0.0
      do iVertex = 1, nVertices
       do i=1,grid % vertexDegree
          iCell = cellsOnVertex(i,iVertex)
          do k = 1,nVertLevels
             pv_cell(k,iCell) = pv_cell(k,iCell) + kiteAreasOnVertex(i, iVertex) * pv_vertex(k, iVertex) / areaCell(iCell)
          end do
       end do
      end do
      ! tdr

      ! tdr
      !
      ! Compute gradient of PV in normal direction
      !   (tdr: 2009-10-02: this is not correct because the pv_cell in the halo is not correct)
      !
      gradPVn(:,:) = 0.0
      do iEdge = 1,nEdges
         do k = 1,nVertLevels
            gradPVn(k,iEdge) = (pv_cell(k,cellsOnEdge(2,iEdge)) - pv_cell(k,cellsOnEdge(1,iEdge))) / &
                                 dcEdge(iEdge)
         end do
      end do
      ! tdr

      ! Modify PV edge with upstream bias.
      !
     do iEdge = 1,nEdges
        do k = 1,nVertLevels
          pv_edge(k,iEdge) = pv_edge(k,iEdge) - config_apvm_upwinding * u(k,iEdge) * dt * gradPVn(k,iEdge)
        end do
     end do


   end subroutine atmh_compute_solve_diagnostics


   subroutine atmh_compute_w (s_new, s_old, grid, dt )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Compute (diagnose) vertical velocity (used by physics)
   !
   ! Input: s_new - current model state
   !        s_old - previous model state
   !        grid - grid metadata
   !        dt - timestep between new and old
   !
   ! Output: w  (m/s)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (state_type), intent(inout) :: s_new
      type (state_type), intent(in) :: s_old
      type (mesh_type), intent(inout) :: grid

      real (kind=RKIND), intent(in) :: dt

      real (kind=RKIND), dimension(:,:), pointer :: geo_new, geo_old, u, ww, h, w
      real (kind=RKIND), dimension(:), pointer :: dvEdge, rdnw, fnm, fnp

      integer :: iEdge, iCell, k, cell1, cell2
      real (kind=RKIND), dimension( grid % nVertlevels + 1 ) :: wdwn
      real (kind=RKIND) :: flux

      geo_new => s_new % geopotential % array      
      geo_old => s_old % geopotential % array      
      u => s_new % u % array 
      ww => s_new % ww % array
      h => s_new % h % array
      w => s_new % w % array
      dvEdge => grid % dvEdge % array
      rdnw => grid % rdnw % array
      fnm => grid % fnm % array
      fnp => grid % fnp % array

      w = 0.

      do iCell=1, grid % nCellsSolve
        wdwn(1) = 0.
        do k=2,grid % nVertlevels + 1
          wdwn(k) = (0.5*(ww(k,iCell)+ww(k-1,iCell))/h(k-1,iCell))   &
                     *rdnw(k-1)*(geo_new(k,iCell)-geo_new(k-1,iCell))
        enddo
        w(1,iCell) = 0.
        do k=2, grid % nVertLevels
          w(k,iCell) = (((geo_new(k,iCell)-geo_old(k,iCell))/dt)+ &
                          fnm(k)*wdwn(k+1)+fnp(k)*wdwn(k))/gravity
        enddo
        k = grid % nVertLevels + 1
        w(k,iCell) = ((geo_new(k,iCell)-geo_old(k,iCell))/dt)/gravity
      enddo

      do iEdge=1, grid % nEdges
        cell1 = grid % cellsOnEdge % array (1,iEdge)
        cell2 = grid % cellsOnEdge % array (2,iEdge)
        if (cell1 <= grid % nCellsSolve .or. cell2 <= grid % nCellsSolve ) then
          do k=2, grid % nVertLevels
            flux = 0.25*(u(k,iEdge)+u(k-1,iEdge)*(geo_new(k,cell2)-geo_new(k,cell1))) * dvEdge(iEdge) / gravity
            w(k,cell1) = w(k,cell1) + flux/ grid % AreaCell % array (cell1)
            w(k,cell2) = w(k,cell2) + flux/ grid % AreaCell % array (cell2)
          enddo
          k = 1
          flux = 0.5*(1.5*u(1,iEdge)-0.5*u(2,iEdge)*(geo_new(k,cell2)-geo_new(k,cell1))) * dvEdge(iEdge) / gravity
          w(k,cell1) = w(k,cell1) + flux/ grid % AreaCell % array (cell1)
          w(k,cell2) = w(k,cell2) + flux/ grid % AreaCell % array (cell2)

          k = grid % nVertLevels + 1
          flux = 0.5*(1.5*u(k-1,iEdge)-0.5*u(k-2,iEdge)*(geo_new(k,cell2)-geo_new(k,cell1))) * dvEdge(iEdge) / gravity
          w(k,cell1) = w(k,cell1) + flux/ grid % AreaCell % array (cell1)
          w(k,cell2) = w(k,cell2) + flux/ grid % AreaCell % array (cell2)

        end if
      end do

   end subroutine atmh_compute_w

end module atmh_time_integration
