!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_split
!
!> \brief MPAS ocean split explicit time integration scheme
!> \author Doug Jacobsen
!> \date   26 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the routine for the split explicit
!>  time integration scheme
!
!-----------------------------------------------------------------------


module ocn_time_integration_split

   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_dmpar
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer

   use ocn_tendency

   use ocn_equation_of_state
   use ocn_vmix


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

   public :: ocn_time_integrator_split

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_split
!
!> \brief MPAS ocean split explicit time integration scheme
!> \author Doug Jacobsen
!> \date   26 September 2011
!> \version SVN:$Id:$
!> \details
!>  This routine integrates a single time step (dt) using a
!>  split explicit time integrator.
!
!-----------------------------------------------------------------------

subroutine ocn_time_integrator_split(domain, dt)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Advance model state forward in time by the specified time step using 
   !   Split_Explicit timestepping scheme
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
      integer :: iCell, i,k,j, iEdge, cell1, cell2, split_explicit_step, split, &
        eoe, oldBtrSubcycleTime, newBtrSubcycleTime, uPerpTime, BtrCorIter, &
        n_bcl_iter(config_n_ts_iter), &
        vertex1, vertex2, iVertex

      type (block_type), pointer :: block
      real (kind=RKIND) :: uhSum, hSum, sshEdge, flux, &
         uPerp, uCorr, tracerTemp, coef, FBtr_coeff, sshCell1, sshCell2
      real (kind=RKIND), dimension(:), pointer :: sshNew

      integer :: num_tracers, ucorr_coef, err
      real (kind=RKIND), dimension(:,:), pointer :: &
        u, h, h_edge, ke_edge, vertViscTopOfEdge, vertDiffTopOfCell
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      integer, dimension(:), pointer :: & 
        maxLevelCell, maxLevelEdgeTop
      real (kind=RKIND), dimension(:), allocatable:: A,C,uTemp, hNew
      real (kind=RKIND), dimension(:,:), allocatable:: tracersTemp

      call mpas_timer_start("split_explicit_timestep")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !  Prep variables before first iteration
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      block => domain % blocklist
      do while (associated(block))

         do iEdge=1,block % mesh % nEdges

            ! The baroclinic velocity needs be recomputed at the beginning of a 
            ! timestep because the implicit vertical mixing is conducted on the
            ! total u.  We keep uBtr from the previous timestep.
              block % state % time_levs(1) % state % uBcl % array(:,iEdge) &
            = block % state % time_levs(1) % state % u % array(:,iEdge) &
            - block % state % time_levs(1) % state % uBtr % array(iEdge)

              block % state % time_levs(2) % state % u % array(:,iEdge) &
            = block % state % time_levs(1) % state % u % array(:,iEdge)

              block % state % time_levs(2) % state % uBcl % array(:,iEdge) &
            = block % state % time_levs(1) % state % uBcl % array(:,iEdge)

         enddo ! iEdge

         ! Initialize * variables that are used compute baroclinic tendencies below.
           block % state % time_levs(2) % state % ssh % array(:) &
         = block % state % time_levs(1) % state % ssh % array(:)

           block % state % time_levs(2) % state % h_edge % array(:,:) &
         = block % state % time_levs(1) % state % h_edge % array(:,:)

         do iCell=1,block % mesh % nCells  ! couple tracers to h
           ! change to maxLevelCell % array(iCell) ?
           do k=1,block % mesh % nVertLevels

                block % state % time_levs(2) % state % tracers % array(:,k,iCell) & 
              = block % state % time_levs(1) % state % tracers % array(:,k,iCell) 
            end do

         end do

         block => block % next
      end do


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN large iteration loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      n_bcl_iter = config_n_bcl_iter_mid
      n_bcl_iter(1) = config_n_bcl_iter_beg
      n_bcl_iter(config_n_ts_iter) = config_n_bcl_iter_end

      do split_explicit_step = 1, config_n_ts_iter
! ---  update halos for diagnostic variables

        block => domain % blocklist
        do while (associated(block))

           call mpas_dmpar_exch_halo_field2d_real(domain % dminfo, block % state % time_levs(2) % state % pv_edge % array(:,:), &
                                            block % mesh % nVertLevels, block % mesh % nEdges, &
                                            block % parinfo % edgesToSend, block % parinfo % edgesToRecv)

           if (config_h_mom_eddy_visc4 > 0.0) then
              call mpas_dmpar_exch_halo_field2d_real(domain % dminfo, block % state % time_levs(2) % state % divergence % array(:,:), &
                                               block % mesh % nVertLevels, block % mesh % nCells, &
                                               block % parinfo % cellsToSend, block % parinfo % cellsToRecv)
              call mpas_dmpar_exch_halo_field2d_real(domain % dminfo, block % state % time_levs(2) % state % vorticity % array(:,:), &
                                               block % mesh % nVertLevels, block % mesh % nVertices, &
                                               block % parinfo % verticesToSend, block % parinfo % verticesToRecv)
           end if

           block => block % next
        end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !  Stage 1: Baroclinic velocity (3D) prediction, explicit with long timestep
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! compute velocity tendencies, T(u*,w*,p*)

      block => domain % blocklist
      do while (associated(block))
         if (.not.config_implicit_vertical_mix) then
            call ocn_vmix_coefs(block % mesh, block % state % time_levs(2) % state, block % diagnostics, err)
         end if
         call ocn_tend_u(block % tend, block % state % time_levs(2) % state , block % diagnostics, block % mesh)
         call enforce_boundaryEdge(block % tend, block % mesh)
         block => block % next
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN baroclinic iterations on linear Coriolis term
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j=1,n_bcl_iter(split_explicit_step)

         ! Use this G coefficient to avoid an if statement within the iEdge loop.
         if     (trim(config_time_integration) == 'unsplit_explicit') then
            split = 0
         elseif (trim(config_time_integration) == 'split_explicit') then
            split = 1
         endif

         block => domain % blocklist
         do while (associated(block))
           allocate(uTemp(block % mesh % nVertLevels))

           ! Put f*uBcl^{perp} in uNew as a work variable
           call ocn_fuperp(block % state % time_levs(2) % state , block % mesh)

           do iEdge=1,block % mesh % nEdges
              cell1 = block % mesh % cellsOnEdge % array(1,iEdge)
              cell2 = block % mesh % cellsOnEdge % array(2,iEdge)

              uTemp = 0.0  ! could put this after with uTemp(maxleveledgetop+1:nvertlevels)=0
              do k=1,block % mesh % maxLevelEdgeTop % array(iEdge)

                 ! uBclNew = uBclOld + dt*(-f*uBclPerp + T(u*,w*,p*) + g*grad(SSH*) )
                 ! Here uNew is a work variable containing -fEdge(iEdge)*uBclPerp(k,iEdge)
                 uTemp(k) &
                 = block % state % time_levs(1) % state % uBcl % array(k,iEdge) &
                 + dt * (block % tend % u % array (k,iEdge) &
                      + block % state % time_levs(2) % state % u % array (k,iEdge) &  ! this is f*uBcl^{perp}
                      + split*gravity &
                        *(  block % state % time_levs(2) % state % ssh % array(cell2) &
                          - block % state % time_levs(2) % state % ssh % array(cell1) ) &
                          /block % mesh % dcEdge % array(iEdge) )
              enddo

              ! Compute GBtrForcing, the vertically averaged forcing
              sshEdge = 0.5*( &
                  block % state % time_levs(1) % state % ssh % array(cell1) & 
                + block % state % time_levs(1) % state % ssh % array(cell2) ) 

              uhSum = (sshEdge + block % mesh % hZLevel % array(1)) * uTemp(1)
              hSum  =  sshEdge + block % mesh % hZLevel % array(1)

              do k=2,block % mesh % maxLevelEdgeTop % array(iEdge)
                 uhSum = uhSum + block % mesh % hZLevel % array(k) * uTemp(k)
                 hSum  =  hSum + block % mesh % hZLevel % array(k)
              enddo
              block % state % time_levs(1) % state % GBtrForcing % array(iEdge) = split*uhSum/hSum/dt


              do k=1,block % mesh % maxLevelEdgeTop % array(iEdge)
                 ! These two steps are together here:
                 !{\bf u}'_{k,n+1} = {\bf u}'_{k,n} - \Delta t {\overline {\bf G}}
                 !{\bf u}'_{k,n+1/2} = \frac{1}{2}\left({\bf u}^{'}_{k,n} +{\bf u}'_{k,n+1}\right) 
                 ! so that uBclNew is at time n+1/2
                block % state % time_levs(2) % state % uBcl % array(k,iEdge) &
                  = 0.5*( &
                  block % state % time_levs(1) % state % uBcl % array(k,iEdge) &
                  + uTemp(k) - dt * block % state % time_levs(1) % state % GBtrForcing % array(iEdge))
              enddo

           enddo ! iEdge

           deallocate(uTemp)

           block => block % next
         end do

        block => domain % blocklist
        do while (associated(block))
           call mpas_dmpar_exch_halo_field2d_real(domain % dminfo, block % state % time_levs(2) % state % uBcl % array(:,:), &
                                            block % mesh % nVertLevels, block % mesh % nEdges, &
                                            block % parinfo % edgesToSend, block % parinfo % edgesToRecv)

           block => block % next
        end do

      enddo  ! do j=1,config_n_bcl_iter
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END baroclinic iterations on linear Coriolis term
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !  Stage 2: Barotropic velocity (2D) prediction, explicitly subcycled
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      oldBtrSubcycleTime = 1
      newBtrSubcycleTime = 2

      if (trim(config_time_integration) == 'unsplit_explicit') then

         block => domain % blocklist
         do while (associated(block))

            ! For Split_Explicit unsplit, simply set uBtrNew=0, uBtrSubcycle=0, and uNew=uBclNew
            block % state % time_levs(2) % state % uBtr % array(:) = 0.0

            block % state % time_levs(newBtrSubcycleTime) % state % uBtrSubcycle % array(:) = 0.0

               block % state % time_levs(2) % state % u    % array(:,:) &
             = block % state % time_levs(2) % state % uBcl % array(:,:) 

            block => block % next
         end do  ! block

      elseif (trim(config_time_integration) == 'split_explicit') then

         ! Initialize variables for barotropic subcycling
         block => domain % blocklist
         do while (associated(block))

        if (config_filter_btr_mode) then
          block % state % time_levs(1) % state % GBtrForcing % array(:) = 0.0
        endif

            do iCell=1,block % mesh % nCells
              ! sshSubcycleOld = sshOld  
                block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(iCell) & 
              = block % state % time_levs(1) % state % ssh % array(iCell)  

              ! sshNew = sshOld  This is the first for the summation
                block % state % time_levs(2) % state % ssh % array(iCell) & 
              = block % state % time_levs(1) % state % ssh % array(iCell)  
            enddo

           do iEdge=1,block % mesh % nEdges

              ! uBtrSubcycleOld = uBtrOld 
                block % state % time_levs(oldBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge) &
              = block % state % time_levs(1) % state % uBtr % array(iEdge) 

              ! uBtrNew = BtrOld  This is the first for the summation
                block % state % time_levs(2) % state % uBtr % array(iEdge) &
              = block % state % time_levs(1) % state % uBtr % array(iEdge) 

              ! FBtr = 0  
              block % state % time_levs(1) % state % FBtr % array(iEdge) = 0.0
            enddo

            block => block % next
         end do  ! block

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN Barotropic subcycle loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do j=1,config_n_btr_subcycles*config_btr_subcycle_loop_factor

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Barotropic subcycle: initial solve for velecity
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (config_btr_gam1_uWt1>1.0e-12) then  ! only do this part if it is needed in next SSH solve
            uPerpTime = oldBtrSubcycleTime

            block => domain % blocklist
            do while (associated(block))

         do iEdge=1,block % mesh % nEdges

               cell1 = block % mesh % cellsOnEdge % array(1,iEdge)
               cell2 = block % mesh % cellsOnEdge % array(2,iEdge)

            ! Compute -f*uPerp
            uPerp = 0.0
            do i = 1,block % mesh % nEdgesOnEdge % array(iEdge)
               eoe = block % mesh % edgesOnEdge % array(i,iEdge)
               uPerp = uPerp + block % mesh % weightsOnEdge % array(i,iEdge) &
                  * block % state % time_levs(uPerpTime) % state % uBtrSubcycle % array(eoe) &
                  * block % mesh % fEdge  % array(eoe)
            end do

          ! mrp 110606 efficiency note: could make this a 1D integer factor instead of an if statement.
          if (block % mesh % boundaryEdge % array(1,iEdge).eq.1) then
                block % state % time_levs(newBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge) = 0.0
          else

             ! uBtrNew = uBtrOld + dt*(-f*uBtroldPerp - g*grad(SSH) + G)
                block % state % time_levs(newBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge) &
              = block % state % time_levs(oldBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge) &
              + dt/config_n_btr_subcycles *( &
                        uPerp &
                      - gravity &
                        *(  block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell2) &
                          - block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell1) ) &
                          /block % mesh % dcEdge % array(iEdge) &
                      + block % state % time_levs(1) % state % GBtrForcing % array(iEdge) )

          endif

         end do

         !  Implicit solve for barotropic momentum decay
         if ( config_btr_mom_decay) then
            !
            !  Add term to RHS of momentum equation: -1/gamma u
            !
            !  This changes the solve to:
            !  u^{n+1} = u_provis^{n+1}/(1+dt/gamma)
            !
            coef = 1.0/(1.0 + dt/config_n_btr_subcycles/config_btr_mom_decay_time)
            do iEdge=1,block % mesh % nEdges
                block % state % time_levs(newBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge) &
              = block % state % time_levs(newBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge) &
              * coef
            end do

          endif


               block => block % next
            end do  ! block


            !   boundary update on uBtrNew
            block => domain % blocklist
            do while (associated(block))

           call mpas_dmpar_exch_halo_field1d_real(domain % dminfo, &
              block % state % time_levs(newBtrSubcycleTime) % state % uBtrSubcycle % array(:), &
              block % mesh % nEdges, &
              block % parinfo % edgesToSend, block % parinfo % edgesToRecv)

               block => block % next
            end do  ! block

          endif ! config_btr_gam1_uWt1>1.0e-12


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Barotropic subcycle: Compute thickness flux and new SSH: PREDICTOR
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            block => domain % blocklist
            do while (associated(block))

           block % tend % ssh % array(:) = 0.0

           if (config_btr_solve_SSH2) then
             ! If config_btr_solve_SSH2=.true., then do NOT accumulate FBtr in this SSH predictor 
             ! section, because it will be accumulated in the SSH corrector section.
             FBtr_coeff = 0.0
           else
             ! otherwise, DO accumulate FBtr in this SSH predictor section
             FBtr_coeff = 1.0
           endif

           ! config_btr_gam1_uWt1 sets the forward weighting of velocity in the SSH computation
           ! config_btr_gam1_uWt1=  1     flux = uBtrNew*H
           ! config_btr_gam1_uWt1=0.5     flux = 1/2*(uBtrNew+uBtrOld)*H
           ! config_btr_gam1_uWt1=  0     flux = uBtrOld*H
           do iEdge=1,block % mesh % nEdges
              cell1 = block % mesh % cellsOnEdge % array(1,iEdge)
              cell2 = block % mesh % cellsOnEdge % array(2,iEdge)

              sshEdge = 0.5 &
                 *(  block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell1) &
                   + block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell2) )
              hSum = sum(block % mesh % hZLevel % array (1:block % mesh % maxLevelEdgeTop % array(iEdge)))

              flux = ((1.0-config_btr_gam1_uWt1) & 
                        * block % state % time_levs(oldBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge) &
                       + config_btr_gam1_uWt1 &
                        * block % state % time_levs(newBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge)) &
                    * (sshEdge + hSum)

               block % tend % ssh % array(cell1) = block % tend % ssh % array(cell1) &
                 - flux * block % mesh % dvEdge % array(iEdge) 
               block % tend % ssh % array(cell2) = block % tend % ssh % array(cell2) &
                 + flux * block % mesh % dvEdge % array(iEdge) 

               block % state % time_levs(1) % state % FBtr % array(iEdge) &
             = block % state % time_levs(1) % state % FBtr % array(iEdge) &
             + FBtr_coeff*flux
         end do
 
         ! SSHnew = SSHold + dt/J*(-div(Flux))
         do iCell=1,block % mesh % nCells 

                block % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle % array(iCell) & 
              = block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(iCell) & 
              + dt/config_n_btr_subcycles &
                * block % tend % ssh % array(iCell) / block % mesh % areaCell % array (iCell)

         end do

               block => block % next
            end do  ! block

            !   boundary update on SSHnew
            block => domain % blocklist
            do while (associated(block))

!              block % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle % array(:), &

           call mpas_dmpar_exch_halo_field1d_real(domain % dminfo, &
              block % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle % array(:), &
              block % mesh % nCells, &
              block % parinfo % cellsToSend, block % parinfo % cellsToRecv)

               block => block % next
            end do  ! block


! mrp 110801 begin
! This whole section, bounded by 'mrp 110801', may be deleted later if it is found
! that barotropic del2 is not useful.
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Barotropic subcycle: compute btr_divergence and btr_vorticity for del2(u_btr) 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            block => domain % blocklist
            do while (associated(block))
      block % state % time_levs(1) % state % u_diffusionBtr % array(:) = 0.0
      if ( config_btr_mom_eddy_visc2 > 0.0 ) then
      !
      ! Compute circulation and relative vorticity at each vertex
      !
      block % state % time_levs(1) % state % circulationBtr % array(:) = 0.0
      do iEdge=1,block % mesh % nEdges
         vertex1 = block % mesh % verticesOnEdge % array(1,iEdge)
         vertex2 = block % mesh % verticesOnEdge % array(2,iEdge)
             block % state % time_levs(1) % state % circulationBtr % array(vertex1) &
           = block % state % time_levs(1) % state % circulationBtr % array(vertex1) &
           - block % mesh % dcEdge % array (iEdge) &
            *block % state % time_levs(oldBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge)

             block % state % time_levs(1) % state % circulationBtr % array(vertex2) &
           = block % state % time_levs(1) % state % circulationBtr % array(vertex2) &
           + block % mesh % dcEdge % array (iEdge) &
            *block % state % time_levs(oldBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge)
      end do 
      do iVertex=1,block % mesh % nVertices
            block % state % time_levs(1) % state % vorticityBtr % array(iVertex) &
          = block % state % time_levs(1) % state % circulationBtr % array(iVertex) / block % mesh % areaTriangle % array (iVertex)
      end do

      !
      ! Compute the divergence at each cell center
      !
      block % state % time_levs(1) % state % divergenceBtr % array(:) = 0.0
      do iEdge=1,block % mesh % nEdges
         cell1 = block % mesh % cellsOnEdge % array(1,iEdge)
         cell2 = block % mesh % cellsOnEdge % array(2,iEdge)
             block % state % time_levs(1) % state % divergenceBtr % array (cell1) &
           = block % state % time_levs(1) % state % divergenceBtr % array (cell1) &
           + block % state % time_levs(oldBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge) &
            *block % mesh % dvEdge % array(iEdge)

             block % state % time_levs(1) % state % divergenceBtr % array (cell2) &
           = block % state % time_levs(1) % state % divergenceBtr % array (cell2) &
           - block % state % time_levs(oldBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge) &
            *block % mesh % dvEdge % array(iEdge)
      end do
      do iCell = 1,block % mesh % nCells
         block % state % time_levs(1) % state % divergenceBtr % array(iCell) &
       = block % state % time_levs(1) % state % divergenceBtr % array(iCell) &
        /block % mesh % areaCell % array(iCell)
      enddo

      !
      ! Compute Btr diffusion
      !
         do iEdge=1,block % mesh % nEdgesSolve
            cell1 = block % mesh % cellsOnEdge % array(1,iEdge)
            cell2 = block % mesh % cellsOnEdge % array(2,iEdge)
            vertex1 = block % mesh % verticesOnEdge % array(1,iEdge)
            vertex2 = block % mesh % verticesOnEdge % array(2,iEdge)

               ! Here -( vorticityBtr(vertex2) - vorticityBtr(vertex1) ) / dvEdge % array (iEdge)
               ! is - \nabla vorticity pointing from vertex 2 to vertex 1, or equivalently 
               !    + k \times \nabla vorticity pointing from cell1 to cell2.

               block % state % time_levs(1) % state % u_diffusionBtr % array(iEdge) = block % mesh % meshScalingDel2 % array (iEdge) * config_btr_mom_eddy_visc2 * &
                   (( block % state % time_levs(1) % state % divergenceBtr % array(cell2)  - block % state % time_levs(1) % state % divergenceBtr % array(cell1) ) / block % mesh % dcEdge % array (iEdge)  &
                  -( block % state % time_levs(1) % state % vorticityBtr % array(vertex2) - block % state % time_levs(1) % state % vorticityBtr % array(vertex1) ) / block % mesh % dvEdge % array (iEdge))

         end do
      end if
               block => block % next
            end do  ! block
! mrp 110801 end


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Barotropic subcycle: Final solve for velocity.  Iterate for Coriolis term.
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do BtrCorIter=1,config_n_btr_cor_iter

          uPerpTime = newBtrSubcycleTime

          block => domain % blocklist
          do while (associated(block))

         do iEdge=1,block % mesh % nEdges 

               cell1 = block % mesh % cellsOnEdge % array(1,iEdge)
               cell2 = block % mesh % cellsOnEdge % array(2,iEdge)

            ! Compute -f*uPerp
            uPerp = 0.0
            do i = 1,block % mesh % nEdgesOnEdge % array(iEdge)
               eoe = block % mesh % edgesOnEdge % array(i,iEdge)
               uPerp = uPerp + block % mesh % weightsOnEdge % array(i,iEdge) &
                  * block % state % time_levs(uPerpTime) % state % uBtrSubcycle % array(eoe) &
                  * block % mesh % fEdge  % array(eoe) 
            end do

          ! mrp 110606 efficiency note: could make this a 1D integer factor instead of an if statement.
          if (block % mesh % boundaryEdge % array(1,iEdge).eq.1) then
                block % state % time_levs(newBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge) = 0.0
          else

             ! uBtrNew = uBtrOld + dt*(-f*uBtroldPerp - g*grad(SSH) + G)

             sshCell1 = &
               (1-config_btr_gam2_SSHWt1)*block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell1) &
                + config_btr_gam2_SSHWt1 *block % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle % array(cell1)

             sshCell2 = &
               (1-config_btr_gam2_SSHWt1)*block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell2) &
                + config_btr_gam2_SSHWt1 *block % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle % array(cell2)

                block % state % time_levs(newBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge) & 
              = block % state % time_levs(oldBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge) & 
              + dt/config_n_btr_subcycles *( &
                        uPerp &
                      - gravity &
                        *(  sshCell2 &
                          - sshCell1 )&
                          /block % mesh % dcEdge % array(iEdge) &
                      + block % state % time_levs(1) % state % GBtrForcing % array(iEdge) &
                      + block % state % time_levs(1) % state % u_diffusionBtr % array(iEdge))
                      ! added del2 diffusion to btr solve

          endif

         end do

            !  Implicit solve for barotropic momentum decay
         if ( config_btr_mom_decay) then
            !  Add term to RHS of momentum equation: -1/gamma u
            !
            !  This changes the solve to:
            !  u^{n+1} = u_provis^{n+1}/(1+dt/gamma)
            !
            coef = 1.0/(1.0 + dt/config_n_btr_subcycles/config_btr_mom_decay_time)
            do iEdge=1,block % mesh % nEdges 
                block % state % time_levs(newBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge) & 
              = block % state % time_levs(newBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge) &
              * coef
            end do

         endif

               block => block % next
            end do  ! block


            !   boundary update on uBtrNew
            block => domain % blocklist
            do while (associated(block))

               call mpas_dmpar_exch_halo_field1d_real(domain % dminfo, &
                  block % state % time_levs(newBtrSubcycleTime) % state % uBtrSubcycle % array(:), &
                  block % mesh % nEdges, &
                  block % parinfo % edgesToSend, block % parinfo % edgesToRecv)

               block => block % next
            end do  ! block

       end do !do BtrCorIter=1,config_n_btr_cor_iter

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Barotropic subcycle: Compute thickness flux and new SSH: CORRECTOR
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (config_btr_solve_SSH2) then

            block => domain % blocklist
            do while (associated(block))

           block % tend % ssh % array(:) = 0.0

           ! config_btr_gam3_uWt2 sets the forward weighting of velocity in the SSH computation
           ! config_btr_gam3_uWt2=  1     flux = uBtrNew*H
           ! config_btr_gam3_uWt2=0.5     flux = 1/2*(uBtrNew+uBtrOld)*H
           ! config_btr_gam3_uWt2=  0     flux = uBtrOld*H

           do iEdge=1,block % mesh % nEdges
              cell1 = block % mesh % cellsOnEdge % array(1,iEdge)
              cell2 = block % mesh % cellsOnEdge % array(2,iEdge)

              sshEdge = 0.5 &
                 *(  block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell1) &
                   + block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell2) )
              hSum = sum(block % mesh % hZLevel % array (1:block % mesh % maxLevelEdgeTop % array(iEdge)))

              flux = ((1.0-config_btr_gam3_uWt2) & 
                        * block % state % time_levs(oldBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge) &
                       + config_btr_gam3_uWt2 &
                        * block % state % time_levs(newBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge)) &
                    * (sshEdge + hSum)

               block % tend % ssh % array(cell1) = block % tend % ssh % array(cell1) &
                 - flux * block % mesh % dvEdge % array(iEdge) 
               block % tend % ssh % array(cell2) = block % tend % ssh % array(cell2) &
                 + flux * block % mesh % dvEdge % array(iEdge) 

               block % state % time_levs(1) % state % FBtr % array(iEdge) &
             = block % state % time_levs(1) % state % FBtr % array(iEdge) &
             + flux


         end do
 
         ! SSHnew = SSHold + dt/J*(-div(Flux))
         do iCell=1,block % mesh % nCells 

                block % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle % array(iCell) & 
              = block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(iCell) & 
              + dt/config_n_btr_subcycles &
                * block % tend % ssh % array(iCell) / block % mesh % areaCell % array (iCell)

         end do

               block => block % next
            end do  ! block

            !   boundary update on SSHnew
            block => domain % blocklist
            do while (associated(block))

           call mpas_dmpar_exch_halo_field1d_real(domain % dminfo, &
              block % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle % array(:), &
              block % mesh % nCells, &
              block % parinfo % cellsToSend, block % parinfo % cellsToRecv)

               block => block % next
            end do  ! block

        endif ! config_btr_solve_SSH2

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Barotropic subcycle: Accumulate running sums, advance timestep pointers
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            block => domain % blocklist
            do while (associated(block))

         ! Accumulate SSH in running sum over the subcycles.
         do iCell=1,block % mesh % nCells 
                block % state % time_levs(2) % state % ssh % array(iCell) & 
              = block % state % time_levs(2) % state % ssh % array(iCell) & 
              + block % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle % array(iCell)  
         end do

            ! uBtrNew = uBtrNew + uBtrSubcycleNEW
            ! This accumulates the sum.
            ! If the Barotropic Coriolis iteration is limited to one, this could 
            ! be merged with the above code.
         do iEdge=1,block % mesh % nEdges 

                block % state % time_levs(2) % state % uBtr % array(iEdge) & 
              = block % state % time_levs(2) % state % uBtr % array(iEdge) & 
              + block % state % time_levs(newBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge)  

            end do  ! iEdge
               block => block % next
         end do  ! block

            ! advance time pointers
            oldBtrSubcycleTime = mod(oldBtrSubcycleTime,2)+1
            newBtrSubcycleTime = mod(newBtrSubcycleTime,2)+1

         end do ! j=1,config_n_btr_subcycles
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END Barotropic subcycle loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            ! Normalize Barotropic subcycle sums: ssh, uBtr, and F
            block => domain % blocklist
            do while (associated(block))

         do iEdge=1,block % mesh % nEdges
               block % state % time_levs(1) % state % FBtr % array(iEdge) &
             = block % state % time_levs(1) % state % FBtr % array(iEdge) &
             / (config_n_btr_subcycles*config_btr_subcycle_loop_factor)

                block % state % time_levs(2) % state % uBtr % array(iEdge) & 
              = block % state % time_levs(2) % state % uBtr % array(iEdge) & 
             / (config_n_btr_subcycles*config_btr_subcycle_loop_factor + 1)
         end do

        if (config_SSH_from=='avg_of_SSH_subcycles') then
         do iCell=1,block % mesh % nCells 
                block % state % time_levs(2) % state % ssh % array(iCell) & 
              = block % state % time_levs(2) % state % ssh % array(iCell) & 
             / (config_n_btr_subcycles*config_btr_subcycle_loop_factor + 1)
         end do
        elseif (config_SSH_from=='avg_flux') then
           ! see below
        else
         write(0,*) 'Abort: Unknown config_SSH_from option: '&
           //trim(config_SSH_from)
         call mpas_dmpar_abort(dminfo)
        endif

               block => block % next
            end do  ! block


            ! boundary update on F
            block => domain % blocklist
            do while (associated(block))

           call mpas_dmpar_exch_halo_field1d_real(domain % dminfo, &
              block % state % time_levs(1) % state % FBtr % array(:), &
              block % mesh % nEdges, &
              block % parinfo % edgesToSend, block % parinfo % edgesToRecv)

               block => block % next
            end do  ! block


            ! Check that you can compute SSH using the total sum or the individual increments
            ! over the barotropic subcycles.
            ! efficiency: This next block of code is really a check for debugging, and can 
            ! be removed later.
            block => domain % blocklist
            do while (associated(block))

               allocate(uTemp(block % mesh % nVertLevels))

        if (config_SSH_from=='avg_flux') then
           ! Accumulate fluxes in the tend % ssh variable
           block % tend % ssh % array(:) = 0.0
           do iEdge=1,block % mesh % nEdges
              cell1 = block % mesh % cellsOnEdge % array(1,iEdge)
              cell2 = block % mesh % cellsOnEdge % array(2,iEdge)

                 block % tend % ssh % array(cell1) &
               = block % tend % ssh % array(cell1) &
               - block % state % time_levs(1) % state % FBtr % array(iEdge) &
                    * block % mesh % dvEdge % array(iEdge)


                 block % tend % ssh % array(cell2) &
               = block % tend % ssh % array(cell2) &
               + block % state % time_levs(1) % state % FBtr % array(iEdge) &
                    * block % mesh % dvEdge % array(iEdge)

          end do

         do iCell=1,block % mesh % nCells 
 
             ! SSHnew = SSHold + dt*(-div(Flux))
                block % state % time_levs(2) % state % ssh % array(iCell) & 
              = block % state % time_levs(1) % state % ssh % array(iCell) & 
              + dt &
                * block % tend % ssh % array(iCell) / block % mesh % areaCell % array (iCell)
         end do
       endif

         ! Correction velocity    uCorr = (Flux - Sum(h u*))/H
         ! or, for the full latex version:
         !u^{corr} = \left( {\overline {\bf F}} 
         !  - \sum_{k=1}^{N^{edge}} \left(\zeta_{k,n}^{*\;edge}+\Delta z_k\right)  u_k^* \right)
         !\left/ \sum_{k=1}^{N^{edge}} \left(\zeta_{k,n}^{*\;edge}+\Delta z_k\right)   \right. 

          if (config_u_correction) then
             ucorr_coef = 1
          else
             ucorr_coef = 0
          endif

           do iEdge=1,block % mesh % nEdges
              cell1 = block % mesh % cellsOnEdge % array(1,iEdge)
              cell2 = block % mesh % cellsOnEdge % array(2,iEdge)

              sshEdge = 0.5 &
                 *(  block % state % time_levs(2) % state % ssh % array(cell1) &
                   + block % state % time_levs(2) % state % ssh % array(cell2) )

             ! This is u*
               uTemp(:) &
             = block % state % time_levs(2) % state % uBtr % array(iEdge) &
             + block % state % time_levs(2) % state % uBcl % array(:,iEdge) 

              uhSum = (sshEdge + block % mesh % hZLevel % array(1)) * uTemp(1)
              hSum  =  sshEdge + block % mesh % hZLevel % array(1)

              do k=2,block % mesh % maxLevelEdgeTop % array(iEdge)
                 uhSum = uhSum + block % mesh % hZLevel % array(k) * uTemp(k)
                 hSum  =  hSum + block % mesh % hZLevel % array(k)
              enddo

            uCorr =   ucorr_coef*(( block % state % time_levs(1) % state % FBtr % array(iEdge) &
                      - uhSum)/hSum)

              ! put u^{tr}, the velocity for tracer transport, in uNew
          ! mrp 060611 not sure if boundary enforcement is needed here.  
          if (block % mesh % boundaryEdge % array(1,iEdge).eq.1) then
              block % state % time_levs(2) % state % u % array(:,iEdge) = 0.0
          else
            do k=1,block % mesh % maxLevelEdgeTop % array(iEdge)
              block % state % time_levs(2) % state % u % array(k,iEdge) = uTemp(k) + uCorr
            enddo
            do k=block % mesh % maxLevelEdgeTop % array(iEdge)+1,block % mesh % nVertLevels
              block % state % time_levs(2) % state % u % array(k,iEdge) = 0.0
            enddo
          endif

         ! Put new sshEdge values in h_edge array, for the OcnTendScalar call below.
             block % state % time_levs(2) % state % h_edge % array(1,iEdge) &
           = sshEdge + block % mesh % hZLevel % array(1)

           do k=2,block % mesh % nVertLevels
             block % state % time_levs(2) % state % h_edge % array(k,iEdge) &
           = block % mesh % hZLevel % array(k)
           enddo

          end do ! iEdge

         ! Put new SSH values in h array, for the OcnTendScalar call below.
         do iCell=1,block % mesh % nCells 
             block % state % time_levs(2) % state % h % array(1,iCell) &
           = block % state % time_levs(2) % state % ssh % array(iCell) &
           + block % mesh % hZLevel % array(1)

           ! mrp 110601 efficiency note: Since h just moves back and forth between pointers,
           ! this is not necessary once initialized.
           do k=2,block % mesh % nVertLevels
             block % state % time_levs(2) % state % h % array(k,iCell) &
           = block % mesh % hZLevel % array(k)
           enddo
         enddo ! iCell

           deallocate(uTemp)

               block => block % next
            end do  ! block


      endif ! split_explicit

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !  Stage 3: Tracer, density, pressure, vertical velocity prediction
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         block => domain % blocklist
         do while (associated(block))

           call ocn_wtop(block % state % time_levs(2) % state, block % mesh)

      if (trim(config_time_integration) == 'unsplit_explicit') then
           call ocn_tend_h(block % tend, block % state % time_levs(2) % state , block % diagnostics, block % mesh)
      endif

           call ocn_tend_scalar(block % tend, block % state % time_levs(2) % state , block % diagnostics, block % mesh)

           block => block % next
         end do

        ! ---  update halos for prognostic variables

        block => domain % blocklist
        do while (associated(block))
           call mpas_dmpar_exch_halo_field2d_real(domain % dminfo, block % tend % h % array(:,:), &
                                            block % mesh % nVertLevels, block % mesh % nCells, &
                                            block % parinfo % cellsToSend, block % parinfo % cellsToRecv)
           call mpas_dmpar_exch_halo_field3d_real(domain % dminfo, block % tend % tracers % array(:,:,:), &
                                            block % tend % num_tracers, block % mesh % nVertLevels, block % mesh % nCells, &
                                            block % parinfo % cellsToSend, block % parinfo % cellsToRecv)
           block => block % next
        end do


        block => domain % blocklist
        do while (associated(block))
           allocate(hNew(block % mesh % nVertLevels))

        if (trim(config_new_btr_variables_from) == 'last_subcycle') then
           ! This points to the last barotropic SSH subcycle
           sshNew => block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array
        elseif (trim(config_new_btr_variables_from) == 'btr_avg') then
           ! This points to the tendency variable SSH*
           sshNew => block % state % time_levs(2) % state % ssh % array
        endif

      if (trim(config_time_integration) == 'unsplit_explicit') then

         do iCell=1,block % mesh % nCells
           ! this is h_{n+1}
             block % state % time_levs(2) % state % h % array(:,iCell) &
           = block % state % time_levs(1) % state % h % array(:,iCell) &
           + dt* block % tend % h % array(:,iCell) 

            ! this is only for the hNew computation below, so there is the correct
            ! value in the ssh variable for unsplit_explicit case.
            block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(iCell) &
          = block % state % time_levs(2) % state % h % array(1,iCell) &
          - block % mesh % hZLevel % array(1)
           end do ! iCell

      endif ! unsplit_explicit

           ! Only need T & S for earlier iterations,
           ! then all the tracers needed the last time through.
         if (split_explicit_step < config_n_ts_iter) then

           hNew(:) = block % mesh % hZLevel % array(:)
           do iCell=1,block % mesh % nCells
              ! sshNew is a pointer, defined above.
              hNew(1) =  sshNew(iCell) + block % mesh % hZLevel % array(1)
              do k=1,block % mesh % maxLevelCell % array(iCell)
                 do i=1,2
                ! This is Phi at n+1
                tracerTemp &
                = (  block % state % time_levs(1) % state % tracers % array(i,k,iCell) &
                   * block % state % time_levs(1) % state % h % array(k,iCell) &
                 + dt * block % tend % tracers % array(i,k,iCell) &
                  ) / hNew(k)

                ! This is Phi at n+1/2
                   block % state % time_levs(2) % state % tracers % array(i,k,iCell)  &
                 = 0.5*( &
                   block % state % time_levs(1) % state % tracers % array(i,k,iCell) &
                 + tracerTemp )
                 enddo
              end do
           end do ! iCell


          if (trim(config_time_integration) == 'unsplit_explicit') then

            ! compute h*, which is h at n+1/2 and put into array hNew
            ! on last iteration, hNew remains at n+1
           do iCell=1,block % mesh % nCells
                 block % state % time_levs(2) % state % h % array(1,iCell) &
                 = 0.5*( &
                 block % state % time_levs(2) % state % h % array(1,iCell) &
               + block % state % time_levs(1) % state % h % array(1,iCell) )

           end do ! iCell
          endif ! unsplit_explicit

          ! compute u*, the velocity for tendency terms.  Put in uNew.
          ! uBclNew is at time n+1/2 here.
          ! This overwrites u^{tr}, the tracer transport velocity, which was in uNew.
          ! The following must occur after  call OcnTendScalar
           do iEdge=1,block % mesh % nEdges
               block % state % time_levs(2) % state % u    % array(:,iEdge) &
             = block % state % time_levs(2) % state % uBtr % array(iEdge) &
             + block % state % time_levs(2) % state % uBcl % array(:,iEdge) 
           end do ! iEdge

         ! mrp 110512  I really only need this to compute h_edge, density, pressure.
         ! I can par this down later.
         call ocn_diagnostic_solve(dt, block % state % time_levs(2) % state, block % mesh)
 

         elseif (split_explicit_step == config_n_ts_iter) then

           hNew(:) = block % mesh % hZLevel % array(:)
           do iCell=1,block % mesh % nCells
              ! sshNew is a pointer, defined above.
              hNew(1) =  sshNew(iCell) + block % mesh % hZLevel % array(1)
              do k=1,block % mesh % maxLevelCell % array(iCell)
                 do i=1,block % state % time_levs(1) % state % num_tracers
                ! This is Phi at n+1
                   block % state % time_levs(2) % state % tracers % array(i,k,iCell)  &
                = (  block % state % time_levs(1) % state % tracers % array(i,k,iCell) &
                   * block % state % time_levs(1) % state % h % array(k,iCell) &
                 + dt * block % tend % tracers % array(i,k,iCell) &
                  ) / hNew(k)

                 enddo
              end do
           end do

         endif ! split_explicit_step
           deallocate(hNew)

         block => block % next
       end do

      end do  ! split_explicit_step = 1, config_n_ts_iter
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END large iteration loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !
      !  A little clean up at the end: decouple new scalar fields and compute diagnostics for new state
      !
      block => domain % blocklist
      do while (associated(block))

        if (trim(config_new_btr_variables_from) == 'last_subcycle') then
         do iEdge=1,block % mesh % nEdges
               ! uBtrNew = uBtrSubcycleNew  (old here is because counter already flipped)
               ! This line is not needed if u is resplit at the beginning of the timestep.
                block % state % time_levs(2) % state % uBtr % array(iEdge) & 
              = block % state % time_levs(oldBtrSubcycleTime) % state % uBtrSubcycle % array(iEdge)  
         enddo ! iEdges
        elseif (trim(config_new_btr_variables_from) == 'btr_avg') then
               ! uBtrNew from u*.  this is done above, so u* is already in
               ! block % state % time_levs(2) % state % uBtr % array(iEdge) 
        else
         write(0,*) 'Abort: Unknown config_new_btr_variables_from: '&
           //trim(config_time_integration)
         call mpas_dmpar_abort(dminfo)
       endif

         ! Recompute final u to go on to next step.
         ! u_{n+1} = uBtr_{n+1} + uBcl_{n+1} 
         ! Right now uBclNew is at time n+1/2, so back compute to get uBcl at time n+1
         !   using uBcl_{n+1/2} = 1/2*(uBcl_n + u_Bcl_{n+1})
         ! so the following lines are
         ! u_{n+1} = uBtr_{n+1} + 2*uBcl_{n+1/2} - uBcl_n
         ! note that uBcl is recomputed at the beginning of the next timestep due to Imp Vert mixing,
         ! so uBcl does not have to be recomputed here.

         do iEdge=1,block % mesh % nEdges
            do k=1,block % mesh % maxLevelEdgeTop % array(iEdge)

               block % state % time_levs(2) % state % u % array(k,iEdge) & 
            =  block % state % time_levs(2) % state % uBtr % array(iEdge) &
            +2*block % state % time_levs(2) % state % uBcl % array(k,iEdge) &
            -  block % state % time_levs(1) % state % uBcl % array(k,iEdge)
            enddo
            ! mrp 110607 zero out velocity below land edges. efficiency: this may not be required.
            do k=block % mesh % maxLevelEdgeTop % array(iEdge) + 1, block % mesh % nVertLevels
               block % state % time_levs(2) % state % u % array(k,iEdge) = 0.0
            enddo

         enddo ! iEdges

        if (trim(config_time_integration) == 'split_explicit') then

        if (trim(config_new_btr_variables_from) == 'last_subcycle') then
         do iCell=1,block % mesh % nCells
         ! SSH for the next step is from the end of the barotropic subcycle.
               block % state % time_levs(2) % state % ssh % array(iCell) & 
            =  block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(iCell) 
         end do ! iCell
        elseif (trim(config_new_btr_variables_from) == 'btr_avg') then
               ! sshNew from ssh*.  This is done above, so ssh* is already in
               ! block % state % time_levs(2) % state % ssh % array(iCell) 
        endif

         do iCell=1,block % mesh % nCells
         ! Put new SSH values in h array, for the OcnTendScalar call below.
             block % state % time_levs(2) % state % h % array(1,iCell) &
           = block % state % time_levs(2) % state % ssh % array(iCell) &
           + block % mesh % hZLevel % array(1)

           ! mrp 110601 efficiency note: Since h just moves back and forth between pointers,
           ! this is not necessary once initialized.
           do k=2,block % mesh % nVertLevels
             block % state % time_levs(2) % state % h % array(k,iCell) &
           = block % mesh % hZLevel % array(k)
           end do
         end do ! iCell
       end if ! split_explicit
 
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       !  Implicit vertical mixing, done after timestep is complete
       !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         u           => block % state % time_levs(2) % state % u % array
         tracers     => block % state % time_levs(2) % state % tracers % array
         h           => block % state % time_levs(2) % state % h % array
         h_edge      => block % state % time_levs(2) % state % h_edge % array
         ke_edge     => block % state % time_levs(2) % state % ke_edge % array
         num_tracers = block % state % time_levs(2) % state % num_tracers
         vertViscTopOfEdge => block % diagnostics % vertViscTopOfEdge % array
         vertDiffTopOfCell => block % diagnostics % vertDiffTopOfCell % array
         maxLevelCell    => block % mesh % maxLevelCell % array
         maxLevelEdgeTop => block % mesh % maxLevelEdgeTop % array

         if (config_implicit_vertical_mix) then
            allocate(A(block % mesh % nVertLevels),C(block % mesh % nVertLevels),uTemp(block % mesh % nVertLevels), &
               tracersTemp(num_tracers,block % mesh % nVertLevels))

            call ocn_vmix_coefs(block % mesh, block % state % time_levs(2) % state, block % diagnostics, err)

            !
            !  Implicit vertical solve for momentum
            !

            call ocn_vel_vmix_tend_implicit(block % mesh, dt, ke_edge, vertvisctopofedge, h, h_edge, u, err)

            !
            !  Implicit vertical solve for tracers
            !
            call ocn_tracer_vmix_tend_implicit(block % mesh, dt, vertdifftopofcell, h, tracers, err)
         end if

         ! mrp 110725 adding momentum decay term
         if (config_mom_decay) then

            !
            !  Implicit solve for momentum decay
            !
            !  Add term to RHS of momentum equation: -1/gamma u
            !
            !  This changes the solve to:
            !  u^{n+1} = u_provis^{n+1}/(1+dt/gamma)
            !
            coef = 1.0/(1.0 + dt/config_mom_decay_time)
            do iEdge=1,block % mesh % nEdges
               do k=1,maxLevelEdgeTop(iEdge)
                  u(k,iEdge) = coef*u(k,iEdge) 
               end do
            end do

         end if

         if (config_test_case == 1) then    ! For case 1, wind field should be fixed
            block % state % time_levs(2) % state % u % array(:,:) = block % state % time_levs(1) % state % u % array(:,:)
         end if

         call ocn_diagnostic_solve(dt, block % state % time_levs(2) % state, block % mesh)

         call mpas_reconstruct(block % mesh, block % state % time_levs(2) % state % u % array,          &
                          block % state % time_levs(2) % state % uReconstructX % array,            &
                          block % state % time_levs(2) % state % uReconstructY % array,            &
                          block % state % time_levs(2) % state % uReconstructZ % array,            &
                          block % state % time_levs(2) % state % uReconstructZonal % array,        &
                          block % state % time_levs(2) % state % uReconstructMeridional % array    &
                         )

         block => block % next
      end do
      call mpas_timer_stop("split_explicit_timestep")

   end subroutine ocn_time_integrator_split!}}}

   subroutine filter_btr_mode_tend_u(tend, s, d, grid)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Filter and remove barotropic mode from the tendencies
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

! mrp 110512 I just split compute_tend into compute_tend_u and compute_tend_h.
!  Some of these variables can be removed, but at a later time.
      integer :: iEdge, iCell, iVertex, k, cell1, cell2, &
        vertex1, vertex2, eoe, i, j

      integer :: nCells, nEdges, nVertices, nVertLevels, nEdgesSolve
      real (kind=RKIND) :: vertSum, uhSum, hSum, sshEdge
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

      call mpas_timer_start("filter_btr_mode_tend_u")

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

           do iEdge=1,grid % nEdges

              ! I am using hZLevel here.  This assumes that SSH is zero everywhere already,
              ! which should be the case if the barotropic mode is filtered.
              ! The more general case is to use sshEdge or h_edge.
              uhSum = (grid % hZLevel % array(1)) * tend_u(1,iEdge)
              hSum  =  grid % hZLevel % array(1)

              do k=2,grid % maxLevelEdgeTop % array(iEdge)
                 uhSum = uhSum + grid % hZLevel % array(k) * tend_u(k,iEdge)
                 hSum  =  hSum + grid % hZLevel % array(k)
              enddo

              vertSum = uhSum/hSum

              do k=1,grid % maxLevelEdgeTop % array(iEdge)
                 tend_u(k,iEdge) = tend_u(k,iEdge) - vertSum
              enddo

           enddo ! iEdge

      call mpas_timer_stop("filter_btr_mode_tend_u")

   end subroutine filter_btr_mode_tend_u!}}}

   subroutine filter_btr_mode_u(s, grid)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Filter and remove barotropic mode.
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
      real (kind=RKIND) :: vertSum, uhSum, hSum, sshEdge
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

      call mpas_timer_start("filter_btr_mode_u")

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
! mrp 110516 cleanup fvertex fedge not used in this subroutine
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

      u_src => grid % u_src % array

           do iEdge=1,grid % nEdges

              ! I am using hZLevel here.  This assumes that SSH is zero everywhere already,
              ! which should be the case if the barotropic mode is filtered.
              ! The more general case is to use sshedge or h_edge.
              uhSum = (grid % hZLevel % array(1)) * u(1,iEdge)
              hSum  =  grid % hZLevel % array(1)

              do k=2,grid % maxLevelEdgeTop % array(iEdge)
                 uhSum = uhSum + grid % hZLevel % array(k) * u(k,iEdge)
                 hSum  =  hSum + grid % hZLevel % array(k)
              enddo

              vertSum = uhSum/hSum
              do k=1,grid % maxLevelEdgeTop % array(iEdge)
                 u(k,iEdge) = u(k,iEdge) - vertSum
              enddo

           enddo ! iEdge

      call mpas_timer_stop("filter_btr_mode_u")

   end subroutine filter_btr_mode_u!}}}

   subroutine enforce_boundaryEdge(tend, grid)!{{{
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

      call mpas_timer_start("enforce_boundaryEdge")

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
      call mpas_timer_stop("enforce_boundaryEdge")

   end subroutine enforce_boundaryEdge!}}}

end module ocn_time_integration_split

! vim: foldmethod=marker
