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
!>  tendency terms for the ocean primitive equations.
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

   use ocn_vel_coriolis
   use ocn_vel_pressure_grad
   use ocn_vel_vadv
   use ocn_vel_hmix
   use ocn_vel_forcing
   use ocn_vmix

   use ocn_tracer_hadv
   use ocn_tracer_vadv
   use ocn_tracer_hmix
   use ocn_restoring

   implicit none
   private
   save

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
             ocn_tend_tracer, &
             ocn_tendency_init

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

      if(config_disable_h_all_tend) return

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
        tend_u, circulation, vorticity, viscosity, ke, ke_edge, Vor_edge, &
        MontPot, wTop, divergence, vertViscTopOfEdge

      real (kind=RKIND), dimension(:,:), pointer :: u_src

      integer :: err

      call mpas_timer_start("ocn_tend_u")

      u           => s % u % array
      rho         => s % rho % array
      wTop        => s % wTop % array
      zMid        => s % zMid % array
      h_edge      => s % h_edge % array
      viscosity   => s % viscosity % array
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

      if(config_disable_u_all_tend) return

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
      if (config_pressure_gradient_type.eq.'MontgomeryPotential') then
          call ocn_vel_pressure_grad_tend(grid, MontPot,  zMid, rho, tend_u, err)
      else
          call ocn_vel_pressure_grad_tend(grid, pressure, zMid, rho, tend_u, err)
      end if
      call mpas_timer_stop("pressure grad", velPgradTimer)

      !
      ! velocity tendency: del2 dissipation, \nu_2 \nabla^2 u
      !   computed as \nu( \nabla divergence + k \times \nabla vorticity )
      !   strictly only valid for config_mom_del2 == constant
      !
      call mpas_timer_start("hmix", .false., velHmixTimer)
      call ocn_vel_hmix_tend(grid, divergence, vorticity, viscosity, tend_u, err)
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
!  routine ocn_tend_tracer
!
!> \brief   Computes tracer tendency
!> \author  Doug Jacobsen
!> \date    23 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes tracer tendencies for the ocean
!
!-----------------------------------------------------------------------
   subroutine ocn_tend_tracer(tend, s, d, grid, dt)!{{{
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

      call mpas_timer_start("ocn_tend_tracer")

      uTransport  => s % uTransport % array
      h           => s % h % array
      wTop        => s % wTop % array
      tracers     => s % tracers % array
      h_edge      => s % h_edge % array
      vertDiffTopOfCell => d % vertDiffTopOfCell % array

      tend_tr     => tend % tracers % array
      tend_h      => tend % h % array

      !
      ! initialize tracer tendency (RHS of tracer equation) to zero.
      !
      tend_tr(:,:,:) = 0.0

      if(config_disable_tr_all_tend) return

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
      call mpas_timer_stop("ocn_tend_tracer")

      deallocate(uh)

   end subroutine ocn_tend_tracer!}}}

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

    end subroutine ocn_tendency_init!}}}

!***********************************************************************

end module ocn_tendency

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
