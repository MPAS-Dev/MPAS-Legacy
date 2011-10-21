!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_hmix_del4
!
!> \brief Ocean horizontal mixing - biharmonic parameterization
!> \author Phil Jones, Doug Jacobsen
!> \date   15 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains routines and variables for computing 
!>  horizontal mixing tendencies using a biharmonic formulation. 
!
!-----------------------------------------------------------------------

module ocn_vel_hmix_del4

   use mpas_grid_types
   use mpas_configure
   use mpas_timer

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

   public :: ocn_vel_hmix_del4_tend, &
             ocn_vel_hmix_del4_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: &
      hmixDel4On       !< local flag to determine whether del4 chosen

   real (kind=RKIND) :: &
      eddyVisc4,        &!< base eddy diffusivity for biharmonic
      viscVortCoef


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_hmix_del4_tend
!
!> \brief   Computes tendency term for biharmonic horizontal momentum mixing
!> \author  Phil Jones, Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the horizontal mixing tendency for momentum
!>  based on a biharmonic form for the mixing.  This mixing tendency
!>  takes the form  \f$-\nu_4 \nabla^4 u\f$
!>  but is computed as 
!>  \f$\nabla^2 u = \nabla divergence + k \times \nabla vorticity\f$
!>  applied recursively.
!>  This formulation is only valid for constant \f$\nu_4\f$ .
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hmix_del4_tend(grid, divergence, vorticity, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         divergence      !< Input: velocity divergence

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vorticity       !< Input: vorticity

      type (mesh_type), intent(in) :: &
         grid           !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend       !< Input/Output: velocity tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------
    
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, nEdges, cell1, cell2, vertex1, vertex2, k
      integer :: iCell, iVertex
      integer :: nVertices, nVertLevels, nCells

      integer, dimension(:), pointer :: maxLevelEdgeTop, maxLevelVertexBot, &
            maxLevelCell
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge


      real (kind=RKIND) :: u_diffusion, r
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaTriangle, &
            meshScalingDel4, areaCell

      real (kind=RKIND), dimension(:,:), allocatable :: delsq_divergence, &
            delsq_u, delsq_circulation, delsq_vorticity

      err = 0

      if(.not.hmixDel4On) return

      call mpas_timer_start("compute_tend-horiz mix-del4")

      nCells = grid % nCells
      nEdges = grid % nEdges
      nVertices = grid % nVertices
      nVertLevels = grid % nVertLevels
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      maxLevelVertexBot => grid % maxLevelVertexBot % array
      maxLevelCell => grid % maxLevelCell % array
      cellsOnEdge => grid % cellsOnEdge % array
      verticesOnEdge => grid % verticesOnEdge % array
      dcEdge => grid % dcEdge % array
      dvEdge => grid % dvEdge % array
      areaTriangle => grid % areaTriangle % array
      areaCell => grid % areaCell % array
      meshScalingDel4 => grid % meshScalingDel4 % array

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
               -viscVortCoef &
               *( vorticity(k,vertex2) - vorticity(k,vertex1)) / dvEdge(iEdge)

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
            delsq_u(k,iEdge) = & 
               ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
              -( vorticity(k,vertex2) - vorticity(k,vertex1)) / dvEdge(iEdge)

            u_diffusion = (  delsq_divergence(k,cell2) &
                           - delsq_divergence(k,cell1) ) / dcEdge(iEdge)  &
                         -viscVortCoef &
                         *(  delsq_vorticity(k,vertex2) &
                           - delsq_vorticity(k,vertex1) ) / dvEdge(iEdge)

            u_diffusion = meshScalingDel4(iEdge) * eddyVisc4 * u_diffusion

            tend(k,iEdge) = tend(k,iEdge) - u_diffusion
         end do
      end do

      deallocate(delsq_divergence)
      deallocate(delsq_u)
      deallocate(delsq_circulation)
      deallocate(delsq_vorticity)

      call mpas_timer_stop("compute_tend-horiz mix-del4")

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del4_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_hmix_del4_init
!
!> \brief   Initializes ocean momentum biharmonic horizontal mixing
!> \author  Phil Jones, Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  biharmonic horizontal tracer mixing in the ocean.  
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hmix_del4_init(err)!{{{

   integer, intent(out) :: err !< Output: error flag

   !--------------------------------------------------------------------
   !
   ! set some local module variables based on input config choices
   !
   !--------------------------------------------------------------------

   err = 0

   hmixDel4On = .false.

   if ( config_h_mom_eddy_visc4 > 0.0 ) then
      hmixDel4On = .true.
      eddyVisc4 = config_h_mom_eddy_visc4
      if (config_visc_vorticity_term) then
         viscVortCoef = 1.0
      else
         viscVortCoef = 0.0
      endif

   endif

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del4_init!}}}

!***********************************************************************

end module ocn_vel_hmix_del4

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
