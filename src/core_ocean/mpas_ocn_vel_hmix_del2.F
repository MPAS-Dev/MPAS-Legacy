!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_hmix_del2
!
!> \brief Ocean horizontal mixing - Laplacian parameterization 
!> \author Phil Jones, Doug Jacobsen
!> \date   15 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains routines for computing horizontal mixing 
!>  tendencies using a Laplacian formulation.
!
!-----------------------------------------------------------------------

module ocn_vel_hmix_del2

   use mpas_grid_types
   use mpas_configure

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

   public :: ocn_vel_hmix_del2_tend, &
             ocn_vel_hmix_del2_init

   !-------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical ::  hmixDel2On  !< integer flag to determine whether del2 chosen

   real (kind=RKIND) :: &
      eddyVisc2,        &!< base eddy diffusivity for Laplacian
      viscVortCoef

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_hmix_del2_tend
!
!> \brief   Computes tendency term for Laplacian horizontal momentum mixing
!> \author  Phil Jones, Doug Jacobsen
!> \date    22 August 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the horizontal mixing tendency for momentum
!>  based on a Laplacian form for the mixing, \f$\nu_2 \nabla^2 u\f$
!>  This tendency takes the
!>  form \f$\nu( \nabla divergence + k \times \nabla relativeVorticity )\f$,
!>  where \f$\nu\f$ is a viscosity and \f$k\f$ is the vertical unit vector.
!>  This form is strictly only valid for constant \f$\nu\f$ .
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hmix_del2_tend(grid, divergence, relativeVorticity, viscosity, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         divergence      !< Input: velocity divergence

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         relativeVorticity       !< Input: relative vorticity

      type (mesh_type), intent(in) :: &
         grid            !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend             !< Input/Output: velocity tendency

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         viscosity       !< Input: viscosity

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

      integer :: iEdge, nEdgesSolve, cell1, cell2, vertex1, vertex2, k
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, edgeMask

      real (kind=RKIND) :: u_diffusion, invLength1, invLength2, visc2
      real (kind=RKIND), dimension(:), pointer :: meshScalingDel2, &
              dcEdge, dvEdge

      !-----------------------------------------------------------------
      !
      ! exit if this mixing is not selected
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.hmixDel2On) return

      nEdgesSolve = grid % nEdgesSolve
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      cellsOnEdge => grid % cellsOnEdge % array
      verticesOnEdge => grid % verticesOnEdge % array
      meshScalingDel2 => grid % meshScalingDel2 % array
      edgeMask => grid % edgeMask % array
      dcEdge => grid % dcEdge % array
      dvEdge => grid % dvEdge % array

      do iEdge=1,nEdgesSolve
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)

         invLength1 = 1.0 / dcEdge(iEdge)
         invLength2 = 1.0 / dvEdge(iEdge)

         do k=1,maxLevelEdgeTop(iEdge)

            ! Here -( relativeVorticity(k,vertex2) - relativeVorticity(k,vertex1) ) / dvEdge(iEdge)
            ! is - \nabla relativeVorticity pointing from vertex 2 to vertex 1, or equivalently 
            !    + k \times \nabla relativeVorticity pointing from cell1 to cell2.

            u_diffusion = ( divergence(k,cell2)  - divergence(k,cell1) ) * invLength1 &
                          -viscVortCoef &
                          *( relativeVorticity(k,vertex2) - relativeVorticity(k,vertex1) ) * invLength2

            visc2 = meshScalingDel2(iEdge) * eddyVisc2

            tend(k,iEdge) = tend(k,iEdge) + edgeMask(k, iEdge) * visc2 * u_diffusion

            viscosity(k,iEdge) = viscosity(k,iEdge) + visc2

         end do
      end do

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del2_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_hmix_del2_init
!
!> \brief   Initializes ocean momentum Laplacian horizontal mixing
!> \author  Phil Jones, Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  Laplacian horizontal momentum mixing in the ocean.  
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hmix_del2_init(err)!{{{


   integer, intent(out) :: err !< Output: error flag

   !--------------------------------------------------------------------
   !
   ! set some local module variables based on input config choices
   !
   !--------------------------------------------------------------------

   err = 0

   hmixDel2On = .false.

   if ( config_mom_del2 > 0.0 ) then
      hmixDel2On = .true.
      eddyVisc2 = config_mom_del2

      if (config_visc_vorticity_term) then
         viscVortCoef = config_vorticity_del2_scale
      else
         viscVortCoef = 0.0
      endif
   endif

   if(.not.config_use_mom_del2) hmixDel2On = .false.


   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del2_init!}}}

!***********************************************************************

end module ocn_vel_hmix_del2

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
