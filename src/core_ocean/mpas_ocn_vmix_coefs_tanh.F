!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vmix_coefs_tanh
!
!> \brief MPAS ocean vertical mixing coefficients
!> \author Doug Jacobsen
!> \date   19 September 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the routines for computing 
!>  tanhant vertical mixing coefficients.  
!>
!
!-----------------------------------------------------------------------

module ocn_vmix_coefs_tanh

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

   public :: ocn_vmix_coefs_tanh_build, &
             ocn_vmix_coefs_tanh_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: tanhViscOn, tanhDiffOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vmix_coefs_tanh_build
!
!> \brief   Computes coefficients for vertical mixing
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the vertical mixing coefficients for momentum
!>  and tracers based user choices of mixing parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_coefs_tanh_build(grid, s, d, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (state_type), intent(inout) :: &
         s             !< Input/Output: state information

      type (diagnostics_type), intent(inout) :: &
         d             !< Input/Output: diagnostic information

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

      integer :: err1, err2

      real (kind=RKIND), dimension(:,:), pointer :: &
        vertViscTopOfEdge, vertDiffTopOfCell

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0
      if((.not.tanhViscOn) .and. (.not.tanhDiffOn)) return

      vertViscTopOfEdge => d % vertViscTopOfEdge % array
      vertDiffTopOfCell => d % vertDiffTopOfCell % array

      call ocn_vel_vmix_coefs_tanh(grid, vertViscTopOfEdge, err1)
      call ocn_tracer_vmix_coefs_tanh(grid, vertDiffTopOfCell, err2)

      err = ior(err1, err2)

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs_tanh_build!}}}

!***********************************************************************
!
!  routine ocn_vel_vmix_coefs_tanh
!
!> \brief   Computes coefficients for vertical momentum mixing
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the tanh vertical mixing coefficients for momentum
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vmix_coefs_tanh(grid, vertViscTopOfEdge, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: vertViscTopOfEdge !< Output: Vertical viscosity

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: k, nVertLevels

      real (kind=RKIND), dimension(:), pointer :: zTopZLevel

      err = 0

      if(.not.tanhViscOn) return

      nVertLevels = grid % nVertLevels
      zTopZLevel => grid % zTopZLevel % array

      do k=1,nVertLevels+1
          vertViscTopOfEdge(k,:) = -(config_max_visc_tanh-config_min_visc_tanh)/2.0 &
            *tanh(-(zTopZLevel(k)-config_ZMid_tanh) &
                  /config_zWidth_tanh) &
            + (config_max_visc_tanh+config_min_visc_tanh)/2
      end do


   !--------------------------------------------------------------------

   end subroutine ocn_vel_vmix_coefs_tanh!}}}

!***********************************************************************
!
!  routine ocn_tracer_vmix_coefs_tanh
!
!> \brief   Computes coefficients for vertical tracer mixing
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the tanh vertical mixing coefficients for tracers
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_vmix_coefs_tanh(grid, vertDiffTopOfCell, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: vertDiffTopOfCell !< Output: Vertical diffusion

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: k, nVertLevels

      real (kind=RKIND), dimension(:), pointer :: zTopZLevel

      err = 0

      if(.not.tanhDiffOn) return

      nVertLevels = grid % nVertLevels
      zTopZLevel => grid % zTopZLevel % array

      do k=1,nVertLevels+1
         vertDiffTopOfCell(k,:) = -(config_max_diff_tanh-config_min_diff_tanh)/2.0 &
            *tanh(-(zTopZLevel(k)-config_ZMid_tanh) &
                  /config_zWidth_tanh) &
            + (config_max_diff_tanh+config_min_diff_tanh)/2
      end do


   !--------------------------------------------------------------------

   end subroutine ocn_tracer_vmix_coefs_tanh!}}}


!***********************************************************************
!
!  routine ocn_vmix_coefs_tanh_init
!
!> \brief   Initializes ocean vertical mixing quantities
!> \author  Doug Jacobsen
!> \date    19 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  tanh vertical mixing in the ocean. 
!
!-----------------------------------------------------------------------


   subroutine ocn_vmix_coefs_tanh_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      tanhViscOn = .false.
      tanhDiffOn = .false.

      if (config_vert_visc_type.eq.'tanh') then
          tanhViscOn = .true.
      endif

      if (config_vert_diff_type.eq.'tanh') then
          tanhDiffOn = .true.
      endif

      if(tanhViscOn .or. tanhDiffOn) then
         if (config_vert_grid_type.ne.'zlevel') then
            write(0,*) 'Abort: config_vert_diff_type.eq.tanh may only', &
                       ' use config_vert_grid_type of zlevel at this time'
            err = 1
         endif
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs_tanh_init!}}}

!***********************************************************************

end module ocn_vmix_coefs_tanh

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
