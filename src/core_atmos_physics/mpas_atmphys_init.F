!=============================================================================================
 module mpas_atmphys_init
 use mpas_grid_types
 use mpas_configure, only: config_do_restart,          &
                           config_lsm_scheme,          &
                           config_microp_scheme,       &
                           config_conv_deep_scheme,    &
                           config_sfclayer_scheme,     &
                           config_radt_lw_scheme,      &
                           config_radt_sw_scheme
 use mpas_timekeeping

 use mpas_atmphys_driver_convection_deep, only: init_convection_deep
 use mpas_atmphys_driver_lsm,only: init_lsm
 use mpas_atmphys_driver_microphysics
 use mpas_atmphys_driver_radiation_lw, only: init_radiation_lw
 use mpas_atmphys_driver_radiation_sw, only: init_radiation_sw
 use mpas_atmphys_driver_sfclayer

 use mpas_atmphys_landuse

 implicit none
 private
 public:: physics_init

 contains

!=============================================================================================
 subroutine physics_init(dminfo,clock,config_do_restart,mesh,state_1,state_2, &
                         diag_physics,sfc_input)
!=============================================================================================

!input arguments:
 logical,intent(in):: config_do_restart
 type (dm_info), intent(in):: dminfo
 type(mesh_type),intent(in):: mesh
 type(MPAS_Clock_type),intent(in):: clock

!inout arguments:
 type(state_type),intent(inout):: state_1,state_2
 type(diag_physics_type),intent(inout):: diag_physics
 type(sfc_input_type) ,intent(inout):: sfc_input

!local variables:
 type(MPAS_Time_Type):: currTime

 logical:: init_done
 integer:: iCell,iLag,ierr,julday 

!---------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine physics_init:'

 currTime = mpas_get_clock_time(clock,MPAS_NOW,ierr)
 call mpas_get_time(curr_time=currTime,DoY=julday,ierr=ierr)

!initialization of east-north directions to convert u-tendencies from cell centers to cell
!edges:
 call init_dirs_forphys(mesh)

!initialization of temperatures needed for updating the deep soil temperature:
 do iCell = 1, mesh % nCellsSolve
    diag_physics % nsteps_accum % array(iCell) = 0
    diag_physics % ndays_accum  % array(iCell) = 1

    diag_physics % tday_accum  % array(iCell) = sfc_input % tmn % array(iCell)
    diag_physics % tyear_mean  % array(iCell) = sfc_input % tmn % array(iCell)
    diag_physics % tyear_accum % array(iCell) = sfc_input % tmn % array(iCell)
    do iLag = 1, mesh % nLags
       diag_physics % tlag % array(iLag,iCell) = sfc_input % tmn % array(iCell)
    enddo
 enddo

!initialization of global surface properties. set here for now, but may be moved when time
!manager is implemented:
 if(.not. config_do_restart) &
    call landuse_init_forMPAS(dminfo,julday,mesh,diag_physics,sfc_input)

!initialization of parameterized deep convective processes:
 if(config_conv_deep_scheme .ne. 'off') &
    call init_convection_deep(config_do_restart,mesh,diag_physics)

!initialization of cloud microphysics processes:
 if(config_microp_scheme .ne. 'off') call microphysics_init

!initialization of surface layer processes:
 if(config_sfclayer_scheme .ne. 'off') call init_sfclayer

!initialization of land-surface model:
!if(.not. config_do_restart) then
!   if(config_lsm_scheme .ne. 'off') call init_lsm(dminfo,mesh,diag_physics,sfc_input)
!endif
 if(config_lsm_scheme .ne. 'off') call init_lsm(dminfo,mesh,diag_physics,sfc_input)

!initialization of shortwave radiation processes:
 init_done = .false.
 if(config_radt_sw_scheme.ne.'off') then
    if(trim(config_radt_sw_scheme) .eq. 'cam_sw') then
       call init_radiation_sw(dminfo,mesh,state_1,state_2)
       init_done = .true.
    else
       call init_radiation_sw(dminfo)
    endif
 endif

!initialization of longwave radiation processes: if we run the CAM radiation codes, the initia
!lization of the longwave and shortwave parameterizations is the same, and needs to be called
!only once:
 if(config_radt_lw_scheme.ne.'off') then    
    if(trim(config_radt_lw_scheme) .eq. 'cam_lw') then
       if(.not. init_done) then
          call init_radiation_lw(dminfo,mesh,state_1,state_2)
       else
          write(0,*)
          write(0,*) '--- camrad lw initialization done above'
       endif
    else
       call init_radiation_lw(dminfo)
    endif
 endif

 write(0,*)
 write(0,*) '--- end subroutine physics_init'
 write(0,*)

 end subroutine physics_init

!=============================================================================================
 subroutine init_dirs_forphys(mesh)
!=============================================================================================

!inout arguments:
!----------------
 type(mesh_type),intent(in):: mesh

!local variables:
 integer:: iCell
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:,:),pointer:: east,north

!---------------------------------------------------------------------------------------------

 latCell => mesh % latCell % array
 lonCell => mesh % lonCell % array
 east    => mesh % east  % array
 north   => mesh % north % array

!Compute unit vectors in east and north directions for each cell:
 do iCell = 1, mesh % nCells

    east(1,iCell) = -sin(lonCell(iCell))
    east(2,iCell) =  cos(lonCell(iCell))
    east(3,iCell) =  0.0
    call r3_normalize(east(1,iCell), east(2,iCell), east(3,iCell))

    north(1,iCell) = -cos(lonCell(iCell))*sin(latCell(iCell))
    north(2,iCell) = -sin(lonCell(iCell))*sin(latCell(iCell))
    north(3,iCell) =  cos(latCell(iCell))
    call r3_normalize(north(1,iCell), north(2,iCell), north(3,iCell))

 end do

 end subroutine init_dirs_forphys

!=============================================================================================
 subroutine r3_normalize(ax, ay, az)
!=============================================================================================
!normalizes the vector (ax, ay, az)

 real (kind=RKIND), intent(inout) :: ax, ay, az
 real (kind=RKIND) :: mi

!---------------------------------------------------------------------------------------------

 mi = 1.0 / sqrt(ax**2 + ay**2 + az**2)
 ax = ax * mi
 ay = ay * mi
 az = az * mi

 end subroutine r3_normalize 
!=============================================================================================
 end module mpas_atmphys_init
!=============================================================================================
