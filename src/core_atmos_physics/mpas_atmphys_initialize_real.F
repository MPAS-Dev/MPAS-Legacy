!=============================================================================================
 module mpas_atmphys_initialize_real
 use mpas_kind_types
 use mpas_configure, only: config_met_prefix,  &
                           config_frac_seaice, &
                           config_input_sst,   &
                           config_nsoillevels, &
                           config_start_time,  &
                           config_sst_prefix
 use mpas_grid_types
 use init_atm_hinterp
 use init_atm_llxy
 use init_atm_read_met

 use mpas_atmphys_date_time
 use mpas_atmphys_utilities

 implicit none
 private
 public:: physics_initialize_real

 contains

!=============================================================================================
 subroutine physics_initialize_sst(mesh,fg)
!=============================================================================================

!input arguments:
 type(mesh_type),intent(in) :: mesh

!inout arguments:
 type(fg_type),intent(inout):: fg 

!local variables:
 character(len=32):: timeString
 integer:: i,j,iCell,istatus
 integer,dimension(5) :: interp_list

 type(met_data) :: field
 type(proj_info):: proj

 real(kind=RKIND):: lat, lon, x, y
 real(kind=RKIND),allocatable,dimension(:,:):: slab_r8

!---------------------------------------------------------------------------------------------

 write(0,*) '--- enter subroutine physics_initialize_sst:'

!set interpolation sequence to be used for SST/SEAICE field:
 interp_list(1) = FOUR_POINT
 interp_list(2) = SEARCH
 interp_list(3) = 0

!open intermediate file:
 call read_met_init(trim(config_sst_prefix),.false.,config_start_time(1:13),istatus)
 if(istatus /= 0) &
    write(0,*) 'Error reading ',trim(config_sst_prefix)//':'//config_start_time(1:13)
 write(0,*) 'Processing ',trim(config_sst_prefix)//':'//config_start_time(1:13)

!scan through all the fields in the file:
 call read_next_met_field(field,istatus)
 do while (istatus == 0)
 
    write(0,*) field % field
    if(index(field % field,'SST'   ) /= 0 .or. &
       index(field % field,'SEAICE') /= 0 .or. &
       index(field % field,'ALBEDO') /= 0 .or. &
       index(field % field,'VEGFRA') /= 0 ) then
       write(0,*) field % field

       !Interpolation routines use real(kind=RKIND), so copy from default real array
       allocate(slab_r8(field % nx, field % ny))
       do j=1,field % ny
       do i=1,field % nx 
          slab_r8(i,j) = field % slab(i,j)
       end do
       end do

       !
       !Set up map projection
       !
       call map_init(proj)
             
       if(field % iproj == PROJ_LATLON) then
          call map_set(PROJ_LATLON, proj, &
                       latinc = real(field % deltalat,RKIND), &
                       loninc = real(field % deltalon,RKIND), &
                       knowni = 1.0_RKIND, &
                       knownj = 1.0_RKIND, &
                       lat1 = real(field % startlat,RKIND), &
                       lon1 = real(field % startlon,RKIND))
       end if
   
       !Interpolate field to each MPAS grid cell:
       do iCell=1,mesh % nCells
          lat = mesh % latCell % array(iCell) * DEG_PER_RAD
          lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
          call latlon_to_ij(proj, lat, lon, x, y)
          if (y < 0.5) then
              y = 1.0
          else if (y >= real(field%ny)+0.5) then
              y = real(field % ny)
          endif 
          if (x < 0.5) then
              lon = lon + 360.0
              call latlon_to_ij(proj, lat, lon, x, y)
          else if (x >= real(field%nx)+0.5) then
              lon = lon - 360.0
              call latlon_to_ij(proj, lat, lon, x, y)
          end if

          if(index(field % field,'SST') /= 0) then
             fg % sst % array(iCell) = interp_sequence(x,y,1,slab_r8,1,field%nx, &
                                              1,field%ny,1,1,-1.e30_RKIND,interp_list,1)
          elseif(index(field % field,'SEAICE') /= 0) then
             fg % xice % array(iCell) = interp_sequence(x,y,1,slab_r8,1,field%nx, &
                                              1,field%ny,1,1,-1.e30_RKIND,interp_list,1)
          endif          
       end do

       deallocate(slab_r8)
       deallocate(field % slab)
!      exit
    end if
    call read_next_met_field(field,istatus)

 enddo
 write(0,*) '--- end subroutine physics_initialize_sst:'

 end subroutine physics_initialize_sst

!=============================================================================================
 subroutine physics_initialize_real(mesh,fg)
!=============================================================================================
!input arguments:
 type(mesh_type),intent(in) :: mesh

!inout arguments:
 type(fg_type),intent(inout):: fg 

!local variables:
 character(len=24):: initial_date

 integer:: iCell,nCellsSolve
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:)  ,pointer:: sfc_albbck
 real(kind=RKIND),dimension(:,:),pointer:: albedo12m

 real(kind=RKIND),dimension(:),pointer:: seaice,xice,xland
 real(kind=RKIND),dimension(:),pointer:: snoalb 
 real(kind=RKIND),dimension(:),pointer:: vegfra,shdmin,shdmax
 real(kind=RKIND),dimension(:),pointer:: snow,snowc,snowh
 real(kind=RKIND),dimension(:,:),pointer:: greenfrac

 real(kind=RKIND),dimension(:),pointer:: skintemp,sst
 
!---------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter physics_initialize_real:'

 nCellsSolve = mesh % nCellsSolve

 landmask   => mesh % landmask   % array
 albedo12m  => mesh % albedo12m  % array
 greenfrac  => mesh % greenfrac  % array
 shdmin     => mesh % shdmin     % array
 shdmax     => mesh % shdmax     % array
 snoalb     => mesh % snoalb     % array

 sfc_albbck => fg % sfc_albbck % array
 vegfra     => fg % vegfra     % array
 snow       => fg % snow       % array
 snowc      => fg % snowc      % array
 snowh      => fg % snowh      % array
 skintemp   => fg % skintemp   % array
 sst        => fg % sst        % array
 seaice     => fg % seaice     % array
 xice       => fg % xice       % array
 xland      => fg % xland      % array

!initialization of the sea-surface temperature and seaice if they are read from a separate
!input file. calling this subroutine will overwrite the arrays sst and seaice already read
!in the file defined by config_input_name:
 if(config_input_sst) then
    call physics_initialize_sst(mesh,fg)

    do iCell = 1, nCellsSolve
       !recalculate the sea-ice flag:
       if(xice(iCell) .gt. 0.) then
          seaice(iCell) = 1
       else
          seaice(iCell) = 0
       endif

       !set the skin temperature to the sea-surface temperature over the oceans:
       if(landmask(iCell).eq.0 .and. sst(iCell).gt.170. .and. sst(iCell).lt.400.) &
          skintemp(iCell) = sst(iCell)
    enddo
 endif

!initialization of the surface background albedo: interpolation of the monthly values to the
!initial date:
 initial_date = trim(config_start_time)
 call monthly_interp_to_date(nCellsSolve,initial_date,albedo12m,sfc_albbck)

 do iCell = 1, nCellsSolve
    sfc_albbck(iCell) = sfc_albbck(iCell) / 100.
    if(landmask(iCell) .eq. 0) sfc_albbck(iCell) = 0.08
 enddo

!initialization of the green-ness (vegetation) fraction: interpolation of the monthly values to
!the initial date. get the min/max for each cell for the monthly green-ness fraction:
!initial_date = trim(config_init_date)
 initial_date = trim(config_start_time)
 call monthly_interp_to_date(nCellsSolve,initial_date,greenfrac,vegfra)

!calculates the maximum and minimum green-ness (vegetation) fraction:
 call monthly_min_max(nCellsSolve,greenfrac,shdmin,shdmax)

!limit the annual maximum snow albedo to 0.08 over open-ocean and to 0.75 over sea-ice cells::
 do iCell = 1, nCellsSolve
    if(landmask(iCell) .eq. 0 .and. seaice(iCell) .eq. 0.) then
       snoalb(iCell) = 0.08
    elseif(landmask(iCell) .eq. 0 .and. seaice(iCell) .eq. 1.) then
       snoalb(iCell) = 0.75
    endif
 enddo

!initialization of the flag indicating the presence of snow (0 or 1) and of the snow depth
!(m) as functions of the input snow water content (kg/m2). we use a 5:1 ratio from liquid
!water equivalent to snow depth:
 do iCell = 1, nCellsSolve
    if(snow(iCell) .ge. 10.) then
       snowc(iCell) = 1.
    else
       snowc(iCell) = 0.
    endif
    snowh(iCell) = snow(iCell) * 5.0 / 1000.
 enddo

!initialization of soil layers properties:
 call init_soil_layers(mesh,fg)

!adjustment of all surface fields for seaice points:
 call init_seaice_points(mesh,fg)

!define xland over land and ocean:
 do iCell = 1, nCellsSolve
    if(landmask(iCell) .eq. 1 .or. (landmask(iCell).eq.0 .and. seaice(iCell).eq.1.)) then
       xland(iCell) = 1.
    elseif(landmask(iCell) .eq. 0) then
       xland(iCell) = 2.
    endif
 enddo

 write(0,*) '--- end physics_initialize_real:'

 end subroutine physics_initialize_real

!=============================================================================================
 subroutine init_soil_layers(mesh,fg)
!=============================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh

!inout arguments:
 type(fg_type),intent(inout):: fg

!---------------------------------------------------------------------------------------------

!adjust the annual mean deep soil temperature:
 call adjust_input_soiltemps(mesh,fg)

!initialize the depth of the soil layers:
 call init_soil_layers_depth(mesh,fg)
 
!initialize the temperature, moisture, and liquid water of the individual soil layers:
 call init_soil_layers_properties(mesh,fg)

 end subroutine init_soil_layers

!=============================================================================================
 subroutine adjust_input_soiltemps(mesh,fg)
!=============================================================================================

!input arguments:
 type(mesh_type),intent(in) :: mesh

!inout arguments:
 type(fg_type),intent(inout):: fg 

!local variables:
 integer:: iCell,ifgSoil
 integer:: nCellsSolve,nFGSoilLevels
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:),pointer  :: soilz,ter
 real(kind=RKIND),dimension(:),pointer  :: skintemp,soiltemp,tmn
 real(kind=RKIND),dimension(:,:),pointer:: st_fg

!---------------------------------------------------------------------------------------------

 nCellsSolve   = mesh % nCellsSolve
 nFGSoilLevels = mesh % nFGSoilLevels

 landmask => mesh % landmask % array
 soiltemp => mesh % soiltemp % array
 ter      => mesh % ter      % array

 skintemp => fg % skintemp % array
 tmn      => fg % tmn      % array
 st_fg    => fg % st_fg    % array
 soilz    => fg % soilz    % array

 do iCell = 1, nCellsSolve
    if(landmask(iCell) .eq. 1) then

       !adjust the annual deep mean soil temperature and skin temperatures over land: 
       tmn(iCell) = soiltemp(iCell) - 0.0065 * ter(iCell)
       skintemp(iCell) = skintemp(iCell) - 0.0065 * (ter(iCell)-soilz(iCell))

       !adjust the soil layer temperatures:
       do ifgSoil = 1, nFGSoilLevels
          st_fg(ifgSoil,iCell) = st_fg(ifgSoil,iCell) - 0.0065 * (ter(iCell)-soilz(iCell))
       enddo

    elseif(landmask(iCell) .eq. 0) then

       tmn(iCell) = skintemp(iCell)
       
    endif
 enddo

 end subroutine adjust_input_soiltemps

!=============================================================================================
 subroutine init_soil_layers_depth(mesh,fg)
!=============================================================================================

!input arguments:
 type(mesh_type),intent(in) :: mesh

!inout arguments:
 type(fg_type),intent(inout):: fg 

!local variables:
 integer:: iCell,iSoil

!---------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter subroutine init_soil_layers_depth:'

 if(config_nsoillevels .ne. 4) &
    call physics_error_fatal('NOAH lsm uses 4 soil layers. Correct config_nsoillevels.')

 do iCell = 1, mesh % nCells

    iSoil = 1
    fg % zs_fg % array(iSoil,iCell) = 0.5 * fg % dzs_fg % array(iSoil,iCell)
!   if(iCell .eq. 1) write(0,101) iSoil,fg % dzs_fg % array(iSoil,iCell), &
!                                 fg % zs_fg % array(iSoil,iCell)
    do iSoil = 2, mesh % nFGSoilLevels
       fg % zs_fg % array(iSoil,iCell) = fg % zs_fg % array(iSoil-1,iCell)        &
                                       + 0.5 * fg % dzs_fg % array(iSoil-1,iCell) &
                                       + 0.5 * fg % dzs_fg % array(iSoil,iCell)
!      if(iCell .eq. 1) write(0,101) iSoil,fg % dzs_fg % array(iSoil,iCell), &
!                                    fg % zs_fg % array(iSoil,iCell)
    enddo

 enddo
 101 format(i4,2(1x,e15.8))

 do iCell = 1, mesh % nCells
    fg % dzs % array(1,iCell) = 0.10
    fg % dzs % array(2,iCell) = 0.30
    fg % dzs % array(3,iCell) = 0.60
    fg % dzs % array(4,iCell) = 1.00

    iSoil = 1
    fg % zs % array(iSoil,iCell)  = 0.5 * fg % dzs % array(iSoil,iCell)
!   if(iCell .eq. 1) write(0,101) iSoil,fg % dzs % array(iSoil,iCell), &
!                                 fg % zs % array(iSoil,iCell)

    do iSoil = 2, mesh % nSoilLevels
       fg % zs % array(iSoil,iCell) = fg % zs % array(iSoil-1,iCell)        &
                                    + 0.5 * fg % dzs % array(iSoil-1,iCell) &
                                    + 0.5 * fg % dzs % array(iSoil,iCell)
!      if(iCell .eq. 1) write(0,101) iSoil,fg % dzs % array(iSoil,iCell),   &
!                                    fg % zs % array(iSoil,iCell)
    enddo

 enddo

 end subroutine init_soil_layers_depth

!=============================================================================================
 subroutine init_soil_layers_properties(mesh,fg)
!=============================================================================================

!input arguments:
 type(mesh_type),intent(in) :: mesh

!inout arguments:
 type(fg_type),intent(inout):: fg 

!local variables:
 integer:: iCell,ifgSoil,iSoil,is
 integer:: nCells,nFGSoilLevels,nSoilLevels
 integer,dimension(:),pointer:: landmask
 
 real(kind=RKIND),dimension(:,:),allocatable:: zhave,sm_input,st_input

 real(kind=RKIND),dimension(:),pointer  :: skintemp,tmn
 real(kind=RKIND),dimension(:,:),pointer:: dzs,zs,tslb,smois,sh2o
 real(kind=RKIND),dimension(:,:),pointer:: sm_fg,st_fg,zs_fg

!---------------------------------------------------------------------------------------------

!write(0,*)
 write(0,*) '--- enter subroutine init_soil_layers_properties:'

 nCells        = mesh % nCells
 nSoilLevels   = mesh % nSoilLevels
 nFGSoilLevels = mesh % nFGSoilLevels
 write(0,*) 'nSoilLevels   =',nSoilLevels
 write(0,*) 'nFGSoilLevels =',nFGSoilLevels

 landmask => mesh % landmask % array

 zs_fg => fg % zs_fg % array
 st_fg => fg % st_fg % array
 sm_fg => fg % sm_fg % array

 zs    => fg % zs  % array
 dzs   => fg % dzs % array 
 sh2o  => fg % sh2o  % array
 smois => fg % smois % array 
 tslb  => fg % tslb  % array
 skintemp => fg % skintemp % array
 tmn      => fg % tmn      % array

 if(config_nsoillevels .ne. 4) &
    call physics_error_fatal('NOAH lsm uses 4 soil layers. Correct config_nsoillevels.')

 if(.not.allocated(zhave)   ) allocate(zhave(nFGSoilLevels+2,nCells)   )
 if(.not.allocated(st_input)) allocate(st_input(nFGSoilLevels+2,nCells))
 if(.not.allocated(sm_input)) allocate(sm_input(nFGSoilLevels+2,nCells))

 do iCell = 1, nCells

    ifgSoil = 1
    zhave(ifgSoil,iCell)    = 0.
    st_input(ifgSoil,iCell) = skintemp(iCell)
    sm_input(ifgSoil,iCell) = sm_fg(ifgSoil+1,iCell)

    do ifgSoil = 1, nFGSoilLevels
       zhave(ifgSoil+1,iCell) = zs_fg(ifgSoil,iCell) / 100.
       st_input(ifgSoil+1,iCell) = st_fg(ifgSoil,iCell)
       sm_input(ifgSoil+1,iCell) = sm_fg(ifgSoil,iCell)
    enddo

    zhave(nFGSoilLevels+2,iCell) = 300./100.
    st_input(nFGSoilLevels+2,iCell) = tmn(iCell)
    sm_input(nFGSoilLevels+2,iCell) = sm_input(nFGSoilLevels,iCell)

    if(iCell .eq. 1) then
       do ifgSoil = 1,nFGSoilLevels+2
          write(0,101) ifgSoil,zhave(ifgSoil,iCell)
       enddo
    endif

 enddo

!... interpolate the soil temperature, soil moisture, and soil liquid temperature to the four
!    layers used in the NOAH land surface scheme:

 do iCell = 1, nCells

    if(landmask(iCell) .eq. 1) then

       noah: do iSoil = 1 , nSoilLevels
          input: do ifgSoil = 1 , nFGSoilLevels+2-1
             if(iCell .eq. 1) write(0,102) iSoil,ifgSoil,zs(iSoil,iCell), &
                zhave(ifgSoil,iCell),zhave(ifgSoil+1,iCell)

             if(zs(iSoil,iCell).ge.zhave(ifgSoil,iCell) .and. &
                zs(iSoil,iCell).le.zhave(ifgSoil+1,iCell)) then

                tslb(iSoil,iCell) = &
                      (st_input(ifgSoil,iCell) * (zhave(ifgSoil+1,iCell)-zs(iSoil,iCell))    &
                    +  st_input(ifgSoil+1,iCell) * (zs(iSoil,iCell)-zhave(ifgSoil,iCell)))   &
                         / (zhave(ifgSoil+1,iCell)-zhave(ifgSoil,iCell))
                if(iCell .eq. 1) write(6,102) iSoil,ifgSoil,zs(iSoil,iCell), &
                   zhave(ifgSoil,iCell),zhave(ifgSoil+1,iCell)
                         
                smois(iSoil,iCell) = &
                       (sm_input(ifgSoil,iCell) * (zhave(ifgSoil+1,iCell)-zs(iSoil,iCell))   &
                    +  sm_input(ifgSoil+1,iCell) * (zs(iSoil,iCell)-zhave(ifgSoil,iCell)))   &
                    / (zhave(ifgSoil+1,iCell)-zhave(ifgSoil,iCell))

                sh2o(iSoil,iCell)  = 0.

                exit input
             endif
          enddo input
          if(iCell.eq. 1) write(0,*)
       enddo noah

    elseif(landmask(iCell) .eq. 0) then

       !fill the soil temperatures with the skin temperatures over oceans:
       do iSoil = 1, nSoilLevels
          tslb(iSoil,iCell)  = skintemp(iCell)
          smois(iSoil,iCell) = 1.0
          sh2o(iSoil,iCell)  = 1.0
       enddo

    endif

 enddo

!... final checks:

 do iCell = 1, nCells

    if(landmask(iCell).eq. 1 .and. tslb(1,iCell).gt.170. .and. tslb(1,iCell).lt.400. .and. &
       smois(1,iCell).lt.0.005) then
       do iSoil = 1, nSoilLevels
          smois(iSoil,iCell) = 0.005
       enddo
    endif

 enddo

!formats:
 101 format(i4,4(1x,e15.8))
 102 format(2i5,5(1x,e15.8))

 if(allocated(zhave)   ) deallocate(zhave )
 if(allocated(st_input)) deallocate(st_input)
 if(allocated(sm_input)) deallocate(sm_input)

 end subroutine init_soil_layers_properties

!=============================================================================================
 subroutine init_seaice_points(mesh,fg)
!=============================================================================================

!input arguments:
 type(mesh_type),intent(in) :: mesh

!inout arguments:
 type(fg_type),intent(inout):: fg 

!local variables:
 character(len=128):: mess
 integer:: iCell,iSoil,nCellsSolve,nSoilLevels
 integer:: num_seaice_changes
 integer,dimension(:),pointer:: landmask,isltyp,ivgtyp

 real(kind=RKIND):: xice_threshold
 real(kind=RKIND):: mid_point_depth
 real(kind=RKIND),dimension(:),pointer  :: vegfra
 real(kind=RKIND),dimension(:),pointer  :: seaice,xice
 real(kind=RKIND),dimension(:),pointer  :: skintemp,sst,tmn
 real(kind=RKIND),dimension(:,:),pointer:: tslb,smois,sh2o

!note that this threshold is also defined in module_physics_vars.F.It is defined here to avoid
!adding "use module_physics_vars" since this subroutine is only used for the initialization of
!a "real" forecast with $CORE = init_nhyd_atmos.
 real(kind=RKIND),parameter:: xice_tsk_threshold = 271.
 real(kind=RKIND),parameter:: total_depth        = 3.   ! 3-meter soil depth.

!---------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter init_seaice_points:'
 write(0,*) '--- config_frac_seaice      :', config_frac_seaice

 nCellsSolve = mesh % nCellsSolve
 nSoilLevels = mesh % nSoilLevels

 landmask => mesh % landmask    % array
 isltyp   => mesh % soilcat_top % array
 ivgtyp   => mesh % lu_index    % array

 seaice   => fg   % seaice      % array
 xice     => fg   % xice        % array
 vegfra   => fg   % vegfra      % array

 skintemp => fg   % skintemp    % array
 sst      => fg   % sst         % array
 tmn      => fg   % tmn         % array

 tslb     => fg   % tslb        % array
 smois    => fg   % smois       % array
 sh2o     => fg   % sh2o        % array

 if(.not. config_frac_seaice) then
    xice_threshold = 0.5
 elseif(config_frac_seaice) then
    xice_threshold = 0.02
 endif

!make sure that all the cells flagged as sea-ice cells are defined as ocean cells:
 num_seaice_changes = 0
 do iCell = 1, nCellsSolve
    if((landmask(iCell).eq.1 .and. xice(iCell).gt.0.) .or. xice(iCell).gt.200.) then
       num_seaice_changes = num_seaice_changes + 1
       seaice(iCell) = 0.
       xice(iCell)   = 0.
    endif
 enddo
 write(mess,fmt='(A,i12)') 'number of seaice cells converted to land cells=', &
       num_seaice_changes
 call physics_message(mess)

 num_seaice_changes = 0
 do iCell =1 , nCellsSolve

    if(xice(iCell) .ge. xice_threshold .or. &
       (landmask(iCell).eq.0 .and. skintemp(iCell).lt.xice_tsk_threshold)) then

       num_seaice_changes = num_seaice_changes + 1
       !sea-ice points are converted to land points:
       if(.not. config_frac_seaice) xice(iCell) = 1.0
       if(landmask(iCell) .eq. 0) tmn(iCell) = 271.4

       ivgtyp(iCell)   = 24 ! (isice = 24)
       isltyp(iCell)   = 16
       vegfra(iCell)   = 0.
       landmask(iCell) = 1.

       do iSoil = 1, nSoilLevels
          mid_point_depth = total_depth/nSoilLevels/2. &
                          + (iSoil-1)*(total_depth/nSoilLevels)
          tslb(iSoil,iCell) = ((total_depth-mid_point_depth) * skintemp(iCell) &
                            +  mid_point_depth * tmn(iCell)) / total_depth
          smois(iSoil,iCell) = 1.0
          sh2o(iSoil,iCell)  = 0.0
       enddo
       
    elseif(xice(iCell) .lt. xice_threshold) then
       xice(iCell) = 0.

    endif

 enddo

 end subroutine init_seaice_points

!=============================================================================================
 end module mpas_atmphys_initialize_real
!=============================================================================================

