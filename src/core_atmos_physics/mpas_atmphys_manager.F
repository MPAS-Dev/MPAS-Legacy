!=============================================================================================
 module mpas_atmphys_manager
 use mpas_configure
 use mpas_grid_types
 use mpas_timekeeping

 use mpas_atmphys_constants
 use mpas_atmphys_vars
 use mpas_atmphys_update_surface
 use mpas_atmphys_utilities

 implicit none
 private
 public:: physics_timetracker,physics_run_init

 integer, public:: year                 !Current year.
 integer, public:: julday               !Initial Julian day.
 real(kind=RKIND), public:: curr_julday !Current Julian day (= 0.0 at 0Z on January 1st).
 real(kind=RKIND), public:: gmt         !Greenwich mean time hour of model start (hr)

 integer, private:: i,k,j,n

!defines alarms for calling the long- and short-wave radiation codes, for calling the convec-
!tion scheme:
 integer, parameter:: radtlwAlarmID = 11
 integer, parameter:: radtswAlarmID = 12
 integer, parameter:: convAlarmID   = 13
 integer, parameter:: pblAlarmID    = 14

!defines alarm to update the surface boundary conditions:
 integer, parameter:: sfcbdyAlarmID = 15

!defines alarm to update the background surface albedo and the greeness fraction:
 integer, parameter:: greenAlarmID  = 16

!defines alarm to update the ozone path length,the trace gas path length,the total emissivity,
!and the total absorptivity in the "CAM" long-wave radiation codes. The default time interval
!between updates is 6 hours and is set with config_camrad_abs_update (00:30:00).
 integer, parameter:: camAlarmID    = 17

 integer :: h, m, s, s_n, s_d, DoY, yr
 real(kind=RKIND) :: utc_h

 contains

!=============================================================================================
 subroutine physics_timetracker(domain,dt,clock,itimestep,xtime_s)
!=============================================================================================

!input arguments:
 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: dt
 real(kind=RKIND),intent(in) :: xtime_s

!inout arguments:
 type(MPAS_Clock_type),intent(inout):: clock
 type(domain_type),intent(inout)    :: domain

!local variables:
 type(block_type),pointer:: block
 type(MPAS_Time_Type):: currTime

 logical:: LeapYear
 character(len=32)   :: timeStamp
 integer:: ierr
 real(kind=RKIND):: utc_s
 real(kind=RKIND):: xtime_m

!=============================================================================================
 write(0,*)
 write(0,*) '--- enter subroutine physics_timetracker: itimestep = ', itimestep

!update the current julian day and current year:
 100 format('     YEAR        =', i5   ,/, &
            '     JULDAY      =', i5   ,/, &
            '     GMT         =', f16.9,/, &
            '     UTC_H       =', f16.9,/, &
            '     CURR_JULDAY =', f16.9,/, &
            '     LEAP_YEAR   =', 1x,l1,/)

 currTime = mpas_get_clock_time(clock,MPAS_NOW,ierr)
 call mpas_get_time(curr_time=currTime,dateTimeString=timeStamp,YYYY=yr,H=h,M=m, &
                    S=s,S_n=s_n,S_d=s_d,DoY=DoY,ierr=ierr)

 utc_h = real(h) + real(m) / 60.0 + real(s + s_n / s_d) / 3600.0
 utc_s = real(s + s_n + s_d) 
 year   = yr
 julday = DoY
 curr_julday = real(julday-1) + utc_h / 24.0
 LeapYear = isLeapYear(year)
 write(0,100) year,julday,gmt,utc_h,curr_julday,LeapYear

 block => domain % blocklist
 do while(associated(block))

    !update the background surface albedo and greeness of vegetation: interpolation of input
    !monthly values to current day:
    if(mpas_is_alarm_ringing(clock,greenAlarmID,ierr=ierr)) then
       call mpas_reset_clock_alarm(clock,greenAlarmID,ierr=ierr)
       write(0,*) '--- update background surface albedo, greeness fraction:', timeStamp
       call physics_update_surface(timeStamp,block%mesh,block%sfc_input)
    endif

    !update surface boundary conditions with input sea-surface temperatures and fractional
    !sea-ice coverage:
    if(mpas_is_alarm_ringing(clock,sfcbdyAlarmID,ierr=ierr)) then
       call mpas_reset_clock_alarm(clock,sfcbdyAlarmID,ierr=ierr)
       if(config_sst_update) &
          call physics_update_sst(block%mesh,block%sfc_input,block%diag_physics)
    endif   

    !apply a diurnal cycle to the sea-surface temperature:
    if(config_sstdiurn_update) &
       call physics_update_sstskin(dt_dyn,block%mesh,block%diag_physics,block%sfc_input)

    !update the deep soil temperature:
    if(config_deepsoiltemp_update) &
       call physics_update_deepsoiltemp(LeapYear,dt_dyn,curr_julday,block%mesh, &
                    block%sfc_input,block%diag_physics)

    block => block % next
 end do 

!check to see if it is time to run the longwave and shortwave radiation codes:
 if(trim(config_radt_lw_scheme) /= "off") then
    l_radtlw = .false.

    if(config_radtlw_interval /= "none") then
       if(mpas_is_alarm_ringing(clock,radtlwAlarmID,ierr=ierr)) then
          call mpas_reset_clock_alarm(clock,radtlwAlarmID,ierr=ierr)
          l_radtlw = .true.
       endif
    elseif(config_radtlw_interval == "none") then
       l_radtlw = .true.
    endif
    write(0,*)
    write(0,*) '--- time to run the LW radiation scheme L_RADLW =',l_radtlw
 endif

 if(trim(config_radt_sw_scheme) /= "off") then
    l_radtsw = .false.

    if(config_radtsw_interval /= "none") then
       if(mpas_is_alarm_ringing(clock,radtswAlarmID,ierr=ierr)) then
          call mpas_reset_clock_alarm(clock,radtswAlarmID,ierr=ierr)
          l_radtsw = .true.
       endif
    elseif(config_radtsw_interval == "none") then
       l_radtsw = .true.
    endif
    write(0,*) '--- time to run the SW radiation scheme L_RADSW =',l_radtsw
 endif

!check to see if it is time to update the ozone trace gas path lengths,the total emissivity,
!and the total absorptivity in the "CAM" long-wave radiation codes.
 if(trim(config_radt_lw_scheme) .eq. "cam_lw" .or. &
    trim(config_radt_sw_scheme) .eq. "cam_sw" ) then

    doabsems = .false.
    if(mpas_is_alarm_ringing(clock,camAlarmID,ierr=ierr)) then
       call mpas_reset_clock_alarm(clock,camAlarmID,ierr=ierr)
       doabsems = .true.
    endif
    write(0,*) '--- update CAM absorptivity and emissivity arrays DOABSEMS =',doabsems

 endif
 
!formats:
 101 format(3x,'l_radtlw = ',l1,3x,'l_radtsw = ',l1)

 end subroutine physics_timetracker

!=============================================================================================
 subroutine physics_run_init(mesh,state,clock)
!=============================================================================================

!input arguments:
 type(mesh_type),intent(in) :: mesh
 type(state_type),intent(in):: state
 type(MPAS_Clock_type):: clock

!local variables:
 type(MPAS_Time_Type):: startTime,alarmStartTime
 type(MPAS_TimeInterval_Type):: alarmTimeStep
 integer:: ierr

!=============================================================================================
 write(0,*)
 write(0,*) '--- enter subroutine physics_run_init:'

!initialization of gmt, julian day, and alarms:
 101 format('     YEAR        =', i5   ,/, &
            '     JULDAY      =', i5   ,/, &
            '     GMT         =', f16.9,/, &
            '     UTC_H       =', f16.9,/, &
            '     CURR_JULDAY =', f16.9,/)

 102 format('     DT_RADTLW   =', f16.9,/, &
            '     DT_RADTSW   =', f16.9,/, &
            '     DT_CU       =', f16.9,/, &
            '     DT_PBL      =', f16.9,/)
            
 103 format('     IMS =', i4,3x,'IME =', i8,/,  &
            '     JMS =', i4,3x,'JME =', i8,/,  &
            '     KMS =', i4,3x,'KME =', i8,//, &
            '     IDS =', i4,3x,'IDE =', i8,/,  &
            '     JDS =', i4,3x,'JDE =', i8,/,  &
            '     KDS =', i4,3x,'KDE =', i8,//, &
            '     ITS =', i4,3x,'ITE =', i8,/,  &
            '     JTS =', i4,3x,'JTE =', i8,/,  &
            '     KTS =', i4,3x,'KTE =', i8)

 startTime = mpas_get_clock_time(clock, MPAS_START_TIME, ierr)
 call mpas_get_time(curr_time=startTime,YYYY=yr,H=h,M=m,S=s,S_n=s_n,S_d=s_d,DoY=DoY,ierr=ierr)
 utc_h  = real(h) + real(m) / 60.0 + real(s + s_n / s_d) / 3600.0
 year   = yr
 gmt    = utc_h
 julday = DoY
 curr_julday = real(julday-1)
 write(0,101) year,julday,gmt,utc_h,curr_julday

!set alarms for calling the longwave and shortwave radiation schemes, the convection schemes,
!and the PBL schemes at intervals different (greater) than the dynamical time-step:

 if(trim(config_radtlw_interval) /= "none") then
    call mpas_set_timeInterval(alarmTimeStep,timeString=config_radtlw_interval,ierr=ierr)
    alarmStartTime = startTime
    call mpas_add_clock_alarm(clock,radtlwAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_run_init: error creating radtlwAlarmID')

    call mpas_get_timeInterval(interval=alarmTimeStep,dt=dt_radtlw,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_run_init: error defining dt_radtlw')

 elseif(trim(config_radtlw_interval) == "none") then
    dt_radtlw = config_dt

 else
    call physics_error_fatal('subroutine physics_run_init: dt_radtlw is not defined')

 endif

 if(trim(config_radtsw_interval) /= "none") then
    call mpas_set_timeInterval(alarmTimeStep,timeString=config_radtsw_interval,ierr=ierr)
    alarmStartTime = startTime
    call mpas_add_clock_alarm(clock,radtswAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_run_init: error creating alarm radtsw')

    call mpas_get_timeInterval(interval=alarmTimeStep,dt=dt_radtsw,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_run_init: error defining radtswAlarmID')

 elseif(trim(config_radtsw_interval) == "none") then
    dt_radtsw = config_dt

 else
    call physics_error_fatal('subroutine physics_run_init: dt_radtsw is not defined')

 endif

 if(trim(config_conv_interval) /= "none") then
    call mpas_set_timeInterval(alarmTimeStep,timeString=config_conv_interval,ierr=ierr)
    alarmStartTime = startTime
    call mpas_add_clock_alarm(clock,convAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_init: error creating convAlarmID')

    call mpas_get_timeInterval(interval=alarmTimeStep,dt=dt_cu,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_run_init: error defining dt_cu')

 elseif(trim(config_conv_interval) == "none") then
    dt_cu = config_dt

 else
    call physics_error_fatal('subroutine physics_run_init: dt_cu is not defined')

 endif

 if(trim(config_pbl_interval) /= "none") then
    call mpas_set_timeInterval(alarmTimeStep,timeString=config_pbl_interval,ierr=ierr)
    alarmStartTime = startTime
    call mpas_add_clock_alarm(clock,pblAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_init: error creating pblAlarmID')

    call mpas_get_timeInterval(interval=alarmTimeStep,dt=dt_pbl,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_run_init: error defining dt_pbl')

 elseif(trim(config_conv_interval) == "none") then
    dt_pbl = config_dt

 else
    call physics_error_fatal('subroutine physics_run_init: dt_pbl is not defined')

 endif

!set alarm for updating the background surface albedo and the greeness fraction:
 call mpas_set_timeInterval(alarmTimeStep,timeString=config_greeness_update,ierr=ierr)
 alarmStartTime = startTime
 call mpas_add_clock_alarm(clock,greenAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_init: error creating alarm greeness')

!set alarm for updating the surface boundary conditions:
 if(trim(config_sfc_update_interval) /= "none") then
    call mpas_set_timeInterval(alarmTimeStep,timeString=config_sfc_update_interval,ierr=ierr)
    alarmStartTime = startTime
    call mpas_add_clock_alarm(clock,sfcbdyAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_init: error creating alarm sfcbdy')
 endif

!set alarm to update the ozone path length, the trace gas path length, the total emissivity,
!and the total absorptivity in the "CAM" long-wave radiation codes.
 if(trim(config_radt_lw_scheme) .eq. "cam_lw" .or. &
    trim(config_radt_sw_scheme) .eq. "cam_sw" ) then
    call mpas_set_timeInterval(alarmTimeStep,timeString=config_camrad_abs_update,ierr=ierr)
    alarmStartTime = startTime
    call mpas_add_clock_alarm(clock,camAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_init: error creating alarm CAM')
 endif

 write(0,102) dt_radtlw,dt_radtsw,dt_cu,dt_pbl

!initialization of physics dimensions to mimic a rectangular grid:
 ims=1   ; ime = mesh % nCellsSolve
 jms=1   ; jme=1
 kms=1   ; kme = mesh % nVertLevels+1

 ids=ims ; ide=ime
 jds=jms ; jde=jme
 kds=kms ; kde=kme 

 its=ims ; ite = ime 
 jts=jms ; jte = jme
 kts=kms ; kte = kme-1

 write(0,103) ims,ime,jms,jme,kms,kme, &
              ids,ide,jds,jde,kds,kde, &
              its,ite,jts,jte,kts,kte

!initialization:
 num_months = mesh % nMonths

!initialization of physics time-steps:
 dt_dyn     = config_dt
 n_microp   = config_n_microp
 n_cu       = config_n_conv

 dt_microp  = dt_dyn/n_microp  !for now.

!write(0,*) 'mod =',mod(dt_dyn,dt_radtsw)
!write(0,*) 'mod =',mod(dt_dyn,dt_microp)
!stop

!cloud microphysics scheme:
 microp_scheme       = trim(config_microp_scheme)
 conv_deep_scheme    = trim(config_conv_deep_scheme)
 conv_shallow_scheme = trim(config_conv_shallow_scheme)
 sfclayer_scheme     = trim(config_sfclayer_scheme)
 pbl_scheme          = trim(config_pbl_scheme)
 radt_cld_scheme     = trim(config_radt_cld_scheme)
 radt_lw_scheme      = trim(config_radt_lw_scheme)
 radt_sw_scheme      = trim(config_radt_sw_scheme)

!CAM radiation schemes:
 if(trim(config_radt_lw_scheme) .eq. "cam_lw" .or. &
    trim(config_radt_sw_scheme) .eq. "cam_sw" ) then

    !this logical initiates the calculation of "pre-dfined arrays in the longwave
    !radiation code. these arrays are calculated once if it is the beginning of a
    !new run or if it is a restart run:
    doabsems     = .false.
    cam_abs_dim1 = mesh % cam_dim1
    cam_abs_dim2 = mesh % nVertLevels + 1

    !initializes the number of aerosols, and the prescribed vertical dimensions for
    !aerosols and ozone mixing ratios:
    num_aerosols  = state % num_aerosols
    num_aerLevels = mesh % naerLevels
    num_oznLevels = mesh % noznLevels

 endif 

!land-surface scheme:
 lsm_scheme = trim(config_lsm_scheme)
 num_soils  = mesh% nSoilLevels

 if(.not. config_frac_seaice) then
    xice_threshold = 0.5
 elseif(config_frac_seaice) then
    xice_threshold = 0.02
 endif

 end subroutine physics_run_init

!=============================================================================================
 end module mpas_atmphys_manager
!=============================================================================================
 
