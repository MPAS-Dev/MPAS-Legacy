!=============================================================================================
 module mpas_atmphys_driver_convection_deep
 use mpas_grid_types

 use mpas_atmphys_constants
 use mpas_atmphys_vars

!wrf physics:
 use module_cu_kfeta

 implicit none
 private
 public:: allocate_convection_deep,       &
          deallocate_convection_deep,     &
          init_convection_deep,           &
          driver_convection_deep,         &
          update_convection_deep

 integer, private:: i,k,j

 contains

!=============================================================================================
 subroutine allocate_convection_deep
!=============================================================================================

 if(.not.allocated(cu_act_flag)) allocate(cu_act_flag(ims:ime,jms:jme)       )
 if(.not.allocated(rthcuten_p) ) allocate(rthcuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqvcuten_p) ) allocate(rqvcuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqccuten_p) ) allocate(rqccuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqrcuten_p) ) allocate(rqrcuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqicuten_p) ) allocate(rqicuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqscuten_p) ) allocate(rqscuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(cubot_p)    ) allocate(cubot_p(ims:ime,jms:jme)           )
 if(.not.allocated(cutop_p)    ) allocate(cutop_p(ims:ime,jms:jme)           )
 if(.not.allocated(pratec_p)   ) allocate(pratec_p(ims:ime,jms:jme)          )
 if(.not.allocated(raincv_p)   ) allocate(raincv_p(ims:ime,jms:jme)          )

 convection_select: select case(conv_deep_scheme)

    case ("kain_fritsch")
       if(.not.allocated(area_p) ) allocate(area_p(ims:ime,jms:jme)         )
       if(.not.allocated(nca_p)  ) allocate(nca_p(ims:ime,jms:jme)          )
       if(.not.allocated(w0avg_p)) allocate(w0avg_p(ims:ime,kms:kme,jms:jme))

    case default

 end select convection_select

 end subroutine allocate_convection_deep

!=============================================================================================
 subroutine deallocate_convection_deep
!=============================================================================================

 if(allocated(cu_act_flag)) deallocate(cu_act_flag)
 if(allocated(rthcuten_p) ) deallocate(rthcuten_p )
 if(allocated(rqvcuten_p) ) deallocate(rqvcuten_p )
 if(allocated(rqccuten_p) ) deallocate(rqccuten_p )
 if(allocated(rqrcuten_p) ) deallocate(rqrcuten_p )
 if(allocated(rqicuten_p) ) deallocate(rqicuten_p )
 if(allocated(rqscuten_p) ) deallocate(rqscuten_p )
 if(allocated(cubot_p)    ) deallocate(cubot_p    )
 if(allocated(cutop_p)    ) deallocate(cutop_p    )
 if(allocated(pratec_p)   ) deallocate(pratec_p   )
 if(allocated(raincv_p)   ) deallocate(raincv_p   )

 convection_select: select case(conv_deep_scheme)

    case ("kain_fritsch")
       if(allocated(area_p) ) deallocate(area_p )
       if(allocated(nca_p)  ) deallocate(nca_p  )
       if(allocated(w0avg_p)) deallocate(w0avg_p)

    case default

 end select convection_select

 end subroutine deallocate_convection_deep

!=============================================================================================
 subroutine init_convection_deep(config_do_restart,mesh,diag_physics)
!=============================================================================================

!input arguments:
 logical,intent(in):: config_do_restart
 type(mesh_type),intent(in):: mesh

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics

!local variables:
 integer:: iCell

!---------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter convection_deep initialization:'

 convection_select: select case(conv_deep_scheme)

    case ("kain_fritsch")
       write(0,*) '    enter kain-fritsch initialization:'
       if(.not. config_do_restart) then
          do iCell = 1, mesh % nCells
             diag_physics % nca % array(iCell) = -100.
          enddo
       endif
       call kf_lutab(svp1,svp2,svp3,svpt0)
       write(0,*) '    end kain-kritsch initialization'

    case default

 end select convection_select

 write(0,*) '--- end convection_deep initialization'

 end subroutine init_convection_deep

!=============================================================================================
 subroutine driver_convection_deep(itimestep,mesh,diag_physics,tend_physics)
!=============================================================================================

!input and output arguments:
!---------------------------
 integer,intent(in):: itimestep
 type(mesh_type),intent(in):: mesh
 type(diag_physics_type),intent(inout):: diag_physics
 type(tend_physics_type),intent(inout):: tend_physics

!local variables and arrays:
!---------------------------
 logical:: log_convection
 integer:: icount
 real(kind=RKIND):: dx

!variables specific to Kain_Fritsch parameterization:
 logical:: warm_rain,adapt_step_flag
 real(kind=RKIND):: curr_secs
 real(kind=RKIND):: cudt
 
!=============================================================================================
 write(0,*)
 write(0,*) '--- enter convection_driver: dt_cu=',dt_cu

!initialize instantaneous precipitation, and copy convective tendencies from the dynamics to
!the physics grid:

 call convection_from_MPAS(dt_dyn,mesh,diag_physics,tend_physics)

!... convert the convection time-step to minutes:
 cudt = dt_cu/60.

!... call to convection schemes:
!dx = sqrt(maxval(mesh % areaCell % array))

 do j = jts, jte
 do i = its, ite
    cu_act_flag(i,j) = .false.
 enddo
 enddo

 convection_select: select case(conv_deep_scheme)

    case ("kain_fritsch")

       !initialization:
       curr_secs = -1
       adapt_step_flag = .false.
       write(0,*) '--- enter subroutine kf_eta_cps:'
       call  kf_eta_cps ( &
             dt        = dt_dyn     , ktau            = itimestep       ,            &
!            dx        = dx         , cudt            = dt_cu           ,            &
             areaCell  = area_p     , cudt            = dt_cu           ,            &
             curr_secs = curr_secs  , adapt_step_flag = adapt_step_flag ,            &
             rho       = rho_p      , raincv          = raincv_p        ,            &
             pratec    = pratec_p   , nca             = nca_p           ,            &
             u         = u_p        , v               = v_p             ,            &
             th        = th_p       , t               = t_p             ,            &
             w         = w_p        , dz8w            = dz_p            ,            &
             pcps      = pres_p     , pi              = pi_p            ,            &
             w0avg     = w0avg_p    , xlv0            = xlv0            ,            &
             xlv1      = xlv1       , xls0            = xls0            ,            &
             xls1      = xls1       , cp              = cp              ,            &
             r         = r_d        , g               = g               ,            &
             ep1       = ep_1       , ep2             = ep_2            ,            &
             svp1      = svp1       , svp2            = svp2            ,            &
             svp3      = svp3       , svpt0           = svpt0           ,            &
             stepcu    = n_cu       , cu_act_flag     = cu_act_flag     ,            &
             warm_rain = warm_rain  , cutop           = cutop_p         ,            &
             cubot     = cubot_p    , qv              = qv_p            ,            &
             f_qv      = f_qv       , f_qc            = f_qc            ,            &
             f_qr      = f_qr       , f_qi            = f_qi            ,            &
             f_qs      = f_qs       , rthcuten        = rthcuten_p      ,            &
             rqvcuten  = rqvcuten_p , rqccuten        = rqccuten_p      ,            &
             rqrcuten  = rqrcuten_p , rqicuten        = rqicuten_p      ,            &
             rqscuten  = rqscuten_p ,                                                &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme , &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte   &
                    )
       201 format(i3,i6,1x,l1,5(1x,e15.8))
       write(0,*) '--- end subroutine kf_eta_cps:'
!      write(0,*)
!      write(0,*) '--- deep convection:'
!      do j = jts,jte
!      do i = its,ite
!         if(nca_p(i,j).gt.0. .and. raincv_p(i,j).gt.0.) then
!            write(0,201) j,i,cu_act_flag(i,j),nca_p(i,j),raincv_p(i,j), &
!                         raincv_p(i,j)/dt_dyn,pratec_p(i,j)
!         endif
!      enddo
!      enddo
!      write(0,*) '--- shallow convection:'
!      do j = jts,jte
!      do i = its,ite
!         if(nca_p(i,j).gt.0. .and. raincv_p(i,j).eq.0.) then
!            write(0,201) j,i,cu_act_flag(i,j),nca_p(i,j),raincv_p(i,j), &
!                         raincv_p(i,j)/dt_dyn,pratec_p(i,j)
!         endif
!      enddo
!      enddo

    case default

 end select convection_select

!copy instantaneous and accumulated precipitation, convective tendencies, and "other" arrays
!specific to convection parameterization back to the dynamics grid:

 call convection_to_MPAS(diag_physics,tend_physics)

 write(0,*) '--- end subroutine convection_driver'

 end subroutine driver_convection_deep

!=============================================================================================
 subroutine convection_from_MPAS(dt_dyn,mesh,diag_physics,tend_physics)
!=============================================================================================
!input arguments:
 type(mesh_type),intent(in):: mesh
 type(diag_physics_type),intent(in):: diag_physics
 type(tend_physics_type),intent(in):: tend_physics
 real(kind=RKIND),intent(in):: dt_dyn

!---------------------------------------------------------------------------------------------

 do j = jts,jte
 do i = its,ite
    cubot_p(i,j)  = diag_physics % cubot % array(i)
    cutop_p(i,j)  = diag_physics % cutop % array(i)
    raincv_p(i,j) = diag_physics % raincv % array(i)
    pratec_p(i,j) = diag_physics % cuprec % array(i)
    do k = kts, kte
       rthcuten_p(i,k,j) = tend_physics % rthcuten % array(k,i)
       rqvcuten_p(i,k,j) = tend_physics % rqvcuten % array(k,i)
       rqccuten_p(i,k,j) = tend_physics % rqccuten % array(k,i)
       rqrcuten_p(i,k,j) = tend_physics % rqrcuten % array(k,i)
       rqicuten_p(i,k,j) = tend_physics % rqicuten % array(k,i)
       rqscuten_p(i,k,j) = tend_physics % rqscuten % array(k,i)
    enddo
 enddo
 enddo
 
 convection_select: select case(conv_deep_scheme)

    case ("kain_fritsch")
 
       do j = jts,jte
       do i = its,ite
          !area of grid-cell:       
          area_p(i,j) = mesh % areaCell % array(i)

          !decreases the characteristic time period that convection remains active. When nca_p
          !becomes less than the convective timestep, convective tendencies and precipitation
          !are reset to zero (note that this is also done in subroutine kf_eta_cps).
          nca_p(i,j) = diag_physics % nca % array(i)

          if(nca_p(i,j) .gt. 0.) then
             nca_p(i,j) = nca_p(i,j) - dt_dyn
             
             if(nca_p(i,j) .lt. 0.5*dt_dyn) then
                do k = kts,kte
                   rthcuten_p(i,k,j) = 0.
                   rqvcuten_p(i,k,j) = 0.
                   rqccuten_p(i,k,j) = 0.
                   rqrcuten_p(i,k,j) = 0.
                   rqicuten_p(i,k,j) = 0.
                   rqscuten_p(i,k,j) = 0.
                enddo
                raincv_p(i,j) = 0.
                pratec_p(i,j) = 0.
                cubot_p(i,j)  = kte+1
                cutop_p(i,j)  = kts
             endif
          endif

          do k = kts, kte
             w0avg_p(i,k,j) = diag_physics % w0avg % array(k,i)
          enddo
       enddo
       enddo

    case default

 end select convection_select
 
 end subroutine convection_from_MPAS
 
!=============================================================================================
 subroutine convection_to_MPAS(diag_physics,tend_physics)
!=============================================================================================
!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics
 type(tend_physics_type),intent(inout):: tend_physics

!---------------------------------------------------------------------------------------------

 do j = jts,jte
 do i = its,ite
    diag_physics % raincv % array(i) = raincv_p(i,j)
    diag_physics % cuprec % array(i) = pratec_p(i,j)
    diag_physics % cubot % array(i)  = cubot_p(i,j)
    diag_physics % cutop % array(i)  = cutop_p(i,j)
    do k = kts, kte
       tend_physics % rthcuten % array(k,i) = rthcuten_p(i,k,j)
       tend_physics % rqvcuten % array(k,i) = rqvcuten_p(i,k,j)
       tend_physics % rqccuten % array(k,i) = rqccuten_p(i,k,j)
       tend_physics % rqrcuten % array(k,i) = rqrcuten_p(i,k,j)
       tend_physics % rqicuten % array(k,i) = rqicuten_p(i,k,j)
       tend_physics % rqscuten % array(k,i) = rqscuten_p(i,k,j)
    enddo
 enddo
 enddo

 convection_select: select case(conv_deep_scheme)

    case ("kain_fritsch")
       do j = jts,jte
       do i = its,ite
          diag_physics % nca   % array(i) = nca_p(i,j)
          do k = kts, kte
             diag_physics % w0avg % array(k,i) = w0avg_p(i,k,j)
          enddo                          
       enddo
       enddo

    case default

 end select convection_select

 end subroutine convection_to_MPAS

!=============================================================================================
 subroutine update_convection_deep(dt_dyn,mesh,diag_physics)
!=============================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh
 real(kind=RKIND),intent(in):: dt_dyn

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics

!local variables:
 integer:: iCell

!---------------------------------------------------------------------------------------------

!update the accumuluted precipitation rate at the end of each dynamic time-step:
 do iCell = 1, mesh % nCells
    diag_physics % rainc % array(iCell) = diag_physics % rainc % array(iCell) &
                                    + diag_physics % cuprec % array(iCell) * dt_dyn
 enddo

 end subroutine update_convection_deep

!=============================================================================================
 end module mpas_atmphys_driver_convection_deep
!=============================================================================================