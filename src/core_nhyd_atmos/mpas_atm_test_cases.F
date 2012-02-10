module atm_test_cases

   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_dmpar
   use atm_advection
#ifdef DO_PHYSICS
   use mpas_atmphys_control
#endif


   contains


   subroutine atm_setup_test_case(domain)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Configure grid metadata and model state for the hydrostatic test case
   !   specified in the namelist
   !
   ! Output: block - a subset (not necessarily proper) of the model domain to be
   !                 initialized
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain

      integer :: i
      type (block_type), pointer :: block_ptr

      if (config_test_case == 0) then
         write(0,*) ' Using initial conditions from input file'
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            do i=2,nTimeLevs
               call mpas_copy_state(block_ptr % state % time_levs(i) % state, block_ptr % state % time_levs(1) % state)
            end do
            block_ptr => block_ptr % next
         end do

      else if ((config_test_case == 1) .or. (config_test_case == 2) .or. (config_test_case == 3)) then
         write(0,*) ' Jablonowski and Williamson baroclinic wave test case '
         if (config_test_case == 1) write(0,*) ' no initial perturbation '
         if (config_test_case == 2) write(0,*) ' initial perturbation included '
         if (config_test_case == 3) write(0,*) ' normal-mode perturbation included '
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            write(0,*) ' calling test case setup '
            call atm_test_case_jw(block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % diag, &
                                   block_ptr % diag_physics ,config_test_case)
            write(0,*) ' returned from test case setup '
            do i=2,nTimeLevs
               call mpas_copy_state(block_ptr % state % time_levs(i) % state, block_ptr % state % time_levs(1) % state)
            end do

            block_ptr => block_ptr % next
         end do

      else if ((config_test_case == 4) .or. (config_test_case ==5)) then

         write(0,*) ' squall line - super cell test case '
         if (config_test_case == 4) write(0,*) ' squall line test case' 
         if (config_test_case == 5) write(0,*) ' supercell test case'
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            write(0,*) ' calling test case setup '
            call atm_test_case_squall_line(domain % dminfo, block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % diag, config_test_case)
            write(0,*) ' returned from test case setup '
            do i=2,nTimeLevs
               call mpas_copy_state(block_ptr % state % time_levs(i) % state, block_ptr % state % time_levs(1) % state)
            end do

            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 6 ) then

         write(0,*) ' mountain wave test case '
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            write(0,*) ' calling test case setup '
            call atm_test_case_mtn_wave(block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % diag, config_test_case)
            write(0,*) ' returned from test case setup '
            do i=2,nTimeLevs
               call mpas_copy_state(block_ptr % state % time_levs(i) % state, block_ptr % state % time_levs(1) % state)
            end do

            block_ptr => block_ptr % next
         end do

      else


         write(0,*) ' Only test case 1, 2, 3, 4, 5 and 6 are currently supported for nonhydrostatic core '
         stop
      end if


#ifdef DO_PHYSICS
      !initialization of surface input variables technically not needed to run our current set of
      !idealized test cases:
      if (config_test_case > 0)  then

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call physics_idealized_init(block_ptr % mesh, block_ptr % sfc_input)
            block_ptr => block_ptr % next
         end do

      endif
#endif

   end subroutine atm_setup_test_case

!----------------------------------------------------------------------------------------------------------

   subroutine atm_test_case_jw(grid, state, diag, diag_physics, test_case)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup baroclinic wave test case from Jablonowski and Williamson 2008 (QJRMS)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mesh_type), intent(inout) :: grid
      type (state_type), intent(inout) :: state
      type (diag_type), intent(inout) :: diag
      type (diag_physics_type), intent(inout) :: diag_physics
      integer, intent(in) :: test_case

      real (kind=RKIND), parameter :: u0 = 35.0
      real (kind=RKIND), parameter :: alpha_grid = 0.  ! no grid rotation
      real (kind=RKIND), parameter :: omega_e = 7.29212e-05
      real (kind=RKIND), parameter :: t0b = 250., t0 = 288., delta_t = 4.8e+05, dtdz = 0.005, eta_t = 0.2
      real (kind=RKIND), parameter :: u_perturbation = 1., pert_radius = 0.1, latitude_pert = 40., longitude_pert = 20.
      real (kind=RKIND), parameter :: theta_c = pii/4.0
      real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
      real (kind=RKIND), parameter :: k_x = 9.           ! Normal mode wave number

      real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp
      real (kind=RKIND), dimension(:), pointer :: surface_pressure
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zx, zz, hx
      real (kind=RKIND), dimension(:,:), pointer :: pressure, ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, dss, t, rt
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      
!.. initialization of moisture:
      integer:: index_qv
      real (kind=RKIND),parameter :: rh_max = 0.40 ! Maximum relative humidity
!      real (kind=RKIND),parameter :: rh_max = 0.70 ! Maximum relative humidity
      real (kind=RKIND),dimension(:,:),  pointer:: qsat, relhum
      real (kind=RKIND),dimension(:,:,:),pointer:: scalars
!.. end initialization of moisture.

      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, nz1, itr, itrp, cell1, cell2, nCellsSolve

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnEdge 
      integer, dimension(:,:), pointer :: edgesOnEdge, CellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: dvEdge, AreaCell 
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      real (kind=RKIND) :: u, v, flux, fluxk, lat1, lat2, eta_v, r_pert, u_pert, lat_pert, lon_pert, r

      real (kind=RKIND) :: ptop, p0, phi
      real (kind=RKIND) :: lon_Edge

      real (kind=RKIND) :: r_earth, etavs, ztemp, zd, zt, dz, gam, delt, str

      real (kind=RKIND) :: es, qvs, xnutr, znut, ptemp 
      integer :: iter

      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: hyai, hybi, znu, znw, znwc, znwv, hyam, hybm
      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: znuc, znuv, bn, divh, dpn

      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: sh, zw, ah
      real (kind=RKIND), dimension(grid % nVertLevels ) :: zu, dzw, rdzwp, rdzwm
      real (kind=RKIND), dimension(grid % nVertLevels ) :: eta, etav, teta, ppi, tt, temperature_1d

      real (kind=RKIND) :: d1, d2, d3, cf1, cf2, cf3, cof1, cof2, psurf

      !  storage for (lat,z) arrays for zonal velocity calculation

      logical, parameter :: rebalance = .true.
      integer, parameter :: nlat=721
      real (kind=RKIND), dimension(grid % nVertLevels) :: flux_zonal
      real (kind=RKIND), dimension(grid % nVertLevels + 1, nlat) :: zgrid_2d
      real (kind=RKIND), dimension(grid % nVertLevels, nlat) :: u_2d, pp_2d, rho_2d, qv_2d, etavs_2d, zz_2d
      real (kind=RKIND), dimension(grid % nVertLevels, nlat-1) :: zx_2d 
      real (kind=RKIND), dimension(nlat) :: lat_2d
      real (kind=RKIND) :: dlat, hx_1d
      real (kind=RKIND) :: z_edge, z_edge3, d2fdx2_cell1, d2fdx2_cell2

      logical, parameter :: moisture = .true.
      !
      ! Scale all distances and areas from a unit sphere to one with radius a
      !
      grid % xCell % array = grid % xCell % array * a
      grid % yCell % array = grid % yCell % array * a
      grid % zCell % array = grid % zCell % array * a
      grid % xVertex % array = grid % xVertex % array * a
      grid % yVertex % array = grid % yVertex % array * a
      grid % zVertex % array = grid % zVertex % array * a
      grid % xEdge % array = grid % xEdge % array * a
      grid % yEdge % array = grid % yEdge % array * a
      grid % zEdge % array = grid % zEdge % array * a
      grid % dvEdge % array = grid % dvEdge % array * a
      grid % dcEdge % array = grid % dcEdge % array * a
      grid % areaCell % array = grid % areaCell % array * a**2.0
      grid % areaTriangle % array = grid % areaTriangle % array * a**2.0
      grid % kiteAreasOnVertex % array = grid % kiteAreasOnVertex % array * a**2.0

      weightsOnEdge     => grid % weightsOnEdge % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array
      dvEdge            => grid % dvEdge % array
      AreaCell          => grid % AreaCell % array
      CellsOnEdge       => grid % CellsOnEdge % array

      deriv_two  => grid % deriv_two % array
      zb  => grid % zb % array
      zb3 => grid % zb3% array
      
      nz1 = grid % nVertLevels
      nz = nz1 + 1
      nCellsSolve = grid % nCellsSolve

      zgrid => grid % zgrid % array
      rdzw => grid % rdzw % array
      dzu => grid % dzu % array
      rdzu => grid % rdzu % array
      fzm => grid % fzm % array
      fzp => grid % fzp % array
      zx => grid % zx % array
      zz => grid % zz % array
      hx => grid % hx % array
      dss => grid % dss % array

      pb => diag % exner_base % array
      rb => diag % rho_base % array
      tb => diag % theta_base % array
      rtb => diag % rtheta_base % array
      p => diag % exner % array

      ppb => diag % pressure_base % array
      pp  => diag % pressure_p % array

      rho_zz => state % rho_zz % array
      rr => diag % rho_p % array
      t => state % theta_m % array      
      rt => diag % rtheta_p % array

      surface_pressure => diag % surface_pressure % array

!.. initialization of moisture:
      scalars => state % scalars % array
      qsat    => diag_physics % qsat % array
      relhum  => diag_physics % relhum % array
      scalars(:,:,:) = 0.0
      qsat(:,:)      = 0.0
      relhum(:,:)    = 0.0
      qv_2d(:,:)     = 0.0
!.. end initialization of moisture.

      surface_pressure(:) = 0.0

      call atm_initialize_advection_rk(grid) 
      call atm_initialize_deformation_weights(grid) 

      index_qv = state % index_qv

      xnutr = 0.
      zd = 12000.
      znut = eta_t

      etavs = (1.-0.252)*pii/2.
      r_earth = a
      p0 = 1.e+05

      write(0,*) ' point 1 in test case setup '

! We may pass in an hx(:,:) that has been precomputed elsewhere.
! For now it is independent of k

      do iCell=1,grid % nCells
        do k=1,nz
          phi = grid % latCell % array (iCell)
          hx(k,iCell) = u0/gravity*cos(etavs)**1.5                                   &
                      *((-2.*sin(phi)**6                                   &
                            *(cos(phi)**2+1./3.)+10./63.)                  &
                            *(u0)*cos(etavs)**1.5                          &
                       +(1.6*cos(phi)**3                                   &
                            *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)
        enddo
      enddo

      !     Metrics for hybrid coordinate and vertical stretching

      str = 1.8
      zt = 45000.
      dz = zt/float(nz1)

      write(0,*) ' hx computation complete '

      do k=1,nz
		
!           sh(k) is the stretching specified for height surfaces

            sh(k) = (real(k-1)*dz/zt)**str 
				
!           to specify specific heights zc(k) for coordinate surfaces,
!           input zc(k) and define sh(k) = zc(k)/zt
!           zw(k) is the hieght of zeta surfaces
!                zw(k) = (k-1)*dz yields constant dzeta
!                        and nonconstant dzeta/dz
!                zw(k) = sh(k)*zt yields nonconstant dzeta
!                        and nearly constant dzeta/dz 

            zw(k) = float(k-1)*dz
!            zw(k) = sh(k)*zt
!
!           ah(k) governs the transition between terrain-following 
!           and pureheight coordinates
!                ah(k) = 0 is a terrain-following coordinate
!                ah(k) = 1 is a height coordinate
 
            ah(k) = 1.-cos(.5*pii*(k-1)*dz/zt)**6
!            ah(k) = 0.
	    write(0,*) ' k, sh, zw, ah ',k,sh(k),zw(k),ah(k)			
      end do
      do k=1,nz1
         dzw (k) = zw(k+1)-zw(k)
         rdzw(k) = 1./dzw(k)
         zu(k  ) = .5*(zw(k)+zw(k+1))
      end do
      do k=2,nz1
         dzu (k)  = .5*(dzw(k)+dzw(k-1))
         rdzu(k)  =  1./dzu(k)
         fzp (k)  = .5* dzw(k  )/dzu(k)
         fzm (k)  = .5* dzw(k-1)/dzu(k)
         rdzwp(k) = dzw(k-1)/(dzw(k  )*(dzw(k)+dzw(k-1)))
         rdzwm(k) = dzw(k  )/(dzw(k-1)*(dzw(k)+dzw(k-1)))
      end do

!**********  how are we storing cf1, cf2 and cf3?

      COF1 = (2.*DZU(2)+DZU(3))/(DZU(2)+DZU(3))*DZW(1)/DZU(2) 
      COF2 =     DZU(2)        /(DZU(2)+DZU(3))*DZW(1)/DZU(3) 
      CF1  = FZP(2) + COF1
      CF2  = FZM(2) - COF1 - COF2
      CF3  = COF2       

!      d1  = .5*dzw(1)
!      d2  = dzw(1)+.5*dzw(2)
!      d3  = dzw(1)+dzw(2)+.5*dzw(3)
!      cf1 = d2*d3*(d3-d2)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf2 = d1*d3*(d1-d3)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf3 = d1*d2*(d2-d1)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))

      write(0,*) ' cf1, cf2, cf3 = ',cf1,cf2,cf3

      grid % cf1 % scalar = cf1
      grid % cf2 % scalar = cf2
      grid % cf3 % scalar = cf3

      do iCell=1,grid % nCells
        do k=1,nz	
          zgrid(k,iCell) = (1.-ah(k))*(sh(k)*(zt-hx(k,iCell))+hx(k,iCell))  &
                         + ah(k) * sh(k)* zt	
        end do
        do k=1,nz1
          zz (k,iCell) = (zw(k+1)-zw(k))/(zgrid(k+1,iCell)-zgrid(k,iCell))
        end do
      end do

      do i=1, grid % nEdges
        iCell1 = grid % CellsOnEdge % array(1,i)
        iCell2 = grid % CellsOnEdge % array(2,i)
        do k=1,nz
          zx (k,i) = (zgrid(k,iCell2)-zgrid(k,iCell1)) / grid % dcEdge % array(i)
        end do
      end do
      do i=1, grid % nCells
        do k=1,nz1
          ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
          dss(k,i) = 0.
          ztemp = zgrid(k,i)
          if(ztemp.gt.zd+.1)  then
             dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
          end if
        end do
      enddo

      !do k=1,nz1
      !  write(0,*) ' k, zgrid(k,1),hx(k,1) ',k,zgrid(k,1),hx(k,1)
      !enddo

      !do k=1,nz1
      !  write(0,*) ' k, zx(k,1) ',k,zx(k,1)
      !enddo

      write(0,*) ' grid metrics setup complete '

!**************  section for 2d (z,lat) calc for zonal velocity

      dlat = 0.5*pii/float(nlat-1)
      do i = 1,nlat

        lat_2d(i) = float(i-1)*dlat
        phi = lat_2d(i)
        hx_1d    = u0/gravity*cos(etavs)**1.5                            &
                   *((-2.*sin(phi)**6                                   &
                         *(cos(phi)**2+1./3.)+10./63.)                  &
                         *(u0)*cos(etavs)**1.5                          &
                    +(1.6*cos(phi)**3                                   &
                         *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)

        do k=1,nz	
          zgrid_2d(k,i) = (1.-ah(k))*(sh(k)*(zt-hx_1d)+hx_1d)  &
                         + ah(k) * sh(k)* zt	
        end do
        do k=1,nz1
          zz_2d (k,i) = (zw(k+1)-zw(k))/(zgrid_2d(k+1,i)-zgrid_2d(k,i))
        end do

        do k=1,nz1
          ztemp    = .5*(zgrid_2d(k+1,i)+zgrid_2d(k,i))
          ppb(k,i) = p0*exp(-gravity*ztemp/(rgas*t0b)) 
          pb (k,i) = (ppb(k,i)/p0)**(rgas/cp)
          rb (k,i) = ppb(k,i)/(rgas*t0b*zz_2d(k,i))
          tb (k,i) = t0b/pb(k,i)
          rtb(k,i) = rb(k,i)*tb(k,i)
          p  (k,i) = pb(k,i)
          pp (k,i) = 0.
          rr (k,i) = 0.
        end do


        do itr = 1,10

          do k=1,nz1
            eta (k) = (ppb(k,i)+pp(k,i))/p0
            etav(k) = (eta(k)-.252)*pii/2.
            if(eta(k).ge.znut)  then
              teta(k) = t0*eta(k)**(rgas*dtdz/gravity)
            else
              teta(k) = t0*eta(k)**(rgas*dtdz/gravity) + delta_t*(znut-eta(k))**5
            end if
          end do

          phi = lat_2d (i)
          do k=1,nz1
            temperature_1d(k) = teta(k)+.75*eta(k)*pii*u0/rgas*sin(etav(k))      &
                            *sqrt(cos(etav(k)))*                   &
                              ((-2.*sin(phi)**6                    &
                                   *(cos(phi)**2+1./3.)+10./63.)   &
                                   *2.*u0*cos(etav(k))**1.5        &
                              +(1.6*cos(phi)**3                    &
                                *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)/(1.+0.61*qv_2d(k,i))


            ztemp   = .5*(zgrid_2d(k,i)+zgrid_2d(k+1,i))
            ptemp   = ppb(k,i) + pp(k,i)

            !get moisture 
            if (moisture) then
              qv_2d(k,i) = env_qv( ztemp, temperature_1d(k), ptemp, rh_max )
            end if

            tt(k) = temperature_1d(k)*(1.+1.61*qv_2d(k,i))
          end do

          do itrp = 1,25
            do k=1,nz1				
              rr(k,i)  = (pp(k,i)/(rgas*zz_2d(k,i)) - rb(k,i)*(tt(k)-t0b))/tt(k)
            end do

            ppi(1) = p0-.5*dzw(1)*gravity                            &
                          *(1.25*(rr(1,i)+rb(1,i))*(1.+qv_2d(1,i))   &
                            -.25*(rr(2,i)+rb(2,i))*(1.+qv_2d(2,i)))

            ppi(1) = ppi(1)-ppb(1,i)
            do k=1,nz1-1
              ppi(k+1) = ppi(k)-.5*dzu(k+1)*gravity*                        &
                            (rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*qv_2d(k  ,i)   &
                            +rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*qv_2d(k+1,i))
            end do

            do k=1,nz1
              pp(k,i) = .2*ppi(k)+.8*pp(k,i)
            end do

          end do  ! end inner iteration loop itrp

        end do  ! end outer iteration loop itr

        do k=1,nz1
          rho_2d(k,i) = rr(k,i)+rb(k,i)
          pp_2d(k,i) = pp(k,i)
          etavs_2d(k,i) = ((ppb(k,i)+pp(k,i))/p0 - 0.252)*pii/2.
          u_2d(k,i) = u0*(sin(2.*lat_2d(i))**2) *(cos(etavs_2d(k,i))**1.5)
        end do

      end do  ! end loop over latitudes for 2D zonal wind field calc

      !SHP-balance:: in case of rebalacing for geostrophic wind component
      if (rebalance) then

        do i=1,nlat-1
          do k=1,nz1
            zx_2d(k,i) = (zgrid_2d(k,i+1)-zgrid_2d(k,i))/(dlat*r_earth)
          end do
        end do

        call atm_recompute_geostrophic_wind(u_2d,rho_2d,pp_2d,qv_2d,lat_2d,zz_2d,zx_2d,     &
                                        cf1,cf2,cf3,fzm,fzp,rdzw,nz1,nlat,dlat)

      end if

!******************************************************************      

!
!---- baroclinc wave initialization ---------------------------------
!
!     reference sounding based on dry isothermal atmosphere
!
      do i=1, grid % nCells
        do k=1,nz1
          ztemp    = .5*(zgrid(k+1,i)+zgrid(k,i))
          ppb(k,i) = p0*exp(-gravity*ztemp/(rgas*t0b)) 
          pb (k,i) = (ppb(k,i)/p0)**(rgas/cp)
          rb (k,i) = ppb(k,i)/(rgas*t0b*zz(k,i))
          tb (k,i) = t0b/pb(k,i)
          rtb(k,i) = rb(k,i)*tb(k,i)
          p  (k,i) = pb(k,i)
          pp (k,i) = 0.
          rr (k,i) = 0.
        end do

!       if(i == 1) then
!         do k=1,nz1
!           write(0,*) ' k, ppb, pb, rb, tb (k,1) ',k,ppb(k,1),pb(k,1),rb(k,1)*zz(k,1),tb(k,1)
!         enddo
!       end if

      200 format(4i6,8(1x,e15.8))
      201 format(3i6,8(1x,e15.8))
      202 format(2i6,10(1x,e15.8))
      203 format(i6,10(1x,e15.8))

!     iterations to converge temperature as a function of pressure
!
        do itr = 1,10

          do k=1,nz1
            eta (k) = (ppb(k,i)+pp(k,i))/p0
            etav(k) = (eta(k)-.252)*pii/2.
            if(eta(k).ge.znut)  then
              teta(k) = t0*eta(k)**(rgas*dtdz/gravity)
            else
              teta(k) = t0*eta(k)**(rgas*dtdz/gravity) + delta_t*(znut-eta(k))**5
            end if
          end do
          phi = grid % latCell % array (i)
          do k=1,nz1
            temperature_1d(k) = teta(k)+.75*eta(k)*pii*u0/rgas*sin(etav(k))      &
                            *sqrt(cos(etav(k)))*                   &
                              ((-2.*sin(phi)**6                    &
                                   *(cos(phi)**2+1./3.)+10./63.)   &
                                   *2.*u0*cos(etav(k))**1.5        &
                              +(1.6*cos(phi)**3                    &
                                *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)/(1.+0.61*scalars(index_qv,k,i))

            ztemp   = .5*(zgrid(k,i)+zgrid(k+1,i))
            ptemp   = ppb(k,i) + pp(k,i)

            !get moisture 
            if (moisture) then
 
                !scalars(index_qv,k,i) = env_qv( ztemp, temperature_1d(k), ptemp, rh_max )

               if(ptemp < 50000.) then
                  relhum(k,i) = 0.0
               elseif(ptemp > p0) then
                  relhum(k,i) = 1.0
               else
                  relhum(k,i) = (1.-((p0-ptemp)/50000.)**1.25)
               endif
               relhum(k,i) = min(rh_max,relhum(k,i))

               !.. calculation of water vapor mixing ratio:
               if (temperature_1d(k) > 273.15) then
                   es  = 1000.*0.6112*exp(17.67*(temperature_1d(k)-273.15)/(temperature_1d(k)-29.65))
               else
                   es  = 1000.*0.6112*exp(21.8745584*(temperature_1d(k)-273.15)/(temperature_1d(k)-7.66))
               end if
               qsat(k,i) = (287.04/461.6)*es/(ptemp-es)
               if(relhum(k,i) .eq. 0.0) qsat(k,i) = 0.0
               scalars(index_qv,k,i) = relhum(k,i)*qsat(k,i)
            end if

            tt(k) = temperature_1d(k)*(1.+1.61*scalars(index_qv,k,i))

          end do
		
          do itrp = 1,25
            do k=1,nz1				
              rr(k,i)  = (pp(k,i)/(rgas*zz(k,i)) - rb(k,i)*(tt(k)-t0b))/tt(k)
            end do

            ppi(1) = p0-.5*dzw(1)*gravity                         &
                          *(1.25*(rr(1,i)+rb(1,i))*(1.+scalars(index_qv,1,i))   &
                            -.25*(rr(2,i)+rb(2,i))*(1.+scalars(index_qv,2,i)))

            ppi(1) = ppi(1)-ppb(1,i)
            do k=1,nz1-1
              ppi(k+1) = ppi(k)-.5*dzu(k+1)*gravity*                     &
                            (rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*scalars(index_qv,k  ,i)   &
                            +rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*scalars(index_qv,k+1,i))
            end do

            do k=1,nz1
              pp(k,i) = .2*ppi(k)+.8*pp(k,i)
            end do

          end do  ! end inner iteration loop itrp

        end do  ! end outer iteration loop itr

        do k=1,nz1	
          p (k,i) = ((ppb(k,i)+pp(k,i))/p0)**(rgas/cp)
          t (k,i) = tt(k)/p(k,i)
          rt (k,i) = t(k,i)*rr(k,i)+rb(k,i)*(t(k,i)-tb(k,i))
          rho_zz (k,i) = rb(k,i) + rr(k,i)
        end do

        !calculation of surface pressure:
        surface_pressure(i) = 0.5*dzw(1)*gravity                                    &
                        * (1.25*(rr(1,i) + rb(1,i)) * (1. + scalars(index_qv,1,i))  &
                        -  0.25*(rr(2,i) + rb(2,i)) * (1. + scalars(index_qv,2,i)))
        surface_pressure(i) = surface_pressure(i) + pp(1,i) + ppb(1,i)

      end do  ! end loop over cells

      !write(0,*)
      !write(0,*) '--- initialization of water vapor:'
      !do iCell = 1, grid % nCells
      !   if(iCell == 1 .or. iCell == grid % nCells) then
      !      do k = nz1, 1, -1
      !         write(0,202) iCell,k,t(k,iCell),relhum(k,iCell),qsat(k,iCell),scalars(index_qv,k,iCell)
      !      enddo
      !      write(0,*)
      !   endif
      !enddo

      lat_pert = latitude_pert*pii/180.
      lon_pert = longitude_pert*pii/180.

      do iEdge=1,grid % nEdges

         vtx1 = grid % VerticesOnEdge % array (1,iEdge)
         vtx2 = grid % VerticesOnEdge % array (2,iEdge)
         lat1 = grid%latVertex%array(vtx1)
         lat2 = grid%latVertex%array(vtx2)
         iCell1 = grid % cellsOnEdge % array(1,iEdge)
         iCell2 = grid % cellsOnEdge % array(2,iEdge)
         flux = (0.5*(lat2-lat1) - 0.125*(sin(4.*lat2) - sin(4.*lat1)))*a/grid % dvEdge % array(iEdge)

         if (config_test_case == 2) then
            r_pert = sphere_distance( grid % latEdge % array (iEdge), grid % lonEdge % array (iEdge), &
                                      lat_pert, lon_pert, 1.0_RKIND)/(pert_radius)
            u_pert = u_perturbation*exp(-r_pert**2)*(lat2-lat1)*a/grid % dvEdge % array(iEdge)

         else if (config_test_case == 3) then
            lon_Edge = grid % lonEdge % array(iEdge)
            u_pert = u_perturbation*cos(k_x*(lon_Edge - lon_pert)) &
                         *(0.5*(lat2-lat1) - 0.125*(sin(4.*lat2) - sin(4.*lat1)))*a/grid % dvEdge % array(iEdge)
         else
            u_pert = 0.0
         end if

         if (rebalance) then

           call atm_calc_flux_zonal(u_2d,etavs_2d,lat_2d,flux_zonal,lat1,lat2,grid % dvEdge % array(iEdge),a,u0,nz1,nlat)
           do k=1,grid % nVertLevels
             fluxk = u0*flux_zonal(k)/(0.5*(rb(k,iCell1)+rb(k,iCell2)+rr(k,iCell1)+rr(k,iCell2)))
             state % u % array(k,iEdge) = fluxk + u_pert
           end do

         else 

           do k=1,grid % nVertLevels
             etavs = (0.5*(ppb(k,iCell1)+ppb(k,iCell2)+pp(k,iCell1)+pp(k,iCell2))/p0 - 0.252)*pii/2.
             fluxk = u0*flux*(cos(etavs)**1.5)
             state % u % array(k,iEdge) = fluxk + u_pert
           end do

         end if

         cell1 = grid % CellsOnEdge % array(1,iEdge)
         cell2 = grid % CellsOnEdge % array(2,iEdge)
         do k=1,nz1
            diag % ru % array (k,iEdge)  = 0.5*(rho_zz(k,cell1)+rho_zz(k,cell2))*state % u % array (k,iEdge)
         end do

      !
      ! Generate rotated Coriolis field
      !

         grid % fEdge % array(iEdge) = 2.0 * omega * &
                                       ( -cos(grid%lonEdge%array(iEdge)) * cos(grid%latEdge%array(iEdge)) * sin(alpha_grid) + &
                                         sin(grid%latEdge%array(iEdge)) * cos(alpha_grid) &
                                       )
      end do

      do iVtx=1,grid % nVertices
         grid % fVertex % array(iVtx) = 2.0 * omega * &
                                         (-cos(grid%lonVertex%array(iVtx)) * cos(grid%latVertex%array(iVtx)) * sin(alpha_grid) + &
                                          sin(grid%latVertex%array(iVtx)) * cos(alpha_grid) &
                                         )
      end do

      !
      !  CALCULATION OF OMEGA, RW = ZX * RU + ZZ * RW
      !

      !
      !     pre-calculation z-metric terms in omega eqn.
      !
      do iEdge = 1,grid % nEdges
         cell1 = CellsOnEdge(1,iEdge)
         cell2 = CellsOnEdge(2,iEdge)

         do k = 1, grid%nVertLevels

            if (config_theta_adv_order == 2) then

               z_edge = (zgrid(k,cell1)+zgrid(k,cell2))/2.

            else if (config_theta_adv_order == 3 .or. config_theta_adv_order ==4) then !theta_adv_order == 3 or 4 

               d2fdx2_cell1 = deriv_two(1,1,iEdge) * zgrid(k,cell1)
               d2fdx2_cell2 = deriv_two(1,2,iEdge) * zgrid(k,cell2)
               do i=1, grid % nEdgesOnCell % array (cell1)
                  d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * zgrid(k,grid % CellsOnCell % array (i,cell1))
                  d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * zgrid(k,grid % CellsOnCell % array (i,cell2))
               end do

               z_edge =  0.5*(zgrid(k,cell1) + zgrid(k,cell2))         &
                             - (grid % dcEdge % array(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.

               if (config_theta_adv_order == 3) then
                  z_edge3 =  - (grid % dcEdge % array(iEdge) **2) * (d2fdx2_cell1 - d2fdx2_cell2) / 12.
               else
                  z_edge3 = 0.
               end if

            end if

               zb(k,1,iEdge) = (z_edge-zgrid(k,cell1))*dvEdge(iEdge)/AreaCell(cell1)
               zb(k,2,iEdge) = (z_edge-zgrid(k,cell2))*dvEdge(iEdge)/AreaCell(cell2)
               zb3(k,1,iEdge)=  z_edge3*dvEdge(iEdge)/AreaCell(cell1)
               zb3(k,2,iEdge)=  z_edge3*dvEdge(iEdge)/AreaCell(cell2)

         end do

      end do

      ! for including terrain
      diag % rw % array = 0.
      state % w % array = 0.
      do iEdge = 1,grid % nEdges

         cell1 = CellsOnEdge(1,iEdge)
         cell2 = CellsOnEdge(2,iEdge)

         do k = 2, grid%nVertLevels
            flux =  (fzm(k)*diag % ru % array(k,iEdge)+fzp(k)*diag % ru % array(k-1,iEdge))
            diag % rw % array(k,cell2) = diag % rw % array(k,cell2) + (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb(k,2,iEdge)*flux
            diag % rw % array(k,cell1) = diag % rw % array(k,cell1) - (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb(k,1,iEdge)*flux

            if (config_theta_adv_order ==3) then 
               diag % rw % array(k,cell2) = diag % rw % array(k,cell2)    &
                                            - sign(1.0_RKIND,diag % ru % array(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb3(k,2,iEdge)*flux
               diag % rw % array(k,cell1) = diag % rw % array(k,cell1)    &
                                            + sign(1.0_RKIND,diag % ru % array(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb3(k,1,iEdge)*flux
            end if

         end do

      end do

      ! Compute w from rho_zz and rw
      do iCell=1,grid%nCells
         do k=2,grid%nVertLevels
            state % w % array(k,iCell) = diag % rw % array(k,iCell)     &
                                       / (fzp(k) * state % rho_zz % array(k-1,iCell) + fzm(k) * state % rho_zz % array(k,iCell))
         end do
      end do


      !
      ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
      !
      diag % v % array(:,:) = 0.0
      do iEdge = 1, grid%nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            do k = 1, grid%nVertLevels
               diag % v % array(k,iEdge) = diag % v %array(k,iEdge) + weightsOnEdge(i,iEdge) * state % u % array(k, eoe)
           end do
         end do
      end do

      do i=1,10
        psurf = (cf1*(ppb(1,i)+pp(1,i)) + cf2*(ppb(2,i)+pp(2,i)) + cf3*(ppb(3,i)+pp(3,i)))/100.

            psurf = (ppb(1,i)+pp(1,i)) + .5*dzw(1)*gravity        &
                          *(1.25*(rr(1,i)+rb(1,i))*(1.+scalars(index_qv,1,i))   &
                            -.25*(rr(2,i)+rb(2,i))*(1.+scalars(index_qv,2,i)))

        write(0,*) ' i, psurf, lat ',i,psurf,grid%latCell%array(i)*180./3.1415828
      enddo

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            diag % rho % array(k,iCell) = state % rho_zz % array(k,iCell) * zz(k,iCell)
            diag % theta % array(k,iCell) = state % theta_m % array(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do

   end subroutine atm_test_case_jw

   subroutine atm_calc_flux_zonal(u_2d,etavs_2d,lat_2d,flux_zonal,lat1_in,lat2_in,dvEdge,a,u0,nz1,nlat)

   implicit none
   integer, intent(in) :: nz1,nlat
   real (kind=RKIND), dimension(nz1,nlat), intent(in) :: u_2d,etavs_2d
   real (kind=RKIND), dimension(nlat), intent(in) :: lat_2d
   real (kind=RKIND), dimension(nz1), intent(out) :: flux_zonal
   real (kind=RKIND), intent(in) :: lat1_in, lat2_in, dvEdge, a, u0

   integer :: k,i
   real (kind=RKIND) :: lat1, lat2, w1, w2
   real (kind=RKIND) :: dlat,da,db

   lat1 = abs(lat1_in)
   lat2 = abs(lat2_in)
   if(lat2 <= lat1) then
     lat1 = abs(lat2_in)
     lat2 = abs(lat1_in)
   end if

   do k=1,nz1
     flux_zonal(k) = 0.
   end do

   do i=1,nlat-1
     if( (lat1 <= lat_2d(i+1)) .and. (lat2 >= lat_2d(i)) ) then

     dlat = lat_2d(i+1)-lat_2d(i)
     da = (max(lat1,lat_2d(i))-lat_2d(i))/dlat
     db = (min(lat2,lat_2d(i+1))-lat_2d(i))/dlat
     w1 = (db-da) -0.5*(db-da)**2
     w2 = 0.5*(db-da)**2

     do k=1,nz1
       flux_zonal(k) = flux_zonal(k) + w1*u_2d(k,i) + w2*u_2d(k,i+1)
     end do

     end if

   end do

!  renormalize for setting cell-face fluxes

   do k=1,nz1
     flux_zonal(k) = sign(1.0_RKIND,lat2_in-lat1_in)*flux_zonal(k)*dlat*a/dvEdge/u0
   end do
     
   end subroutine atm_calc_flux_zonal

   !SHP-balance
   subroutine atm_recompute_geostrophic_wind(u_2d,rho_2d,pp_2d,qv_2d,lat_2d,zz_2d,zx_2d,     &
                                         cf1,cf2,cf3,fzm,fzp,rdzw,nz1,nlat,dlat)

   implicit none
   integer, intent(in) :: nz1,nlat
   real (kind=RKIND), dimension(nz1,nlat), intent(inout) :: u_2d
   real (kind=RKIND), dimension(nz1,nlat), intent(in) :: rho_2d, pp_2d, qv_2d, zz_2d
   real (kind=RKIND), dimension(nz1,nlat-1), intent(in) :: zx_2d
   real (kind=RKIND), dimension(nlat), intent(in) :: lat_2d
   real (kind=RKIND), dimension(nz1), intent(in) :: fzm, fzp, rdzw
   real (kind=RKIND), intent(in) :: cf1, cf2, cf3, dlat

   !local variable
   real (kind=RKIND), dimension(nz1,nlat-1) :: pgrad, ru, u
   real (kind=RKIND), dimension(nlat-1) :: f
   real (kind=RKIND), dimension(nz1+1)  :: dpzx

   real (kind=RKIND), parameter :: omega_e = 7.29212e-05
   real (kind=RKIND) :: rdx, qtot, r_earth, phi
   integer :: k,i, itr

   r_earth  = a
   rdx = 1./(dlat*r_earth)

   do i=1,nlat-1
     do k=1,nz1
       pgrad(k,i) = rdx*(pp_2d(k,i+1)/zz_2d(k,i+1)-pp_2d(k,i)/zz_2d(k,i))
     end do

     dpzx(:) = 0.

     k=1
     dpzx(k) = .5*zx_2d(k,i)*(cf1*(pp_2d(k  ,i+1)+pp_2d(k  ,i))        &
                             +cf2*(pp_2d(k+1,i+1)+pp_2d(k+1,i))        &
                             +cf3*(pp_2d(k+2,i+1)+pp_2d(k+2,i)))
     do k=2,nz1
        dpzx(k) = .5*zx_2d(k,i)*(fzm(k)*(pp_2d(k  ,i+1)+pp_2d(k  ,i))   &
                                +fzp(k)*(pp_2d(k-1,i+1)+pp_2d(k-1,i)))
     end do

     do k=1,nz1
        pgrad(k,i) = pgrad(k,i) - rdzw(k)*(dpzx(k+1)-dpzx(k))
     end do
   end do


   !initial value of v and rv -> that is from analytic sln. 
   do i=1,nlat-1
      do k=1,nz1
         u(k,i) = .5*(u_2d(k,i)+u_2d(k,i+1))
         ru(k,i) = u(k,i)*(rho_2d(k,i)+rho_2d(k,i+1))*.5
      end do
   end do

   write(0,*) "MAX U wind before REBALANCING ---->", maxval(abs(u))

   !re-calculate geostrophic wind using iteration 
   do itr=1,50
   do i=1,nlat-1
      phi = (lat_2d(i)+lat_2d(i+1))/2.
      f(i) = 2.*omega_e*sin(phi)
      do k=1,nz1
         if (f(i).eq.0.) then
           ru(k,i) = 0.
         else
           qtot = .5*(qv_2d(k,i)+qv_2d(k,i+1))
           ru(k,i) = - ( 1./(1.+qtot)*pgrad(k,i) + tan(phi)/r_earth*u(k,i)*ru(k,i) )/f(i)
         end if
           u(k,i) = ru(k,i)*2./(rho_2d(k,i)+rho_2d(k,i+1))
      end do
   end do
   end do

   write(0,*) "MAX U wind after REBALANCING ---->", maxval(abs(u))

   !update 2d ru
   do i=2,nlat-1
     do k=1,nz1
       u_2d(k,i) = (ru(k,i-1)+ru(k,i))*.5
     end do
   end do

   i=1
   do k=1,nz1
      u_2d(k,i) = (3.*u_2d(k,i+1)-u_2d(k,i+2))*.5
   end do
   i=nlat
   do k=1,nz1
      u_2d(k,i) = (3.*u_2d(k,i-1)-u_2d(k,i-2))*.5
   end do


   end subroutine atm_recompute_geostrophic_wind
!----------------------------------------------------------------------------------------------------------

   subroutine atm_test_case_squall_line(dminfo, grid, state, diag, test_case)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup squall line and supercell test case
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (dm_info), intent(in) :: dminfo
      type (mesh_type), intent(inout) :: grid
      type (state_type), intent(inout) :: state
      type (diag_type), intent(inout) :: diag
      integer, intent(in) :: test_case

      real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zx, zz, hx, cqw
      real (kind=RKIND), dimension(:,:), pointer :: ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, dss, t, rt, u, ru
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnEdge 
      integer, dimension(:,:), pointer :: edgesOnEdge
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, nz1, itr, cell1, cell2, nCellsSolve
      integer :: index_qv

      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: znu, znw, znwc, znwv
      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: znuc, znuv

      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: zc, zw, ah
      real (kind=RKIND), dimension(grid % nVertLevels ) :: zu, dzw, rdzwp, rdzwm

      real (kind=RKIND), dimension(grid % nVertLevels, grid % nCells) :: rh, thi, tbi, cqwb

      real (kind=RKIND) ::  r, xnutr
      real (kind=RKIND) ::  ztemp, zd, zt, dz, str

      real (kind=RKIND), dimension(grid % nVertLevels ) :: qvb
      real (kind=RKIND), dimension(grid % nVertLevels ) :: t_init_1d

      real (kind=RKIND) :: d1, d2, d3, cf1, cf2, cf3, cof1, cof2
      real (kind=RKIND) :: ztr, thetar, ttr, thetas, um, us, zts, pitop, pibtop, ptopb, ptop, rcp, rcv, p0
      real (kind=RKIND) :: radx, radz, zcent, xmid, delt, xloc, rad, yloc, ymid, a_scale
      real (kind=RKIND) :: pres, temp, es, qvs

      !
      ! Scale all distances
      !

      a_scale = 1.0

      grid % xCell % array = grid % xCell % array * a_scale
      grid % yCell % array = grid % yCell % array * a_scale
      grid % zCell % array = grid % zCell % array * a_scale
      grid % xVertex % array = grid % xVertex % array * a_scale
      grid % yVertex % array = grid % yVertex % array * a_scale
      grid % zVertex % array = grid % zVertex % array * a_scale
      grid % xEdge % array = grid % xEdge % array * a_scale
      grid % yEdge % array = grid % yEdge % array * a_scale
      grid % zEdge % array = grid % zEdge % array * a_scale
      grid % dvEdge % array = grid % dvEdge % array * a_scale
      grid % dcEdge % array = grid % dcEdge % array * a_scale
      grid % areaCell % array = grid % areaCell % array * a_scale**2.0
      grid % areaTriangle % array = grid % areaTriangle % array * a_scale**2.0
      grid % kiteAreasOnVertex % array = grid % kiteAreasOnVertex % array * a_scale**2.0

      weightsOnEdge     => grid % weightsOnEdge % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array
      
      nz1 = grid % nVertLevels
      nz = nz1 + 1
      nCellsSolve = grid % nCellsSolve

      zgrid => grid % zgrid % array
      rdzw => grid % rdzw % array
      dzu => grid % dzu % array
      rdzu => grid % rdzu % array
      fzm => grid % fzm % array
      fzp => grid % fzp % array
      zx => grid % zx % array
      zz => grid % zz % array
      hx => grid % hx % array
      dss => grid % dss % array

      ppb => diag % pressure_base % array
      pb => diag % exner_base % array
      rb => diag % rho_base % array
      tb => diag % theta_base % array
      rtb => diag % rtheta_base % array
      p => diag % exner % array
      cqw => diag % cqw % array

      rho_zz => state % rho_zz % array

      pp => diag % pressure_p % array
      rr => diag % rho_p % array
      t => state % theta_m % array      
      rt => diag % rtheta_p % array
      u => state % u % array
      ru => diag % ru % array

      scalars => state % scalars % array

      index_qv = state % index_qv

      scalars(:,:,:) = 0.

      call atm_initialize_advection_rk(grid) 
      call atm_initialize_deformation_weights(grid) 

      xnutr = 0.
      zd = 12000.

      p0 = 1.e+05
      rcp = rgas/cp
      rcv = rgas/(cp-rgas)

     write(0,*) ' point 1 in test case setup '

! We may pass in an hx(:,:) that has been precomputed elsewhere.
! For now it is independent of k

      do iCell=1,grid % nCells
        do k=1,nz
          hx(k,iCell) = 0.  ! squall line or supercell on flat plane
        enddo
      enddo

      !     metrics for hybrid coordinate and vertical stretching

      str = 1.0
      zt = 20000.
      dz = zt/float(nz1)

!      write(0,*) ' dz = ',dz
      write(0,*) ' hx computation complete '

      do k=1,nz
		
!           sh(k) is the stretching specified for height surfaces

            zc(k) = zt*(real(k-1)*dz/zt)**str 
				
!           to specify specific heights zc(k) for coordinate surfaces,
!           input zc(k) 
!           zw(k) is the hieght of zeta surfaces
!                zw(k) = (k-1)*dz yields constant dzeta
!                        and nonconstant dzeta/dz
!                zw(k) = sh(k)*zt yields nonconstant dzeta
!                        and nearly constant dzeta/dz 

!            zw(k) = float(k-1)*dz
            zw(k) = zc(k)
!
!           ah(k) governs the transition between terrain-following 
!           and pureheight coordinates
!                ah(k) = 0 is a terrain-following coordinate
!                ah(k) = 1 is a height coordinate
 
!            ah(k) = 1.-cos(.5*pii*(k-1)*dz/zt)**6
            ah(k) = 1.
!	    write(0,*) ' k, zc, zw, ah ',k,zc(k),zw(k),ah(k)			
      end do
      do k=1,nz1
         dzw (k) = zw(k+1)-zw(k)
         rdzw(k) = 1./dzw(k)
         zu(k  ) = .5*(zw(k)+zw(k+1))
      end do
      do k=2,nz1
         dzu (k)  = .5*(dzw(k)+dzw(k-1))
         rdzu(k)  =  1./dzu(k)
         fzp (k)  = .5* dzw(k  )/dzu(k)
         fzm (k)  = .5* dzw(k-1)/dzu(k)
         rdzwp(k) = dzw(k-1)/(dzw(k  )*(dzw(k)+dzw(k-1)))
         rdzwm(k) = dzw(k  )/(dzw(k-1)*(dzw(k)+dzw(k-1)))
      end do

!**********  how are we storing cf1, cf2 and cf3?

      COF1 = (2.*DZU(2)+DZU(3))/(DZU(2)+DZU(3))*DZW(1)/DZU(2) 
      COF2 =     DZU(2)        /(DZU(2)+DZU(3))*DZW(1)/DZU(3) 
      CF1  = FZP(2) + COF1
      CF2  = FZM(2) - COF1 - COF2
      CF3  = COF2       

!      d1  = .5*dzw(1)
!      d2  = dzw(1)+.5*dzw(2)
!      d3  = dzw(1)+dzw(2)+.5*dzw(3)
!      cf1 = d2*d3*(d3-d2)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf2 = d1*d3*(d1-d3)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf3 = d1*d2*(d2-d1)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))

      grid % cf1 % scalar = cf1
      grid % cf2 % scalar = cf2
      grid % cf3 % scalar = cf3

      do iCell=1,grid % nCells
        do k=1,nz	
            zgrid(k,iCell) = ah(k)*(zc(k)*(1.-hx(k,iCell)/zt)+hx(k,iCell)) &
                           + (1.-ah(k)) * zc(k)	
        end do
        do k=1,nz1
          zz (k,iCell) = (zw(k+1)-zw(k))/(zgrid(k+1,iCell)-zgrid(k,iCell))
        end do
      end do

      do i=1, grid % nEdges
        iCell1 = grid % CellsOnEdge % array(1,i)
        iCell2 = grid % CellsOnEdge % array(2,i)
        do k=1,nz
          zx (k,i) = (zgrid(k,iCell2)-zgrid(k,iCell1)) / grid % dcEdge % array(i)
        end do
      end do
      do i=1, grid % nCells
        do k=1,nz1
          ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
          dss(k,i) = 0.
          ztemp = zgrid(k,i)
          if(ztemp.gt.zd+.1)  then
             dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
          end if
        end do
      enddo

!
! convective initialization
!
         ztr    = 12000.
         thetar = 343.
         ttr    = 213.
         thetas = 300.5

!         write(0,*) ' rgas, cp, gravity ',rgas,cp, gravity

      if ( config_test_case == 4) then ! squall line parameters
         um = 12.
         us = 10.
         zts = 2500.
      else if (config_test_case == 5) then !supercell parameters
         um = 30.
         us = 15.
         zts = 5000.
      end if

         do i=1,grid % nCells
            do k=1,nz1
               ztemp = .5*(zgrid(k,i)+zgrid(k+1,i))
               if(ztemp .gt. ztr) then
                  t (k,i) = thetar*exp(9.8*(ztemp-ztr)/(1003.*ttr))
                  rh(k,i) = 0.25
               else
                  t (k,i) = 300.+43.*(ztemp/ztr)**1.25
                  rh(k,i) = (1.-0.75*(ztemp/ztr)**1.25)
                  if(t(k,i).lt.thetas) t(k,i) = thetas
               end if
               tb(k,i) = t(k,i)
               thi(k,i) = t(k,i)
               tbi(k,i) = t(k,i)
               cqw(k,i) = 1.
               cqwb(k,i) = 1.
            end do
         end do

!         rh(:,:) = 0.

!  set the velocity field - we are on a plane here.

         do i=1, grid % nEdges
            cell1 = grid % CellsOnEdge % array(1,i)
            cell2 = grid % CellsOnEdge % array(2,i)
            if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
            do k=1,nz1
               ztemp = .25*( zgrid(k,cell1 )+zgrid(k+1,cell1 )  &
                            +zgrid(k,cell2)+zgrid(k+1,cell2))
               if(ztemp.lt.zts)  then
                  u(k,i) = um*ztemp/zts
               else
                  u(k,i) = um
               end if
               if(i == 1 ) grid % u_init % array(k) = u(k,i) - us
               u(k,i) = cos(grid % angleEdge % array(i)) * (u(k,i) - us)
            end do
            end if
         end do

         call mpas_dmpar_bcast_reals(dminfo, nz1, grid % u_init % array)

!
!    for reference sounding 
!
     do itr=1,30

      pitop = 1.-.5*dzw(1)*gravity*(1.+scalars(index_qv,1,1))/(cp*t(1,1)*zz(1,1))
      pibtop = 1.-.5*dzw(1)*gravity*(1.+qvb(1))/(cp*tb(1,1)*zz(1,1))
      do k=2,nz1
         pitop = pitop-dzu(k)*gravity/(cp*cqw(k,1)*.5*(t(k,1)+t(k-1,1))   &
                                   *.5*(zz(k,1)+zz(k-1,1)))
         pibtop = pibtop-dzu(k)*gravity/(cp*cqwb(k,1)*.5*(tb(k,1)+tb(k-1,1))   &
                                   *.5*(zz(k,1)+zz(k-1,1)))

         !write(0,*) k,pitop,tb(k,1),dzu(k),tb(k,1)
      end do
      pitop = pitop-.5*dzw(nz1)*gravity*(1.+scalars(index_qv,nz1,1))/(cp*t(nz1,1)*zz(nz1,1))
      pibtop = pibtop-.5*dzw(nz1)*gravity*(1.+qvb(nz1))/(cp*tb(nz1,1)*zz(nz1,1))

      call mpas_dmpar_bcast_real(dminfo, pitop)
      call mpas_dmpar_bcast_real(dminfo, pibtop)

      ptopb = p0*pibtop**(1./rcp)
      write(6,*) 'ptopb = ',.01*ptopb

      do i=1, grid % nCells
         pb(nz1,i) = pibtop+.5*dzw(nz1)*gravity*(1.+qvb(nz1))/(cp*tb(nz1,i)*zz(nz1,i))
         p (nz1,i) = pitop+.5*dzw(nz1)*gravity*(1.+scalars(index_qv,nz1,i))/(cp*t (nz1,i)*zz(nz1,i))
         do k=nz1-1,1,-1
            pb(k,i)  = pb(k+1,i) + dzu(k+1)*gravity/(cp*cqwb(k+1,i)*.5*(tb(k,i)+tb(k+1,i))   &
                                           *.5*(zz(k,i)+zz(k+1,i)))
            p (k,i)  = p (k+1,i) + dzu(k+1)*gravity/(cp*cqw(k+1,i)*.5*(t (k,i)+t (k+1,i))   &
                                           *.5*(zz(k,i)+zz(k+1,i)))
         end do
         do k=1,nz1
            rb (k,i) = pb(k,i)**(1./rcv)/((rgas/p0)*tb(k,i)*zz(k,i))
            rtb(k,i) = rb(k,i)*tb(k,i)
            rr (k,i) = p (k,i)**(1./rcv)/((rgas/p0)*t (k,i)*zz(k,i))-rb(k,i)
            ppb(k,i) = p0*(zz(k,i)*rgas*rtb(k,i)/p0)**(cp/cv)
         end do
      end do

     !
     ! update water vapor mixing ratio from humidity profile
     !
      do i= 1,grid%nCells
         do k=1,nz1
            temp     = p(k,i)*thi(k,i)
            pres     = p0*p(k,i)**(1./rcp)
            qvs      = 380.*exp(17.27*(temp-273.)/(temp-36.))/pres
            scalars(index_qv,k,i) = min(0.014_RKIND,rh(k,i)*qvs)
         end do
      end do

      do k=1,nz1
!*********************************************************************
!           QVB = QV INCLUDES MOISTURE IN REFERENCE STATE
!            qvb(k) = scalars(index_qv,k,1)
                                        
!           QVB = 0 PRODUCES DRY REFERENCE STATE
            qvb(k) = 0.
!*********************************************************************
      end do

      do i= 1,grid%nCells
         do k=1,nz1
            t (k,i) = thi(k,i)*(1.+1.61*scalars(index_qv,k,i))
            tb(k,i) = tbi(k,i)*(1.+1.61*qvb(k))
         end do
         do k=2,nz1
            cqw (k,i) = 1./(1.+.5*(scalars(index_qv,k,i)+scalars(index_qv,k-1,i)))
            cqwb(k,i) = 1./(1.+.5*(qvb(k)+qvb(k-1)))
         end do
      end do

      end do !end of iteration loop

      write(0,*) ' base state sounding '
      write(0,*) ' k,     pb,     rb,     tb,     rtb,     t,     rr,      p,    qvb'
      do k=1,grid%nVertLevels
         write (0,'(i2,8(2x,f19.15))') k,pb(k,1),rb(k,1),tb(k,1),rtb(k,1),t(k,1),rr(k,1),p(k,1),qvb(k)
      end do

!
!     potential temperature perturbation
!
!      delt = -10.
!      delt = -0.01
      delt = 3.
      radx  = 10000.
      radz  = 1500.
      zcent = 1500.

      if (config_test_case == 4) then          ! squall line prameters
         call mpas_dmpar_max_real(dminfo, maxval(grid % xCell % array(:)), xmid)
         xmid = xmid * 0.5
         ymid = 0.0          ! Not used for squall line
      else if (config_test_case == 5) then     ! supercell parameters
         call mpas_dmpar_max_real(dminfo, maxval(grid % xCell % array(:)), xmid)
         call mpas_dmpar_max_real(dminfo, maxval(grid % yCell % array(:)), ymid)
         xmid = xmid * 0.5
         ymid = ymid * 0.5
      end if

      do i=1, grid % nCells
        xloc = grid % xCell % array(i) - xmid
        if (config_test_case == 4) then 
           yloc = 0.                            !squall line setting
        else if (config_test_case == 5) then
           yloc = grid % yCell % array(i) - ymid !supercell setting
        end if

        do k = 1,nz1
          ztemp     = .5*(zgrid(k+1,i)+zgrid(k,i))
          rad =sqrt((xloc/radx)**2+(yloc/radx)**2+((ztemp-zcent)/radz)**2)
          if(rad.lt.1)  then
            thi(k,i) = thi(k,i) + delt*cos(.5*pii*rad)**2
          end if
           t (k,i) = thi(k,i)*(1.+1.61*scalars(index_qv,k,i))
        end do
      end do

      do itr=1,30

        pitop = 1.-.5*dzw(1)*gravity*(1.+scalars(index_qv,1,1))/(cp*t(1,1)*zz(1,1))
        do k=2,nz1
          pitop = pitop-dzu(k)*gravity/(cp*cqw(k,1)*.5*(t (k,1)+t (k-1,1)) &
                                                  *.5*(zz(k,1)+zz(k-1,1)))
        end do
        pitop = pitop - .5*dzw(nz1)*gravity*(1.+scalars(index_qv,nz1,1))/(cp*t(nz1,1)*zz(nz1,1))
        ptop = p0*pitop**(1./rcp)
        write(0,*) 'ptop  = ',.01*ptop, .01*ptopb

        call mpas_dmpar_bcast_real(dminfo, ptop)

        do i = 1, grid % nCells

          pp(nz1,i) = ptop-ptopb+.5*dzw(nz1)*gravity*   &
                       (rr(nz1,i)+(rr(nz1,i)+rb(nz1,i))*scalars(index_qv,nz1,i))
          do k=nz1-1,1,-1
!             pp(k,i) = pp(k+1,i)+.5*dzu(k+1)*gravity*                   &
!                            (rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*scalars(index_qv,k  ,i)  &
!                            +rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*scalars(index_qv,k+1,i))
               pp(k,i) = pp(k+1,i)+dzu(k+1)*gravity*(    &
                            fzm(k+1)*(rb(k+1,i)*(scalars(index_qv,k+1,i)-qvb(k+1))    &
                                     +rr(k+1,i)*(1.+scalars(index_qv,k+1,i)))         &
                           +fzp(k+1)*(rb(k  ,i)*(scalars(index_qv,k  ,i)-qvb(k))      &
                                     +rr(k  ,i)*(1.+scalars(index_qv,k  ,i))))
          end do
          if (itr==1.and.i==1) then
          do k=1,nz1
          write(0,*) "pp-check", pp(k,i) 
          end do
          end if
          do k=1,nz1
             rt(k,i) = (pp(k,i)/(rgas*zz(k,i))                   &
                     -rtb(k,i)*(p(k,i)-pb(k,i)))/p(k,i)       
             p (k,i) = (zz(k,i)*(rgas/p0)*(rtb(k,i)+rt(k,i)))**rcv
             rr(k,i) = (rt(k,i)-rb(k,i)*(t(k,i)-tb(k,i)))/t(k,i)
          end do

        end do ! loop over cells

      end do !  iteration loop
!----------------------------------------------------------------------
!
      do k=1,nz1
        grid % qv_init % array(k) = scalars(index_qv,k,1)
      end do

      t_init_1d(:) = t(:,1)
      call mpas_dmpar_bcast_reals(dminfo, nz1, t_init_1d)
      call mpas_dmpar_bcast_reals(dminfo, nz1, grid % qv_init % array)

      do i=1,grid % ncells
         do k=1,nz1
            grid % t_init % array(k,i) = t_init_1d(k)
            rho_zz(k,i) = rb(k,i)+rr(k,i)
         end do
      end do

      do i=1,grid % nEdges
        cell1 = grid % CellsOnEdge % array(1,i)
        cell2 = grid % CellsOnEdge % array(2,i)
        if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
          do k=1,nz1
            ru (k,i)  = 0.5*(rho_zz(k,cell1)+rho_zz(k,cell2))*u(k,i)    
          end do
        end if
      end do


      !
      !  we are assuming w and rw are zero for this initialization
      !  i.e., no terrain
      !
       diag % rw % array = 0.
       state % w % array = 0.

       grid % zb % array = 0.
       grid % zb3% array = 0.

      !
      ! Generate rotated Coriolis field
      !
      do iEdge=1,grid % nEdges
         grid % fEdge % array(iEdge) = 0.
      end do

      do iVtx=1,grid % nVertices
         grid % fVertex % array(iVtx) = 0.
      end do

      !
      ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
      !
      diag % v % array(:,:) = 0.0
      do iEdge = 1, grid%nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            if (eoe > 0) then
               do k = 1, grid%nVertLevels
                 diag % v % array(k,iEdge) = diag % v %array(k,iEdge) + weightsOnEdge(i,iEdge) * state % u % array(k, eoe)
              end do
            end if
         end do
      end do

     ! write(0,*) ' k,u_init, t_init, qv_init '
     ! do k=1,grid%nVertLevels
     !   write(0,'(i2,3(2x,f14.10)') k,grid % u_init % array(k),grid % t_init% array(k),grid % qv_init % array(k)
     ! end do

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            diag % rho % array(k,iCell) = state % rho_zz % array(k,iCell) * zz(k,iCell)
            diag % theta % array(k,iCell) = state % theta_m % array(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do

   end subroutine atm_test_case_squall_line


!----------------------------------------------------------------------------------------------------------


   subroutine atm_test_case_mtn_wave(grid, state, diag, test_case)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup baroclinic wave test case from Jablonowski and Williamson 2008 (QJRMS)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mesh_type), intent(inout) :: grid
      type (state_type), intent(inout) :: state
      type (diag_type), intent(inout) :: diag
      integer, intent(in) :: test_case

      real (kind=RKIND), parameter :: t0=288., hm=250.

      real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zx, zz, hx, cqw
      real (kind=RKIND), dimension(:,:), pointer :: ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, dss, t, rt, u, ru 
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars, deriv_two, zb, zb3

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnEdge 
      integer, dimension(:,:), pointer :: edgesOnEdge, CellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: dvEdge, AreaCell, xCell, yCell 
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, nz1, itr, itrp, cell1, cell2, nCellsSolve
      integer :: index_qv

      real (kind=RKIND) :: ptop, pitop, ptopb, p0, flux, d2fdx2_cell1, d2fdx2_cell2

      real (kind=RKIND) :: ztemp, zd, zt, dz, str

      real (kind=RKIND), dimension(grid % nVertLevels, grid % nCells) :: rh
      real (kind=RKIND) :: es, qvs, xnutr, ptemp
      integer :: iter

      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: zc, zw, ah
      real (kind=RKIND), dimension(grid % nVertLevels ) :: zu, dzw, rdzwp, rdzwm

      real (kind=RKIND) :: d1, d2, d3, cof1, cof2, cf1, cf2, cf3
      real (kind=RKIND) :: um, us,  rcp, rcv
      real (kind=RKIND) :: xmid, temp, pres, a_scale

      real (kind=RKIND) :: xi, xa, xc, xla, zinv, xn2, xn2m, xn2l, sm, dzh, dzht, dzmin, z_edge, z_edge3 

      integer, dimension(grid % nCells, 2) :: next_cell
      real (kind=RKIND),  dimension(grid % nCells) :: hxzt
      logical, parameter :: terrain_smooth = .false. 

      !
      ! Scale all distances
      !

      a_scale = 1.0

      grid % xCell % array = grid % xCell % array * a_scale
      grid % yCell % array = grid % yCell % array * a_scale
      grid % zCell % array = grid % zCell % array * a_scale
      grid % xVertex % array = grid % xVertex % array * a_scale
      grid % yVertex % array = grid % yVertex % array * a_scale
      grid % zVertex % array = grid % zVertex % array * a_scale
      grid % xEdge % array = grid % xEdge % array * a_scale
      grid % yEdge % array = grid % yEdge % array * a_scale
      grid % zEdge % array = grid % zEdge % array * a_scale
      grid % dvEdge % array = grid % dvEdge % array * a_scale
      grid % dcEdge % array = grid % dcEdge % array * a_scale
      grid % areaCell % array = grid % areaCell % array * a_scale**2.0
      grid % areaTriangle % array = grid % areaTriangle % array * a_scale**2.0
      grid % kiteAreasOnVertex % array = grid % kiteAreasOnVertex % array * a_scale**2.0

      weightsOnEdge     => grid % weightsOnEdge % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array  
      dvEdge            => grid % dvEdge % array
      AreaCell          => grid % AreaCell % array
      CellsOnEdge       => grid % CellsOnEdge % array
      deriv_two         => grid % deriv_two % array
      
      nz1 = grid % nVertLevels
      nz = nz1 + 1
      nCellsSolve = grid % nCellsSolve

      zgrid => grid % zgrid % array
      zb => grid % zb % array
      zb3 => grid % zb3 % array
      rdzw => grid % rdzw % array
      dzu => grid % dzu % array
      rdzu => grid % rdzu % array
      fzm => grid % fzm % array
      fzp => grid % fzp % array
      zx => grid % zx % array
      zz => grid % zz % array
      hx => grid % hx % array
      dss => grid % dss % array
 
      xCell => grid % xCell % array
      yCell => grid % yCell % array

      ppb => diag % pressure_base % array
      pb => diag % exner_base % array
      rb => diag % rho_base % array
      tb => diag % theta_base % array
      rtb => diag % rtheta_base % array
      p => diag % exner % array
      cqw => diag % cqw % array

      rho_zz => state % rho_zz % array

      pp => diag % pressure_p % array
      rr => diag % rho_p % array
      t => state % theta_m % array      
      rt => diag % rtheta_p % array
      u => state % u % array
      ru => diag % ru % array

      scalars => state % scalars % array

      index_qv = state % index_qv

      scalars(:,:,:) = 0.

      call atm_initialize_advection_rk(grid) 
      call atm_initialize_deformation_weights(grid) 

      xnutr = 0.1
      zd = 10500.

      p0 = 1.e+05
      rcp = rgas/cp
      rcv = rgas/(cp-rgas)

      ! for hx computation
      xa = 5000. !SHP - should be changed based on grid distance 
      xla = 4000.
      xc = maxval (grid % xCell % array(:))/2. 

      !     metrics for hybrid coordinate and vertical stretching
      str = 1.0
      zt = 21000.
      dz = zt/float(nz1)
!      write(0,*) ' dz = ',dz

      do k=1,nz
		
!           sh(k) is the stretching specified for height surfaces

            zc(k) = zt*(real(k-1)*dz/zt)**str 
				
!           to specify specific heights zc(k) for coordinate surfaces,
!           input zc(k) 
!           zw(k) is the hieght of zeta surfaces
!                zw(k) = (k-1)*dz yields constant dzeta
!                        and nonconstant dzeta/dz
!                zw(k) = sh(k)*zt yields nonconstant dzeta
!                        and nearly constant dzeta/dz 

!            zw(k) = float(k-1)*dz
            zw(k) = zc(k)
!
!           ah(k) governs the transition between terrain-following 
!           and pureheight coordinates
!                ah(k) = 0 is a terrain-following coordinate
!                ah(k) = 1 is a height coordinate
 
!            ah(k) = 1.-cos(.5*pii*(k-1)*dz/zt)**6
            ah(k) = 1.
!	    write(0,*) ' k, zc, zw, ah ',k,zc(k),zw(k),ah(k)			
      end do
      do k=1,nz1
         dzw (k) = zw(k+1)-zw(k)
         rdzw(k) = 1./dzw(k)
         zu(k  ) = .5*(zw(k)+zw(k+1))
      end do
      do k=2,nz1
         dzu (k)  = .5*(dzw(k)+dzw(k-1))
         rdzu(k)  =  1./dzu(k)
         fzp (k)  = .5* dzw(k  )/dzu(k)
         fzm (k)  = .5* dzw(k-1)/dzu(k)
         rdzwp(k) = dzw(k-1)/(dzw(k  )*(dzw(k)+dzw(k-1)))
         rdzwm(k) = dzw(k  )/(dzw(k-1)*(dzw(k)+dzw(k-1)))
      end do

!**********  how are we storing cf1, cf2 and cf3?

      d1  = .5*dzw(1)
      d2  = dzw(1)+.5*dzw(2)
      d3  = dzw(1)+dzw(2)+.5*dzw(3)
      !cf1 = d2*d3*(d3-d2)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
      !cf2 = d1*d3*(d1-d3)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
      !cf3 = d1*d2*(d2-d1)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))

      cof1 = (2.*dzu(2)+dzu(3))/(dzu(2)+dzu(3))*dzw(1)/dzu(2)
      cof2 =     dzu(2)        /(dzu(2)+dzu(3))*dzw(1)/dzu(3)
      cf1  = fzp(2) + cof1
      cf2  = fzm(2) - cof1 - cof2
      cf3  = cof2

      grid % cf1 % scalar = cf1
      grid % cf2 % scalar = cf2
      grid % cf3 % scalar = cf3

! setting for terrain
      do iCell=1,grid % nCells
         xi = grid % xCell % array(iCell)
         !====1. for pure cosine mountain
         ! if(abs(xi-xc).ge.2.*xa)  then
         !    hx(1,iCell) = 0.
         ! else
         !    hx(1,iCell) = hm*cos(.5*pii*(xi-xc)/(2.*xa))**2.
         ! end if

         !====2. for cosine mountain
         !if(abs(xi-xc).lt.xa)  THEN
         !     hx(1,iCell) = hm*cos(pii*(xi-xc)/xla)**2. *cos(.5*pii*(xi-xc)/xa )**2.
         ! else
         !    hx(1,iCell) = 0.
         ! end if

         !====3. for shock mountain 
         hx(1,iCell) = hm*exp(-((xi-xc)/xa)**2)*cos(pii*(xi-xc)/xla)**2.

         hx(nz,iCell) = zt

!***** SHP -> get the temporary point information for the neighbor cell ->> should be changed!!!!! 
         do i=1,grid % nCells 
            !option 1
            !IF(yCell(i).eq.yCell(iCell).and.xCell(i).eq.xCell(iCell)-sqrt(3.)*grid % dcEdge % array(1)) next_cell(iCell,1) = i 
            !IF(yCell(i).eq.yCell(iCell).and.xCell(i).eq.xCell(iCell)+sqrt(3.)*grid % dcEdge % array(1)) next_cell(iCell,2) = i 
            !option 2
            next_cell(iCell,1) = iCell - 8 ! note ny=4
            next_cell(iCell,2) = iCell + 8 ! note ny=4

            if (xCell(iCell).le. 3.*grid % dcEdge % array(1)) then
                next_cell(iCell,1) = 1
            else if (xCell(iCell).ge. maxval(xCell(:))-3.*grid % dcEdge % array(1)) then
                next_cell(iCell,2) = 1
            end if

         end do
      enddo
      
      write(0,*) ' hx computation complete '


! smoothing grid for the upper level >> but not propoer for parallel programing 
      dzmin=.7
      do k=2,nz1
         sm = .25*min((zc(k)-zc(k-1))/dz,1.0_RKIND)
         do i=1,grid % nCells
            hx(k,i) = hx(k-1,i)
         end do

         do iter = 1,20 !iteration for smoothing

            do i=1,grid % nCells
               hxzt(i) = hx(k,i) + sm*(hx(k,next_cell(i,2))-2.*hx(k,i)+hx(k,next_cell(i,1)))
            end do
            dzh = zc(k) - zc(k-1)
            do i=1,grid % nCells
               dzht = zc(k)+hxzt(i) - zc(k-1)-hx(k-1,i)
               if(dzht.lt.dzh)  dzh = dzht
            end do

            if(dzh.gt.dzmin*(zc(k)-zc(k-1)))  then
               do i=1,grid % nCells
                  hx(k,i) = hxzt(i)
               end do
            else 
               goto 99  !SHP - this algorithm should be changed
            end if

         end do !end of iteration for smoothing
99       write(0,*) "PASS-SHP"
      end do

      do iCell=1,grid % nCells
        do k=1,nz
            if (terrain_smooth) then
            zgrid(k,iCell) = ah(k)*(zc(k)*(1.-hx(k,iCell)/zt)+hx(k,iCell)) &
                           + (1.-ah(k)) * zc(k)
            else
            zgrid(k,iCell) = ah(k)*(zc(k)*(1.-hx(1,iCell)/zt)+hx(1,iCell)) &
                           + (1.-ah(k)) * zc(k)
            end if
        end do
        do k=1,nz1
          zz (k,iCell) = (zw(k+1)-zw(k))/(zgrid(k+1,iCell)-zgrid(k,iCell))
        end do
      end do

      do i=1, grid % nEdges
        iCell1 = grid % CellsOnEdge % array(1,i)
        iCell2 = grid % CellsOnEdge % array(2,i)
        do k=1,nz
          zx (k,i) = (zgrid(k,iCell2)-zgrid(k,iCell1)) / grid % dcEdge % array(i)
        end do
      end do
      do i=1, grid % nCells
        do k=1,nz1
          ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
          dss(k,i) = 0.
          ztemp = zgrid(k,i)
          if(ztemp.gt.zd+.1)  then
             dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
          end if
        end do
      enddo

      write(0,*) ' grid metrics setup complete '

!
! mountain wave initialization
!
         !SHP-original
         !zinv = 1000.
         !SHP-schar case
         zinv = 3000.

         xn2  = 0.0001
         xn2m = 0.0000
         xn2l = 0.0001

         um = 10.
         us = 0.

         do i=1,grid % nCells
            do k=1,nz1
               ztemp   = .5*(zgrid(k,i)+zgrid(k+1,i))
               tb(k,i) =  t0*(1. + xn2m/gravity*ztemp) 
               if(ztemp .le. zinv) then
                  t (k,i) = t0*(1.+xn2l/gravity*ztemp)
               else
                  t (k,i) = t0*(1.+xn2l/gravity*zinv+xn2/gravity*(ztemp-zinv)) 
               end if
                  rh(k,i) = 0. 
            end do
         end do

!  set the velocity field - we are on a plane here.

         do i=1, grid % nEdges
            cell1 = grid % CellsOnEdge % array(1,i)
            cell2 = grid % CellsOnEdge % array(2,i)
            if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
            do k=1,nz1
               ztemp = .25*( zgrid(k,cell1 )+zgrid(k+1,cell1 )  &
                            +zgrid(k,cell2)+zgrid(k+1,cell2))
               u(k,i) = um
               if(i == 1 ) grid % u_init % array(k) = u(k,i) - us
               u(k,i) = cos(grid % angleEdge % array(i)) * (u(k,i) - us) 
            end do
            end if
         end do

!
!     reference sounding based on dry atmosphere
!
      pitop = 1.-.5*dzw(1)*gravity/(cp*tb(1,1)*zz(1,1))
      do k=2,nz1
         pitop = pitop-dzu(k)*gravity/(cp*(fzm(k)*tb(k,1)+fzp(k)*tb(k-1,1))   &
                                         *(fzm(k)*zz(k,1)+fzp(k)*zz(k-1,1)))
      end do
      pitop = pitop-.5*dzw(nz1)*gravity/(cp*tb(nz1,1)*zz(nz1,1))
      ptopb = p0*pitop**(1./rcp)
                
      do i=1, grid % nCells
         pb(nz1,i) = pitop+.5*dzw(nz1)*gravity/(cp*tb(nz1,i)*zz(nz1,i))
         p (nz1,i) = pitop+.5*dzw(nz1)*gravity/(cp*t (nz1,i)*zz(nz1,i))
         do k=nz1-1,1,-1
            pb(k,i)  = pb(k+1,i) + dzu(k+1)*gravity/(cp*.5*(tb(k,i)+tb(k+1,i))   &
                                           *.5*(zz(k,i)+zz(k+1,i)))
            p (k,i)  = p (k+1,i) + dzu(k+1)*gravity/(cp*.5*(t (k,i)+t (k+1,i))   &
                                           *.5*(zz(k,i)+zz(k+1,i)))
         end do
         do k=1,nz1
            rb (k,i) = pb(k,i)**(1./rcv)/((rgas/p0)*tb(k,i)*zz(k,i))
            rtb(k,i) = rb(k,i)*tb(k,i)
            rr (k,i) = p (k,i)**(1./rcv)/((rgas/p0)*t (k,i)*zz(k,i))-rb(k,i)
            cqw(k,i) = 1.
         end do
      end do

       write(0,*) ' ***** base state sounding ***** '
       write(0,*) 'k       pb        p         rb         rtb         rr          tb          t'
       do k=1,grid%nVertLevels
          write(0,'(i2,7(2x,f14.9))') k,pb(k,1),p(k,1),rb(k,1),rtb(k,1),rr(k,1),tb(k,1),t(k,1)
       end do
 
       scalars(index_qv,:,:) = 0.

!-------------------------------------------------------------------
!     ITERATIONS TO CONVERGE MOIST SOUNDING
      do itr=1,30
        pitop = 1.-.5*dzw(1)*gravity*(1.+scalars(index_qv,1,1))/(cp*t(1,1)*zz(1,1))

        do k=2,nz1
          pitop = pitop-dzu(k)*gravity/(cp*cqw(k,1)*(fzm(k)*t (k,1)+fzp(k)*t (k-1,1)) &
                                                   *(fzm(k)*zz(k,1)+fzp(k)*zz(k-1,1)))
        end do
        pitop = pitop - .5*dzw(nz1)*gravity*(1.+scalars(index_qv,nz1,1))/(cp*t(nz1,1)*zz(nz1,1))
        ptop = p0*pitop**(1./rcp)

        do i = 1, grid % nCells

           pp(nz1,i) = ptop-ptopb+.5*dzw(nz1)*gravity*   &
                       (rr(nz1,i)+(rr(nz1,i)+rb(nz1,i))*scalars(index_qv,nz1,i))
           do k=nz1-1,1,-1
              pp(k,i) = pp(k+1,i)+dzu(k+1)*gravity*                   &
                            (fzm(k)*(rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*scalars(index_qv,k  ,i))  &
                            +fzp(k)*(rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*scalars(index_qv,k+1,i)))
           end do
           do k=1,nz1
              rt(k,i) = (pp(k,i)/(rgas*zz(k,i))                   &
                      -rtb(k,i)*(p(k,i)-pb(k,i)))/p(k,i)
              p (k,i) = (zz(k,i)*(rgas/p0)*(rtb(k,i)+rt(k,i)))**rcv
              rr(k,i) = (rt(k,i)-rb(k,i)*(t(k,i)-tb(k,i)))/t(k,i)
           end do
!
!     update water vapor mixing ratio from humitidty profile
!
           do k=1,nz1
              temp   = p(k,i)*t(k,i)
              pres   = p0*p(k,i)**(1./rcp)
              qvs    = 380.*exp(17.27*(temp-273.)/(temp-36.))/pres
              scalars(index_qv,k,i) = min(0.014_RKIND,rh(k,i)*qvs)
           end do
                         
           do k=1,nz1
              t (k,i) = t(k,i)*(1.+1.61*scalars(index_qv,k,i))
           end do
           do k=2,nz1
              cqw(k,i) = 1./(1.+.5*( scalars(index_qv,k-1,i)  &
                                    +scalars(index_qv,k  ,i)))
           end do

        end do ! loop over cells

      end do !  iteration loop
!----------------------------------------------------------------------
!
      write(0,*) ' *** sounding for the simulation ***'
      write(0,*) '    z       theta       pres         qv       rho_m        u        rr'
      do k=1,nz1
         write(6,'(8(f14.9,2x))') .5*(zgrid(k,1)+zgrid(k+1,1))/1000.,   &
                       t(k,1)/(1.+1.61*scalars(index_qv,k,1)),        &
                       .01*p0*p(k,1)**(1./rcp),                       &
                       1000.*scalars(index_qv,k,1),                   &
                       (rb(k,1)+rr(k,1))*(1.+scalars(index_qv,k,1)),  &
                       grid % u_init % array(k), rr(k,1)
      end do

      do i=1,grid % ncells
         do k=1,nz1
            rho_zz(k,i) = rb(k,i)+rr(k,i)
         end do

        do k=1,nz1
            grid % t_init % array(k,i) = t(k,i)
        end do
      end do

      do i=1,grid % nEdges
        cell1 = grid % CellsOnEdge % array(1,i)
        cell2 = grid % CellsOnEdge % array(2,i)
        if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
          do k=1,nz1
            ru (k,i)  = 0.5*(rho_zz(k,cell1)+rho_zz(k,cell2))*u(k,i)    
          end do
        end if
      end do

!
!     pre-calculation z-metric terms in omega eqn.
!
      do iEdge = 1,grid % nEdges
         cell1 = CellsOnEdge(1,iEdge)
         cell2 = CellsOnEdge(2,iEdge)
         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then

            do k = 1, grid%nVertLevels

               if (config_theta_adv_order == 2) then

                  z_edge = (zgrid(k,cell1)+zgrid(k,cell2))/2.

               else !theta_adv_order == 3 or 4 

                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * zgrid(k,cell1)
                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * zgrid(k,cell2)
                  do i=1, grid % nEdgesOnCell % array (cell1)
                     if ( grid % CellsOnCell % array (i,cell1) > 0)       &
                     d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * zgrid(k,grid % CellsOnCell % array (i,cell1))
                  end do
                  do i=1, grid % nEdgesOnCell % array (cell2)
                     if ( grid % CellsOnCell % array (i,cell2) > 0)       &
                     d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * zgrid(k,grid % CellsOnCell % array (i,cell2))
                  end do             
             
                  z_edge =  0.5*(zgrid(k,cell1) + zgrid(k,cell2))         &
                                - (grid % dcEdge % array(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12. 

                  if (config_theta_adv_order == 3) then
                     z_edge3 =  - (grid % dcEdge % array(iEdge) **2) * (d2fdx2_cell1 - d2fdx2_cell2) / 12.   
                  else 
                     z_edge3 = 0.
                  end if

               end if

                  zb(k,1,iEdge) = (z_edge-zgrid(k,cell1))*dvEdge(iEdge)/AreaCell(cell1) 
                  zb(k,2,iEdge) = (z_edge-zgrid(k,cell2))*dvEdge(iEdge)/AreaCell(cell2) 
                  zb3(k,1,iEdge)=  z_edge3*dvEdge(iEdge)/AreaCell(cell1) 
                  zb3(k,2,iEdge)=  z_edge3*dvEdge(iEdge)/AreaCell(cell2) 
  
            end do

         end if
       end do

!     for including terrain
      state % w % array(:,:) = 0.0
      diag % rw % array(:,:) = 0.0

!
!     calculation of omega, rw = zx * ru + zz * rw
!

      do iEdge = 1,grid % nEdges

         cell1 = CellsOnEdge(1,iEdge)
         cell2 = CellsOnEdge(2,iEdge)

         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then
         do k = 2, grid%nVertLevels
            flux =  (fzm(k)*ru(k,iEdge)+fzp(k)*ru(k-1,iEdge))  
            diag % rw % array(k,cell2) = diag % rw % array(k,cell2) + (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb(k,2,iEdge)*flux 
            diag % rw % array(k,cell1) = diag % rw % array(k,cell1) - (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb(k,1,iEdge)*flux 

            if (config_theta_adv_order ==3) then
               diag % rw % array(k,cell2) = diag % rw % array(k,cell2)    &
                                            - sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb3(k,2,iEdge)*flux
               diag % rw % array(k,cell1) = diag % rw % array(k,cell1)    &
                                            + sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb3(k,1,iEdge)*flux
            end if

         end do
         end if

      end do

      ! Compute w from rho_zz and rw
      do iCell=1,grid%nCells
         do k=2,grid%nVertLevels
            state % w % array(k,iCell) = diag % rw % array(k,iCell)     & 
                                       / (fzp(k) * state % rho_zz % array(k-1,iCell) + fzm(k) * state % rho_zz % array(k,iCell))
         end do
      end do


      do iEdge=1,grid % nEdges
         grid % fEdge % array(iEdge) = 0.
      end do

      do iVtx=1,grid % nVertices
         grid % fVertex % array(iVtx) = 0.
      end do

      !
      ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
      !
      diag % v % array(:,:) = 0.0
      do iEdge = 1, grid%nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            if (eoe > 0) then
               do k = 1, grid%nVertLevels
                 diag % v % array(k,iEdge) = diag % v %array(k,iEdge) + weightsOnEdge(i,iEdge) * state % u % array(k, eoe)
              end do
            end if
         end do
      end do

!      do k=1,grid%nVertLevels
!        write(0,*) ' k,u_init, t_init, qv_init ',k,grid % u_init % array(k),grid % t_init% array(k),grid % qv_init % array(k)
!      end do

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            diag % rho % array(k,iCell) = state % rho_zz % array(k,iCell) * zz(k,iCell)
            diag % theta % array(k,iCell) = state % theta_m % array(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do

   end subroutine atm_test_case_mtn_wave


!----------------------------------------------------------------------------------------------------------

   real (kind=RKIND) function sphere_distance(lat1, lon1, lat2, lon2, radius)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on a
   !   sphere with given radius.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius

      real (kind=RKIND) :: arg1

      arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
                   cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
      sphere_distance = 2.*radius*asin(arg1)

   end function sphere_distance

!--------------------------------------------------------------------
   real (kind=RKIND) function env_qv( z, temperature, pressure, rh_max )

      implicit none
      real (kind=RKIND) :: z, temperature, pressure, ztr, es, qvs, p0, rh_max

      p0 = 100000.

!      ztr = 5000.
!
!      if(z .gt. ztr) then
!         env_qv = 0.
!      else
!         if(z.lt.2000.) then
!            env_qv = .5
!         else
!            env_qv = .5*(1.-(z-2000.)/(ztr-2000.))
!         end if
!      end if

       if (pressure .lt. 50000. ) then
           env_qv = 0.0
       else
           env_qv = (1.-((p0-pressure)/50000.)**1.25)
       end if

       env_qv = min(rh_max,env_qv)

! env_qv is the relative humidity, turn it into mixing ratio
       if (temperature .gt. 273.15) then
           es  = 1000.*0.6112*exp(17.67*(temperature-273.15)/(temperature-29.65))
       else
           es  = 1000.*0.6112*exp(21.8745584*(temperature-273.16)/(temperature-7.66))
       end if
       qvs = (287.04/461.6)*es/(pressure-es)

       ! qvs =  380.*exp(17.27*(temperature-273.)/(temperature-36.))/pressure

        env_qv = env_qv*qvs

   end function env_qv

end module atm_test_cases
