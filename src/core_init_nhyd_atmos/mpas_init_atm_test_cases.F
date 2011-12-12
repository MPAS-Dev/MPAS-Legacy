module init_atm_test_cases

   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_dmpar
   use atm_advection
   use mpas_sort
   use mpas_timekeeping

   use mpas_atmphys_initialize_real


   contains


   subroutine init_atm_setup_test_case(domain)
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
            call init_atm_test_case_jw(block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % diag, config_test_case)
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
            call init_atm_test_case_squall_line(domain % dminfo, block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % diag, config_test_case)
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
            call init_atm_test_case_mtn_wave(block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % diag, config_test_case)
            write(0,*) ' returned from test case setup '
            do i=2,nTimeLevs
               call mpas_copy_state(block_ptr % state % time_levs(i) % state, block_ptr % state % time_levs(1) % state)
            end do

            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 7 ) then

         write(0,*) ' real-data GFS test case '
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call init_atm_test_case_gfs(domain % dminfo, block_ptr % mesh, block_ptr % fg, block_ptr % state % time_levs(1) % state, &
                                    block_ptr % diag, config_test_case, block_ptr % parinfo)
            if(config_physics_init) &
               call physics_initialize_real(block_ptr % mesh, block_ptr % fg)

            do i=2,nTimeLevs
               call mpas_copy_state(block_ptr % state % time_levs(i) % state, block_ptr % state % time_levs(1) % state)
            end do

            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 8 ) then

         write(0,*) ' real-data surface (SST) update test case '
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call init_atm_test_case_sst(domain, domain % dminfo, block_ptr % mesh, block_ptr % fg, block_ptr % state % time_levs(1) % state, &
                                    block_ptr % diag, config_test_case, block_ptr % parinfo)
            do i=2,nTimeLevs
               call mpas_copy_state(block_ptr % state % time_levs(i) % state, block_ptr % state % time_levs(1) % state)
            end do

            block_ptr => block_ptr % next
         end do

      else


         write(0,*) ' Only test case 1, 2, 3, 4, 5, 6, and 7 are currently supported for nonhydrostatic core '
         stop
      end if

   end subroutine init_atm_setup_test_case

!----------------------------------------------------------------------------------------------------------

   subroutine init_atm_test_case_jw(grid, state, diag, test_case)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup baroclinic wave test case from Jablonowski and Williamson 2008 (QJRMS)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mesh_type), intent(inout) :: grid
      type (state_type), intent(inout) :: state
      type (diag_type), intent(inout) :: diag
      integer, intent(in) :: test_case

      real (kind=RKIND), parameter :: u0 = 35.0
      real (kind=RKIND), parameter :: alpha_grid = 0.  ! no grid rotation
      real (kind=RKIND), parameter :: omega_e = 7.29212e-05
      real (kind=RKIND), parameter :: t0b = 250., t0 = 288., delta_t = 4.8e+05, dtdz = 0.005, eta_t = 0.2
      real (kind=RKIND), parameter :: u_perturbation = 1., pert_radius = 0.1, latitude_pert = 40., longitude_pert = 20.
      real (kind=RKIND), parameter :: theta_c = pii/4.0
      real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
      real (kind=RKIND), parameter :: rh_max = 0.4       ! Maximum relative humidity
      real (kind=RKIND), parameter :: k_x = 9.           ! Normal mode wave number

      real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zx, zz, hx
      real (kind=RKIND), dimension(:,:), pointer :: pressure, ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, dss, t, rt
      real (kind=RKIND), dimension(:,:,:), pointer :: zf, zf3, zb, zb3
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two

      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, nz1, itr, itrp, cell1, cell2, nCellsSolve

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnEdge 
      integer, dimension(:,:), pointer :: edgesOnEdge, CellsOnEdge
      real, dimension(:), pointer :: dvEdge, AreaCell 
      real, dimension(:,:), pointer :: weightsOnEdge

      real (kind=RKIND) :: u, v, flux, fluxk, lat1, lat2, eta_v, r_pert, u_pert, lat_pert, lon_pert, r

      real (kind=RKIND) :: ptop, p0, phi
      real (kind=RKIND) :: lon_Edge

      real (kind=RKIND) :: r_earth, etavs, ztemp, zd, zt, dz, gam, delt, str

      real (kind=RKIND), dimension(grid % nVertLevels, grid % nCells) :: rel_hum, temperature, qv
      real (kind=RKIND) :: ptmp, es, qvs, xnutr, znut, ptemp
      integer :: iter

      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: hyai, hybi, znu, znw, znwc, znwv, hyam, hybm
      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: znuc, znuv, bn, divh, dpn

      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: sh, zw, ah
      real (kind=RKIND), dimension(grid % nVertLevels ) :: zu, dzw, rdzwp, rdzwm
      real (kind=RKIND), dimension(grid % nVertLevels ) :: eta, etav, teta, ppi, tt

      real (kind=RKIND) :: d1, d2, d3, cf1, cf2, cf3, cof1, cof2, psurf

      !  storage for (lat,z) arrays for zonal velocity calculation

      integer, parameter :: nlat=361
      real (kind=RKIND), dimension(grid % nVertLevels + 1) :: zz_1d, zgrid_1d, hx_1d
      real (kind=RKIND), dimension(grid % nVertLevels) :: flux_zonal
      real (kind=RKIND), dimension(nlat, grid % nVertLevels) :: u_2d, etavs_2d
      real (kind=RKIND), dimension(nlat) :: lat_2d
      real (kind=RKIND) :: dlat
      real (kind=RKIND) :: z_edge, z_edge3, d2fdx2_cell1, d2fdx2_cell2

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
      zf  => grid % zf % array
      zf3 => grid % zf3% array
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

      scalars => state % scalars % array

      scalars(:,:,:) = 0.

      call atm_initialize_advection_rk(grid) 
      call atm_initialize_deformation_weights(grid) 

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

      str = 1.5
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

      do k=1,nz1
        write(0,*) ' k, zgrid(k,1),hx(k,1) ',k,zgrid(k,1),hx(k,1)
      enddo

      do k=1,nz1
        write(0,*) ' k, zx(k,1) ',k,zx(k,1)
      enddo

      write(0,*) ' grid metrics setup complete '

!**************  section for 2d (lat,z) calc for zonal velocity

      dlat = 0.5*pii/float(nlat-1)
      do i = 1,nlat

        lat_2d(i) = float(i-1)*dlat
!        write(0,*) ' zonal setup, latitude = ',lat_2d(i)*180./pii

        do k=1,nz
          phi = lat_2d(i)
          hx_1d(k) = u0/gravity*cos(etavs)**1.5                            &
                      *((-2.*sin(phi)**6                                   &
                            *(cos(phi)**2+1./3.)+10./63.)                  &
                            *(u0)*cos(etavs)**1.5                          &
                       +(1.6*cos(phi)**3                                   &
                            *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)
        enddo

        do k=1,nz	
          zgrid_1d(k) = (1.-ah(k))*(sh(k)*(zt-hx_1d(k))+hx_1d(k))  &
                         + ah(k) * sh(k)* zt	
        end do
        do k=1,nz1
          zz_1d (k) = (zw(k+1)-zw(k))/(zgrid_1d(k+1)-zgrid_1d(k))
        end do

        do k=1,nz1
          ztemp    = .5*(zgrid_1d(k+1)+zgrid_1d(k))
          ppb(k,i) = p0*exp(-gravity*ztemp/(rgas*t0b)) 
          pb (k,i) = (ppb(k,i)/p0)**(rgas/cp)
          rb (k,i) = ppb(k,i)/(rgas*t0b*zz_1d(k))
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
          ! phi = grid % latCell % array (i)
          phi = lat_2d (i)
          do k=1,nz1
            tt(k) = 0.
            tt(k) = teta(k)+.75*eta(k)*pii*u0/rgas*sin(etav(k))      &
                            *sqrt(cos(etav(k)))*                   &
                              ((-2.*sin(phi)**6                    &
                                   *(cos(phi)**2+1./3.)+10./63.)   &
                                   *2.*u0*cos(etav(k))**1.5        &
                              +(1.6*cos(phi)**3                    &
                                *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)


            ztemp   = .5*(zgrid_1d(k)+zgrid_1d(k+1))
            ptemp   = ppb(k,i) + pp(k,i)
            qv(k,i) = 0.

          end do
		
          do itrp = 1,25
            do k=1,nz1				
              rr(k,i)  = (pp(k,i)/(rgas*zz_1d(k))  &
                          -rb(k,i)*(tt(k)-t0b))/tt(k)
            end do

            ppi(1) = p0-.5*dzw(1)*gravity                         &
                          *(1.25*(rr(1,i)+rb(1,i))*(1.+qv(1,i))   &
                            -.25*(rr(2,i)+rb(2,i))*(1.+qv(2,i)))

            ppi(1) = ppi(1)-ppb(1,i)
            do k=1,nz1-1
              ppi(k+1) = ppi(k)-.5*dzu(k+1)*gravity*                     &
                            (rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*qv(k  ,i)   &
                            +rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*qv(k+1,i))
            end do

            do k=1,nz1
              pp(k,i) = .2*ppi(k)+.8*pp(k,i)
            end do

          end do  ! end inner iteration loop itrp

        end do  ! end outer iteration loop itr

        do k=1,nz1
          etavs_2d(i,k) = (0.5*(ppb(k,i)+ppb(k,i)+pp(k,i)+pp(k,i))/p0 - 0.252)*pii/2.
!          u_2d(i,k) = u0*(sin(2.*lat_2d(i))**2) *(cos(etavs_2d(i,k))**1.5)
          u_2d(i,k) = u0*(sin(2.*lat_2d(i))**2) *(cos(etavs_2d(i,k))**1.5)*(rb(k,i)+rr(k,i))
        end do

      end do  ! end loop over latitudes for 2D zonal wind field calc

!      do i=1,nlat
!        do k=1,nz1
!          u_2d(i,k) = u_2d(i,k) - u0*(sin(2.*lat_2d(i))**2) *(cos(etavs_2d(nlat/2,k))**1.5)
!        end do
!      end do
!
!      write(22,*) nz1,nlat,u_2d

!******************************************************************      

!
!---- baroclinc wave initialization ---------------------------------
!
!     reference sounding based on dry isothermal atmosphere
!
      do i=1, grid % nCells
        !write(0,*) ' thermodynamic setup, cell ',i
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

        if(i == 1) then
          do k=1,nz1
            write(0,*) ' k, ppb, pb, rb, tb (k,1) ',k,ppb(k,1),pb(k,1),rb(k,1)*zz(k,1),tb(k,1)
          enddo
        end if
!
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
            tt(k) = 0.
            tt(k) = teta(k)+.75*eta(k)*pii*u0/rgas*sin(etav(k))      &
                            *sqrt(cos(etav(k)))*                   &
                              ((-2.*sin(phi)**6                    &
                                   *(cos(phi)**2+1./3.)+10./63.)   &
                                   *2.*u0*cos(etav(k))**1.5        &
                              +(1.6*cos(phi)**3                    &
                                *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)


            !write(0,*) ' k, tt(k) ',k,tt(k)
            ztemp   = .5*(zgrid(k,i)+zgrid(k+1,i))
            ptemp   = ppb(k,i) + pp(k,i)
!            qv(k,i) = env_qv( ztemp, tt(k), ptemp, 0 )
            qv(k,i) = 0.

          end do
!          do k=2,nz1
!            cqw(k,i) = 1./(1.+.5*(qv(k,i)+qv(k-1,i)))
!          end do
		
          do itrp = 1,25
            do k=1,nz1				
              rr(k,i)  = (pp(k,i)/(rgas*zz(k,i))  &
                          -rb(k,i)*(tt(k)-t0b))/tt(k)
            end do

            ppi(1) = p0-.5*dzw(1)*gravity                         &
                          *(1.25*(rr(1,i)+rb(1,i))*(1.+qv(1,i))   &
                            -.25*(rr(2,i)+rb(2,i))*(1.+qv(2,i)))

            ppi(1) = ppi(1)-ppb(1,i)
            do k=1,nz1-1
              ppi(k+1) = ppi(k)-.5*dzu(k+1)*gravity*                     &
                            (rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*qv(k  ,i)   &
                            +rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*qv(k+1,i))
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

        if(i == 1) then
          do k=1,nz1
            write(0,*) ' k, p, t, rt ',k,p(k,1),t(k,1),rt(k,1)
          enddo
        end if

      end do  ! end loop over cells

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
                                      lat_pert, lon_pert, 1.)/(pert_radius)
            u_pert = u_perturbation*exp(-r_pert**2)*(lat2-lat1)*a/grid % dvEdge % array(iEdge)

         else if (config_test_case == 3) then
            lon_Edge = grid % lonEdge % array(iEdge)
            u_pert = u_perturbation*cos(k_x*(lon_Edge - lon_pert)) &
                         *(0.5*(lat2-lat1) - 0.125*(sin(4.*lat2) - sin(4.*lat1)))*a/grid % dvEdge % array(iEdge)
         else
            u_pert = 0.0
         end if

         call init_atm_calc_flux_zonal(u_2d,etavs_2d,lat_2d,flux_zonal,lat1,lat2,grid % dvEdge % array(iEdge),a,u0,nz1,nlat)

         do k=1,grid % nVertLevels
!!           etavs = (0.5*(ppb(k,iCell1)+ppb(k,iCell2)+pp(k,iCell1)+pp(k,iCell2))/p0 - 0.252)*pii/2.
!           etavs = (0.5*(ppb(k,1)+ppb(k,1)+pp(k,1)+pp(k,1))/p0 - 0.252)*pii/2.
           etavs = (0.5*(ppb(k,440)+ppb(k,440)+pp(k,440)+pp(k,440))/p0 - 0.252)*pii/2.  ! 10262 mesh
!           etavs = (0.5*(ppb(k,505)+ppb(k,505)+pp(k,505)+pp(k,505))/p0 - 0.252)*pii/2.  ! 40962 mesh
  
!           fluxk = u0*flux*(cos(etavs)**1.5)

            fluxk = u0*flux_zonal(k)/(0.5*(rb(k,iCell1)+rb(k,iCell2)+rr(k,iCell1)+rr(k,iCell2)))

!           if(k.eq.18) then
!              write(21,*) ' iEdge, u1, u2 ',iEdge,fluxk,u0*flux_zonal(k)
!           end if
!!           fluxk = u0*flux*(cos(znuv(k))**(1.5))
!!           fluxk = u0 * cos(grid % angleEdge % array(iEdge)) * (sin(lat1+lat2)**2) *(cos(etavs)**1.5)
           state % u % array(k,iEdge) = fluxk + u_pert
         end do

         cell1 = grid % CellsOnEdge % array(1,i)
         cell2 = grid % CellsOnEdge % array(2,i)
         if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
           do k=1,nz1
             diag % ru % array (k,iEdge)  = 0.5*(rho_zz(k,cell1)+rho_zz(k,cell2))*state % u % array (k,iEdge)
           end do
         end if

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
         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then

            do k = 1, grid%nVertLevels

               if (config_theta_adv_order == 2) then

                  z_edge = (zgrid(k,cell1)+zgrid(k,cell2))/2.

               else if (config_theta_adv_order == 3 .or. config_theta_adv_order ==4) then !theta_adv_order == 3 or 4 

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

                  if (k /= 1) then
                     zf(k,1,iEdge) = ( fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb(k,1,iEdge)
                     zf(k,2,iEdge) = ( fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb(k,2,iEdge)
                     zf3(k,1,iEdge)= ( fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb3(k,1,iEdge)
                     zf3(k,2,iEdge)= ( fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb3(k,2,iEdge)
                  end if

            end do

         end if
       end do

      ! for including terrain
      diag % rw % array = 0.
      state % w % array = 0.
      do iEdge = 1,grid % nEdges

         cell1 = CellsOnEdge(1,iEdge)
         cell2 = CellsOnEdge(2,iEdge)

         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then
         do k = 2, grid%nVertLevels
            flux =  (fzm(k)*diag % ru % array(k,iEdge)+fzp(k)*diag % ru % array(k-1,iEdge))
            diag % rw % array(k,cell2) = diag % rw % array(k,cell2) + zf(k,2,iEdge)*flux
            diag % rw % array(k,cell1) = diag % rw % array(k,cell1) - zf(k,1,iEdge)*flux

            if (config_theta_adv_order ==3) then 
               diag % rw % array(k,cell2) = diag % rw % array(k,cell2)    &
                                            - sign(1.,diag % ru % array(k,iEdge))*config_coef_3rd_order*zf3(k,2,iEdge)*flux
               diag % rw % array(k,cell1) = diag % rw % array(k,cell1)    &
                                            + sign(1.,diag % ru % array(k,iEdge))*config_coef_3rd_order*zf3(k,1,iEdge)*flux
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

      do i=1,10
        psurf = (cf1*(ppb(1,i)+pp(1,i)) + cf2*(ppb(2,i)+pp(2,i)) + cf3*(ppb(3,i)+pp(3,i)))/100.

            psurf = (ppb(1,i)+pp(1,i)) + .5*dzw(1)*gravity        &
                          *(1.25*(rr(1,i)+rb(1,i))*(1.+qv(1,i))   &
                            -.25*(rr(2,i)+rb(2,i))*(1.+qv(2,i)))

        write(0,*) ' i, psurf, lat ',i,psurf,grid%latCell%array(i)*180./3.1415828
      enddo

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            diag % rho % array(k,iCell) = state % rho_zz % array(k,iCell) * zz(k,iCell)
            diag % theta % array(k,iCell) = state % theta_m % array(k,iCell) / (1.0 + 1.61 * scalars(state % index_qv,k,iCell))
         end do
      end do

   end subroutine init_atm_test_case_jw

   subroutine init_atm_calc_flux_zonal(u_2d,etavs_2d,lat_2d,flux_zonal,lat1_in,lat2_in,dvEdge,a,u0,nz1,nlat)

   implicit none
   integer, intent(in) :: nz1,nlat
   real (kind=RKIND), dimension(nlat,nz1), intent(in) :: u_2d,etavs_2d
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
       flux_zonal(k) = flux_zonal(k) + w1*u_2d(i,k) + w2*u_2d(i+1,k)
     end do

     end if

   end do

!  renormalize for setting cell-face fluxes

   do k=1,nz1
     flux_zonal(k) = sign(1.,lat2_in-lat1_in)*flux_zonal(k)*dlat*a/dvEdge/u0
   end do
     
   end subroutine init_atm_calc_flux_zonal


!----------------------------------------------------------------------------------------------------------

   subroutine init_atm_test_case_squall_line(dminfo, grid, state, diag, test_case)
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
      real, dimension(:,:), pointer :: weightsOnEdge

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
      write(0,*) 'ptopb = ',.01*ptopb

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
            scalars(index_qv,k,i) = amin1(0.014,rh(k,i)*qvs)
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

        call mpas_dmpar_bcast_real(dminfo, pitop)

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
          print *, "pp-check", pp(k,i) 
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

       grid % zf % array = 0.
       grid % zf3% array = 0.
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
            diag % theta % array(k,iCell) = state % theta_m % array(k,iCell) / (1.0 + 1.61 * scalars(state % index_qv,k,iCell))
         end do
      end do

   end subroutine init_atm_test_case_squall_line


!----------------------------------------------------------------------------------------------------------


   subroutine init_atm_test_case_mtn_wave(grid, state, diag, test_case)
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
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars, deriv_two, zf, zf3, zb, zb3

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnEdge 
      integer, dimension(:,:), pointer :: edgesOnEdge, CellsOnEdge
      real, dimension(:), pointer :: dvEdge, AreaCell, xCell, yCell 
      real, dimension(:,:), pointer :: weightsOnEdge

      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, nz1, itr, itrp, cell1, cell2, nCellsSolve
      integer :: index_qv

      real (kind=RKIND) :: ptop, pitop, ptopb, p0, flux, d2fdx2_cell1, d2fdx2_cell2

      real (kind=RKIND) :: ztemp, zd, zt, dz, str

      real (kind=RKIND), dimension(grid % nVertLevels, grid % nCells) :: rh
      real (kind=RKIND) :: ptmp, es, qvs, xnutr, ptemp
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
      zf => grid % zf % array
      zf3 => grid % zf3 % array
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
         sm = .25*min((zc(k)-zc(k-1))/dz,1.)
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
99       print *,"PASS-SHP"
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
              scalars(index_qv,k,i) = amin1(0.014,rh(k,i)*qvs)
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
         write(0,'(8(f14.9,2x))') .5*(zgrid(k,1)+zgrid(k+1,1))/1000.,   &
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
  
                  if (k /= 1) then
                     zf(k,1,iEdge) = ( fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1) )*zb(k,1,iEdge)
                     zf(k,2,iEdge) = ( fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2) )*zb(k,2,iEdge)
                     zf3(k,1,iEdge)= ( fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1) )*zb3(k,1,iEdge)
                     zf3(k,2,iEdge)= ( fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2) )*zb3(k,2,iEdge)
                  end if

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
            diag % rw % array(k,cell2) = diag % rw % array(k,cell2) + zf(k,2,iEdge)*flux 
            diag % rw % array(k,cell1) = diag % rw % array(k,cell1) - zf(k,1,iEdge)*flux 

            if (config_theta_adv_order ==3) then
               diag % rw % array(k,cell2) = diag % rw % array(k,cell2)    &
                                            - sign(1.,ru(k,iEdge))*config_coef_3rd_order*zf3(k,2,iEdge)*flux
               diag % rw % array(k,cell1) = diag % rw % array(k,cell1)    &
                                            + sign(1.,ru(k,iEdge))*config_coef_3rd_order*zf3(k,1,iEdge)*flux
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
            diag % theta % array(k,iCell) = state % theta_m % array(k,iCell) / (1.0 + 1.61 * scalars(state % index_qv,k,iCell))
         end do
      end do

   end subroutine init_atm_test_case_mtn_wave


   subroutine init_atm_test_case_gfs(dminfo, grid, fg, state, diag, test_case, parinfo)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Real-data test case using GFS data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use mpas_dmpar
      use init_atm_read_met
      use init_atm_llxy
      use init_atm_hinterp

      implicit none

      type (dm_info), intent(in) :: dminfo
      type (mesh_type), intent(inout) :: grid
      type (fg_type), intent(inout) :: fg
      type (state_type), intent(inout) :: state
      type (diag_type), intent(inout) :: diag
      integer, intent(in) :: test_case
      type (parallel_info), pointer :: parinfo

      real (kind=RKIND), parameter :: u0 = 35.0
      real (kind=RKIND), parameter :: alpha_grid = 0.  ! no grid rotation
      real (kind=RKIND), parameter :: omega_e = 7.29212e-05
      real (kind=RKIND), parameter :: t0b = 250., t0 = 288., delta_t = 4.8e+05, dtdz = 0.005, eta_t = 0.2
      real (kind=RKIND), parameter :: u_perturbation = 1., pert_radius = 0.1, latitude_pert = 40., longitude_pert = 20.
      real (kind=RKIND), parameter :: theta_c = pii/4.0
      real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
      real (kind=RKIND), parameter :: rh_max = 0.4       ! Maximum relative humidity
      real (kind=RKIND), parameter :: k_x = 9.           ! Normal mode wave number

      type (met_data) :: field
      type (proj_info) :: proj

      real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp
      real (kind=RKIND), dimension(:), pointer :: vert_level, latPoints, lonPoints, ter
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zx, zz, hx
      real (kind=RKIND), dimension(:,:), pointer :: pressure, ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, dss, t, rt
      real (kind=RKIND), dimension(:), pointer :: destField1d
      real (kind=RKIND), dimension(:,:), pointer :: destField2d
      real (kind=RKIND), dimension(:,:,:), pointer :: zf, zf3, zb, zb3
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two

      real (kind=RKIND) :: target_z
      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, nz1, itr, itrp, cell1, cell2, nCellsSolve
      integer :: nInterpPoints, ndims

      integer, dimension(5) :: interp_list
      real (kind=RKIND) :: maskval
      real (kind=RKIND) :: msgval
      real (kind=RKIND) :: fillval
      integer :: masked

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnEdge, nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnEdge, cellsOnEdge, edgesOnCell, cellsOnCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, AreaCell 
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge
      real (kind=RKIND), dimension(:,:), pointer :: v
      real (kind=RKIND), dimension(:,:), pointer :: sorted_arr

      real(kind=RKIND), dimension(:), pointer :: hs, hs1
      real(kind=RKIND) :: hm, zh, dzmin, dzmina, dzmina_global, dzminf, sm
      integer :: nsmterrain, kz, sfc_k
      logical :: hybrid, smooth

      integer :: it
      real (kind=RKIND) :: p_check

      ! For interpolating terrain and land use
      integer :: nx, ny, nzz, iPoint, subx, suby
      integer :: isigned, endian, wordsize, istatus
      integer :: iTileStart, iTileEnd
      integer :: jTileStart, jTileEnd
      integer, allocatable, dimension(:) :: nhs
      integer, allocatable, dimension(:,:) :: ncat
      real (kind=4) :: scalefactor   ! NB: this should be a single-precision real
      real (kind=RKIND) :: lat_pt, lon_pt, lon_pt_o
      real (kind=4), allocatable, dimension(:,:,:) :: rarray    ! NB: this should be a single-precision real array
      real (kind=RKIND), allocatable, dimension(:,:) :: rslab, maskslab
      real (kind=RKIND), allocatable, dimension(:,:) :: maxsnowalb
      real (kind=RKIND), allocatable, dimension(:,:) :: soiltemp_1deg
      real (kind=RKIND), allocatable, dimension(:,:,:) :: vegfra
      integer, dimension(:), pointer :: mask_array
      integer, dimension(grid % nEdges), target :: edge_mask
      character (len=1024) :: fname

      real (kind=RKIND) :: u, flux, fluxk, lat1, lat2, eta_v, r_pert, u_pert, lat_pert, lon_pert, r
      real (kind=RKIND) :: lat, lon, x, y

      real (kind=RKIND) :: ptop, p0, phi
      real (kind=RKIND) :: lon_Edge

      real (kind=RKIND) :: r_earth, etavs, ztemp, zd, zt, dz, gam, delt, str

      real (kind=RKIND), dimension(grid % nVertLevels, grid % nCells) :: rel_hum, temperature, qv
      real (kind=RKIND) :: ptmp, es, rs, rgas_moist, qvs, xnutr, znut, ptemp, rcv
      integer :: iter

      real (kind=RKIND), dimension(grid % nVertLevels + 1) :: hyai, hybi, znu, znw, znwc, znwv, hyam, hybm
      real (kind=RKIND), dimension(grid % nVertLevels + 1) :: znuc, znuv, bn, divh, dpn

      real (kind=RKIND), dimension(grid % nVertLevels + 1) :: sh, zw, ah
      real (kind=RKIND), dimension(grid % nVertLevels) :: zu, dzw, rdzwp, rdzwm
      real (kind=RKIND), dimension(grid % nVertLevels) :: eta, etav, teta, ppi, tt

      real (kind=RKIND) :: d1, d2, d3, cf1, cf2, cf3, cof1, cof2, psurf

      !  storage for (lat,z) arrays for zonal velocity calculation

      integer, parameter :: nlat=361
      real (kind=RKIND), dimension(grid % nVertLevels + 1) :: zz_1d, zgrid_1d, hx_1d
      real (kind=RKIND), dimension(grid % nVertLevels) :: flux_zonal
      real (kind=RKIND), dimension(nlat, grid % nVertLevels) :: u_2d, etavs_2d
      real (kind=RKIND), dimension(grid % nVertLevels + 1) :: fsum
      real (kind=RKIND), dimension(nlat) :: lat_2d
      real (kind=RKIND) :: dlat
      real (kind=RKIND) :: z_edge, z_edge3, d2fdx2_cell1, d2fdx2_cell2

      weightsOnEdge     => grid % weightsOnEdge % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      nEdgesOnCell      => grid % nEdgesOnCell % array
      edgesOnEdge       => grid % edgesOnEdge % array
      edgesOnCell       => grid % edgesOnCell % array
      dvEdge            => grid % dvEdge % array
      dcEdge            => grid % dcEdge % array
      AreaCell          => grid % AreaCell % array
      CellsOnEdge       => grid % CellsOnEdge % array
      cellsOnCell       => grid % cellsOnCell % array

      deriv_two  => grid % deriv_two % array
      zf  => grid % zf % array
      zf3 => grid % zf3% array
      zb  => grid % zb % array
      zb3 => grid % zb3% array

      zgrid => grid % zgrid % array
      rdzw => grid % rdzw % array
      dzu => grid % dzu % array
      rdzu => grid % rdzu % array
      fzm => grid % fzm % array
      fzp => grid % fzp % array
      zx => grid % zx % array
      zz => grid % zz % array
      hx => grid % hx % array
      ter => grid % ter % array
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

      scalars => state % scalars % array
      
      nz1 = grid % nVertLevels
      nz = nz1 + 1
      nCellsSolve = grid % nCellsSolve

      xnutr = 0.
      zd = 12000.
      znut = eta_t

      etavs = (1.-0.252)*pii/2.
      rcv = rgas/(cp-rgas)
      r_earth = a
      p0 = 1.e+05

      interp_list(1) = FOUR_POINT
      interp_list(2) = SEARCH
      interp_list(3) = 0


      !
      ! Scale all distances and areas from a unit sphere to one with radius a
      !

      if (config_static_interp) then

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

      scalars(:,:,:) = 0.

      call atm_initialize_advection_rk(grid) 
      call atm_initialize_deformation_weights(grid) 


      !
      ! Interpolate HGT
      !
!     nx = 126
!     ny = 126
      nx = 1206
      ny = 1206
      nzz = 1
      isigned = 1
      endian = 0
      wordsize = 2
      scalefactor = 1.0
      allocate(rarray(nx,ny,nzz))
      allocate(nhs(grid % nCells))
      nhs(:) = 0
      grid % ter % array(:) = 0.0

      do jTileStart=1,20401,ny-6
!     do jTileStart=1,961,ny-6
         jTileEnd = jTileStart + ny - 1 - 6
         do iTileStart=1,42001,nx-6
!        do iTileStart=1,2041,nx-6
            iTileEnd = iTileStart + nx - 1 - 6
            write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)//'topo_30s/',iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
!           write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)//'topo_10m/',iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
write(0,*) trim(fname)

            call read_geogrid(fname, len_trim(fname), &
                              rarray, &
                              nx, ny, nzz, &
                              isigned, endian, scalefactor, wordsize, istatus)
            call init_atm_check_read_error(istatus, fname, dminfo)

            iPoint = 1
            do j=4,ny-3
            do i=4,nx-3
               lat_pt = -89.99583 + (jTileStart + j - 5) * 0.0083333333
               lon_pt = -179.99583 + (iTileStart + i - 5) * 0.0083333333
!              lat_pt = -89.91667 + (jTileStart + j - 5) * 0.166667
!              lon_pt = -179.91667 + (iTileStart + i - 5) * 0.166667
               lat_pt = lat_pt * pii / 180.0
               lon_pt = lon_pt * pii / 180.0

               iPoint = nearest_cell(lat_pt, lon_pt, &
                                     iPoint, &
                                     grid % nCells, grid % maxEdges, grid % nEdgesOnCell % array, grid % cellsOnCell % array, &
                                     grid % latCell % array, grid % lonCell % array)

               grid % ter % array(iPoint) = grid % ter % array(iPoint) + rarray(i,j,1)
               nhs(iPoint) = nhs(iPoint) + 1

            end do
            end do

          end do
       end do

      do iCell=1, grid % nCells
         grid % ter % array(iCell) = grid % ter % array(iCell) / real(nhs(iCell))
      end do

      deallocate(rarray)
      deallocate(nhs)


      !
      ! Interpolate LU_INDEX
      !
      nx = 1200
      ny = 1200
      nzz = 1
      isigned = 1
      endian = 0
      wordsize = 1
      scalefactor = 1.0
      allocate(rarray(nx,ny,nzz))
      allocate(ncat(24,grid % nCells))
      ncat(:,:) = 0
      grid % lu_index % array(:) = 0.0

      do jTileStart=1,20401,ny
         jTileEnd = jTileStart + ny - 1
         do iTileStart=1,42001,nx
            iTileEnd = iTileStart + nx - 1
            write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)//'/landuse_30s/',iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
write(0,*) trim(fname)

            call read_geogrid(fname, len_trim(fname), &
                              rarray, &
                              nx, ny, nzz, &
                              isigned, endian, scalefactor, wordsize, istatus)
            call init_atm_check_read_error(istatus, fname, dminfo)

            iPoint = 1
            do j=1,ny
            do i=1,nx
               lat_pt = -89.99583 + (jTileStart + j - 2) * 0.0083333333
               lon_pt = -179.99583 + (iTileStart + i - 2) * 0.0083333333
               lat_pt = lat_pt * pii / 180.0
               lon_pt = lon_pt * pii / 180.0

               iPoint = nearest_cell(lat_pt, lon_pt, &
                                     iPoint, &
                                     grid % nCells, grid % maxEdges, grid % nEdgesOnCell % array, grid % cellsOnCell % array, &
                                     grid % latCell % array, grid % lonCell % array)

               ncat(int(rarray(i,j,1)),iPoint) = ncat(int(rarray(i,j,1)),iPoint) + 1

            end do
            end do

         end do
      end do

      do iCell=1, grid % nCells
         grid % lu_index % array(iCell) = 1
         do i=2,24
            if (ncat(i,iCell) > ncat(grid % lu_index % array(iCell),iCell)) then
               grid % lu_index % array(iCell) = i
            end if
         end do
      end do

      deallocate(rarray)
      deallocate(ncat)


      !
      ! Interpolate SOILCAT_TOP
      !
      nx = 1200
      ny = 1200
      nzz = 1
      isigned = 1
      endian = 0
      wordsize = 1
      scalefactor = 1.0
      allocate(rarray(nx,ny,nzz))
      allocate(ncat(16,grid % nCells))
      ncat(:,:) = 0
      grid % soilcat_top % array(:) = 0.0

      do jTileStart=1,20401,ny
         jTileEnd = jTileStart + ny - 1
         do iTileStart=1,42001,nx
            iTileEnd = iTileStart + nx - 1
            write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)//'/soiltype_top_30s/',iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
write(0,*) trim(fname)

            call read_geogrid(fname, len_trim(fname), &
                              rarray, &
                              nx, ny, nzz, &
                              isigned, endian, scalefactor, wordsize, istatus)
            call init_atm_check_read_error(istatus, fname, dminfo)

            iPoint = 1
            do j=1,ny
            do i=1,nx
               lat_pt = -89.99583 + (jTileStart + j - 2) * 0.0083333333
               lon_pt = -179.99583 + (iTileStart + i - 2) * 0.0083333333
               lat_pt = lat_pt * pii / 180.0
               lon_pt = lon_pt * pii / 180.0

               iPoint = nearest_cell(lat_pt, lon_pt, &
                                     iPoint, &
                                     grid % nCells, grid % maxEdges, grid % nEdgesOnCell % array, grid % cellsOnCell % array, &
                                     grid % latCell % array, grid % lonCell % array)

               ncat(int(rarray(i,j,1)),iPoint) = ncat(int(rarray(i,j,1)),iPoint) + 1

            end do
            end do

         end do
      end do

      do iCell=1, grid % nCells
         grid % soilcat_top % array(iCell) = 1
         do i=2,16
            if (ncat(i,iCell) > ncat(grid % soilcat_top % array(iCell),iCell)) then
               grid % soilcat_top % array(iCell) = i
            end if
         end do
      end do

      deallocate(rarray)
      deallocate(ncat)


      !
      ! Interpolate SOILCAT_BOT
      !
      nx = 1200
      ny = 1200
      nzz = 1
      isigned = 1
      endian = 0
      wordsize = 1
      scalefactor = 1.0
      allocate(rarray(nx,ny,nzz))
      allocate(ncat(16,grid % nCells))
      ncat(:,:) = 0
      grid % soilcat_bot % array(:) = 0.0

      do jTileStart=1,20401,ny
         jTileEnd = jTileStart + ny - 1
         do iTileStart=1,42001,nx
            iTileEnd = iTileStart + nx - 1
            write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)//'/soiltype_bot_30s/',iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
write(0,*) trim(fname)

            call read_geogrid(fname, len_trim(fname), &
                              rarray, &
                              nx, ny, nzz, &
                              isigned, endian, scalefactor, wordsize, istatus)
            call init_atm_check_read_error(istatus, fname, dminfo)

            iPoint = 1
            do j=1,ny
            do i=1,nx
               lat_pt = -89.99583 + (jTileStart + j - 2) * 0.0083333333
               lon_pt = -179.99583 + (iTileStart + i - 2) * 0.0083333333
               lat_pt = lat_pt * pii / 180.0
               lon_pt = lon_pt * pii / 180.0

               iPoint = nearest_cell(lat_pt, lon_pt, &
                                     iPoint, &
                                     grid % nCells, grid % maxEdges, grid % nEdgesOnCell % array, grid % cellsOnCell % array, &
                                     grid % latCell % array, grid % lonCell % array)

               ncat(int(rarray(i,j,1)),iPoint) = ncat(int(rarray(i,j,1)),iPoint) + 1

            end do
            end do

         end do
      end do

      do iCell=1, grid % nCells
         grid % soilcat_bot % array(iCell) = 1
         do i=2,16
            if (ncat(i,iCell) > ncat(grid % soilcat_bot % array(iCell),iCell)) then
               grid % soilcat_bot % array(iCell) = i
            end if
         end do
      end do

      deallocate(rarray)
      deallocate(ncat)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! KLUDGE TO FIX SOIL TYPE OVER ANTARCTICA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      where (grid % lu_index % array == 24) grid % soilcat_top % array = 16
      where (grid % lu_index % array == 24) grid % soilcat_bot % array = 16

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CORRECT INCONSISTENT SOIL AND LAND USE DATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do iCell = 1,grid % nCells
         if (grid % lu_index % array(iCell) == 16 .or. &
             grid % soilcat_top % array(iCell) == 14 .or. &
             grid % soilcat_bot % array(iCell) == 14) then
            if (grid % lu_index % array(iCell) /= 16) then
               write(0,*) 'Turning lu_index into water at ', iCell
               grid % lu_index % array(iCell) = 16
            end if
            if (grid % soilcat_top % array(iCell) /= 14) then
               write(0,*) 'Turning soilcat_top into water at ', iCell
               grid % soilcat_top % array(iCell) = 14
            end if
            if (grid % soilcat_bot % array(iCell) /= 14) then
               write(0,*) 'Turning soilcat_bot into water at ', iCell
               grid % soilcat_bot % array(iCell) = 14
            end if
         end if
      end do


      !
      ! Derive LANDMASK
      !
      grid % landmask % array(:) = 0
      do iCell=1, grid % nCells
         if (grid % lu_index % array(iCell) /= 16) grid % landmask % array(iCell) = 1
      end do


      !
      ! Interpolate SOILTEMP:
      !
      nx = 186
      ny = 186
      nzz = 1
      isigned = 0
      endian = 0
      wordsize = 2
      scalefactor = 0.01
      allocate(rarray(nx,ny,nzz))
      allocate(soiltemp_1deg(360,180))
      grid % soiltemp % array(:) = 0.0

      call map_set(PROJ_LATLON, proj, &
                   latinc = 1.0, &
                   loninc = 1.0, &
                   knowni = 1.0, &
                   knownj = 1.0, &
                   lat1 = -89.5, &
                   lon1 = -179.5)

      write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)//'soiltemp_1deg/',1,'-',180,'.',1,'-',180
write(0,*) trim(fname)
      call read_geogrid(fname, len_trim(fname), &
                        rarray, &
                        nx, ny, nzz, &
                        isigned, endian, scalefactor, wordsize, istatus)
      call init_atm_check_read_error(istatus, fname, dminfo)

      soiltemp_1deg(1:180,1:180) = rarray(4:183,4:183,1)

      write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)//'soiltemp_1deg/',181,'-',360,'.',1,'-',180
write(0,*) trim(fname)
      call read_geogrid(fname, len_trim(fname), &
                        rarray, &
                        nx, ny, nzz, &
                        isigned, endian, scalefactor, wordsize, istatus)
      call init_atm_check_read_error(istatus, fname, dminfo)

      soiltemp_1deg(181:360,1:180) = rarray(4:183,4:183,1)

      interp_list(1) = FOUR_POINT
      interp_list(2) = W_AVERAGE4
      interp_list(3) = W_AVERAGE16
      interp_list(4) = SEARCH
      interp_list(5) = 0

      do iCell=1,grid%nCells
  
         if (grid % landmask % array(iCell) == 1) then
            lat = grid % latCell % array(iCell)*DEG_PER_RAD
            lon = grid % lonCell % array(iCell)*DEG_PER_RAD
            call latlon_to_ij(proj, lat, lon, x, y)
            if (x < 0.5) then
               lon = lon + 360.0
               call latlon_to_ij(proj, lat, lon, x, y)
            else if (x >= 360.5) then
               lon = lon - 360.0
               call latlon_to_ij(proj, lat, lon, x, y)
            end if
if (y < 1.0) y = 1.0
if (y > 179.0) y = 179.0
!            grid % soiltemp % array(iCell) = interp_sequence(x, y, 1, soiltemp_1deg, 1, 360, 1, 180, 1, 1, -1.e30, interp_list, 1)
            grid % soiltemp % array(iCell) = interp_sequence(x, y, 1, soiltemp_1deg, 1, 360, 1, 180, 1, 1, 0., interp_list, 1)
         else
            grid % soiltemp % array(iCell) = 0.0
         end if

      end do

      deallocate(rarray)
      deallocate(soiltemp_1deg)


      !
      ! Interpolate SNOALB
      !
      nx = 186
      ny = 186
      nzz = 1
      isigned = 0
      endian = 0
      wordsize = 1
      scalefactor = 1.0
      allocate(rarray(nx,ny,nzz))
      allocate(maxsnowalb(360,180))
      grid % snoalb % array(:) = 0.0

      call map_set(PROJ_LATLON, proj, &
                   latinc = 1.0, &
                   loninc = 1.0, &
                   knowni = 1.0, &
                   knownj = 1.0, &
                   lat1 = -89.5, &
                   lon1 = -179.5)

      write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)//'maxsnowalb/',1,'-',180,'.',1,'-',180
write(0,*) trim(fname)
      call read_geogrid(fname, len_trim(fname), &
                        rarray, &
                        nx, ny, nzz, &
                        isigned, endian, scalefactor, wordsize, istatus)
      call init_atm_check_read_error(istatus, fname, dminfo)

      maxsnowalb(1:180,1:180) = rarray(4:183,4:183,1)

      write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)//'maxsnowalb/',181,'-',360,'.',1,'-',180
write(0,*) trim(fname)
      call read_geogrid(fname, len_trim(fname), &
                        rarray, &
                        nx, ny, nzz, &
                        isigned, endian, scalefactor, wordsize, istatus)
      call init_atm_check_read_error(istatus, fname, dminfo)

      maxsnowalb(181:360,1:180) = rarray(4:183,4:183,1)

      interp_list(1) = FOUR_POINT
      interp_list(2) = W_AVERAGE4
      interp_list(3) = W_AVERAGE16
      interp_list(4) = SEARCH
      interp_list(5) = 0

      do iCell=1,grid%nCells
  
         if (grid % landmask % array(iCell) == 1) then
            lat = grid % latCell % array(iCell)*DEG_PER_RAD
            lon = grid % lonCell % array(iCell)*DEG_PER_RAD
            call latlon_to_ij(proj, lat, lon, x, y)
            if (x < 0.5) then
               lon = lon + 360.0
               call latlon_to_ij(proj, lat, lon, x, y)
            else if (x >= 360.5) then
               lon = lon - 360.0
               call latlon_to_ij(proj, lat, lon, x, y)
            end if
if (y < 1.0) y = 1.0
if (y > 179.0) y = 179.0
!            grid % snoalb % array(iCell) = interp_sequence(x, y, 1, maxsnowalb, 1, 360, 1, 180, 1, 1, -1.e30, interp_list, 1)
            grid % snoalb % array(iCell) = interp_sequence(x, y, 1, maxsnowalb, 1, 360, 1, 180, 1, 1, 0., interp_list, 1)
         else
            grid % snoalb % array(iCell) = 0.0
         end if

      end do

      grid % snoalb % array(:) = grid % snoalb % array(:) / 100.0

      deallocate(rarray)
      deallocate(maxsnowalb)


      !
      ! Interpolate GREENFRAC
      !
      nx = 1256
      ny = 1256
      nzz = 12
      isigned = 0
      endian = 0
      wordsize = 1
      scalefactor = 1.0
      allocate(rarray(nx,ny,nzz))
      allocate(vegfra(2500,1250,12))
!     grid % vegfra % array(:) = 0.0
      grid % greenfrac % array(:,:) = 0.0

      call map_set(PROJ_LATLON, proj, &
                   latinc = 0.144, &
                   loninc = 0.144, &
                   knowni = 1.0, &
                   knownj = 1.0, &
                   lat1 = -89.928, &
                   lon1 = -179.928)

      write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)//'greenfrac/',1,'-',1250,'.',1,'-',1250
write(0,*) trim(fname)
      call read_geogrid(fname, len_trim(fname), &
                        rarray, &
                        nx, ny, nzz, &
                        isigned, endian, scalefactor, wordsize, istatus)
      call init_atm_check_read_error(istatus, fname, dminfo)

      vegfra(1:1250,1:1250,1:12) = rarray(4:1253,4:1253,1:12)

      write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)//'greenfrac/',1251,'-',2500,'.',1,'-',1250
write(0,*) trim(fname)
      call read_geogrid(fname, len_trim(fname), &
                        rarray, &
                        nx, ny, nzz, &
                        isigned, endian, scalefactor, wordsize, istatus)
      call init_atm_check_read_error(istatus, fname, dminfo)

      vegfra(1251:2500,1:1250,1:12) = rarray(4:1253,4:1253,1:12)

      do iCell=1,grid%nCells
         if (grid % landmask % array(iCell) == 1) then
            lat = grid % latCell % array(iCell)*DEG_PER_RAD
            lon = grid % lonCell % array(iCell)*DEG_PER_RAD
            call latlon_to_ij(proj, lat, lon, x, y)
            if (x < 0.5) then
               lon = lon + 360.0
               call latlon_to_ij(proj, lat, lon, x, y)
            else if (x >= 2500.5) then
               lon = lon - 360.0
               call latlon_to_ij(proj, lat, lon, x, y)
            end if
if (y < 1.0) y = 1.0
if (y > 1249.0) y = 1249.0
            do k=1,12
               grid % greenfrac % array(k,iCell) = interp_sequence(x, y, k, vegfra, 1, 2500, 1, 1250, 1, 12, -1.e30, interp_list, 1)
            end do
         else
            grid % greenfrac % array(:,iCell) = 0.0
         end if
         grid % shdmin % array(iCell) = minval(grid % greenfrac % array(:,iCell))
         grid % shdmax % array(iCell) = maxval(grid % greenfrac % array(:,iCell))
      
      end do

      deallocate(rarray)
      deallocate(vegfra)


      !
      ! Interpolate ALBEDO12M
      !
      nx = 1256
      ny = 1256
      nzz = 12
      isigned = 0
      endian = 0
      wordsize = 1
      scalefactor = 1.0
      allocate(rarray(nx,ny,nzz))
      allocate(vegfra(2500,1250,12))
      grid % albedo12m % array(:,:) = 0.0

      call map_set(PROJ_LATLON, proj, &
                   latinc = 0.144, &
                   loninc = 0.144, &
                   knowni = 1.0, &
                   knownj = 1.0, &
                   lat1 = -89.928, &
                   lon1 = -179.928)

      write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)//'albedo_ncep/',1,'-',1250,'.',1,'-',1250
write(0,*) trim(fname)
      call read_geogrid(fname, len_trim(fname), &
                        rarray, &
                        nx, ny, nzz, &
                        isigned, endian, scalefactor, wordsize, istatus)
      call init_atm_check_read_error(istatus, fname, dminfo)

      vegfra(1:1250,1:1250,1:12) = rarray(4:1253,4:1253,1:12)

      write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)//'albedo_ncep/',1251,'-',2500,'.',1,'-',1250
write(0,*) trim(fname)
      call read_geogrid(fname, len_trim(fname), &
                        rarray, &
                        nx, ny, nzz, &
                        isigned, endian, scalefactor, wordsize, istatus)
      call init_atm_check_read_error(istatus, fname, dminfo)

      vegfra(1251:2500,1:1250,1:12) = rarray(4:1253,4:1253,1:12)

      do iCell=1,grid%nCells
         if (grid % landmask % array(iCell) == 1) then
            lat = grid % latCell % array(iCell)*DEG_PER_RAD
            lon = grid % lonCell % array(iCell)*DEG_PER_RAD
            call latlon_to_ij(proj, lat, lon, x, y)
            if (x < 0.5) then
               lon = lon + 360.0
               call latlon_to_ij(proj, lat, lon, x, y)
            else if (x >= 2500.5) then
               lon = lon - 360.0
               call latlon_to_ij(proj, lat, lon, x, y)
            end if
if (y < 1.0) y = 1.0
if (y > 1249.0) y = 1249.0
            do k=1,12
               grid % albedo12m % array(k,iCell) = interp_sequence(x, y, k, vegfra, 1, 2500, 1, 1250, 1, 12, 0.0, interp_list, 1)
            end do
         else
            grid % albedo12m % array(:,iCell) = 8.0
         end if
      end do

      deallocate(rarray)
      deallocate(vegfra)

 
      end if    ! config_static_interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BEGIN ADOPT GFS TERRAIN HEIGHT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if 0
      call read_met_init(trim(config_met_prefix), .false., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         write(0,*) 'Error reading initial met data'
         return
      end if

      call read_next_met_field(field, istatus)
      do while (istatus == 0)
         if (index(field % field, 'SOILHGT') /= 0) then

            !
            ! Set up projection
            !
            call map_init(proj)
          
            if (field % iproj == PROJ_LATLON) then
               call map_set(PROJ_LATLON, proj, &
                            latinc = real(field % deltalat), &
                            loninc = real(field % deltalon), &
                            knowni = 1.0, &
                            knownj = 1.0, &
                            lat1 = real(field % startlat), &
                            lon1 = real(field % startlon))
            end if


            if (index(field % field, 'SOILHGT') /= 0) then
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField1d => grid % ter % array
               ndims = 1
            end if

            do i=1,nInterpPoints
               lat = latPoints(i)*DEG_PER_RAD
               lon = lonPoints(i)*DEG_PER_RAD
               call latlon_to_ij(proj, lat, lon, x, y)
               if (x < 0.5) then
                  lon = lon + 360.0
                  call latlon_to_ij(proj, lat, lon, x, y)
               end if
               if (ndims == 1) then
                  destField1d(i) = interp_sequence(x, y, 1, field % slab, 1, field % nx, 1, field % ny, 1, 1, -1.e30, interp_list, 1)
               else if (ndims == 2) then
                  destField2d(k,i) = interp_sequence(x, y, 1, field % slab, 1, field % nx, 1, field % ny, 1, 1, -1.e30, interp_list, 1)
               end if
            end do
         end if
   
         deallocate(field % slab)
         call read_next_met_field(field, istatus)
      end do

      call read_met_close()
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! END ADOPT GFS TERRAIN HEIGHT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      if (config_vertical_grid) then

      !
      ! Vertical grid setup
      !
      allocate(hs (grid%nCells+1))
      allocate(hs1(grid%nCells+1))

!     Fourth order smoother for terrain

      nsmterrain = config_nsmterrain

      do i=1,nsmterrain

         do iCell=1,grid%nCells
            hs(iCell) = 0.
            do j = 1,nEdgesOnCell(iCell)
               hs(iCell) = hs(iCell) + dvEdge(edgesOnCell(j,iCell))    &
                                     / dcEdge(edgesOnCell(j,iCell))    &
                                     *   (ter(cellsOnCell(j,iCell))-ter(iCell))
            end do
            hs(iCell) = ter(iCell) + 0.25*hs(iCell)
         end do

         do iCell=1,grid %nCells
            ter(iCell) = 0.
            do j = 1,nEdgesOnCell(iCell)
               ter(iCell) = ter(iCell) + dvEdge(edgesOnCell(j,iCell))    &
                                     / dcEdge(edgesOnCell(j,iCell))    &
                                     *    (hs(cellsOnCell(j,iCell))-hs(iCell))
            end do
            ter(iCell) = hs(iCell) - 0.25*ter(iCell)
         end do

         call mpas_dmpar_exch_halo_field1d_real(dminfo, ter(:), &
                                          grid % nCells, &
                                          parinfo % cellsToSend, parinfo % cellsToRecv)


      end do

      do iCell=1,grid % nCells
         hx(:,iCell) = ter(iCell)
      end do

      hm = maxval(ter(:))
      write(0,*) "max ter = ", hm

!     Metrics for hybrid coordinate and vertical stretching

      str = 1.5
!      str = 1.
      zt = config_ztop
      dz = zt/float(nz1)

      do k=1,nz
         zw(k) = (real(k-1)/real(nz1))**str*zt
         if (k > 1) dzw(k-1) = zw(k)-zw(k-1)
      end do

!     ah(k) governs the transition between terrain-following 
!        and pure height coordinates
!           ah(k) = 1           is a smoothed terrain-following coordinate
!           ah(k) = 1.-zw(k)/zt is the basic terrain-following coordinate
!           ah(k) = 0           is a height coordinate
 
      hybrid = .true.
!      hybrid = .false.

      kz = nz
      if (hybrid) then
      
         zh = zt
!         zh = 0.5*zt

         do k=1,nz
            if (zw(k) < zh) then
               ah(k) = cos(.5*pii*zw(k)/zh)**6

!!!               ah(k) = ah(k)*(1.-zw(k)/zt)

            else
               ah(k) = 0.
               kz = min(kz,k)
            end if
         end do

      else
	
         do k=1,nz
            ah(k) = 1.-zw(k)/zt
         end do

      end if

      do k=1,nz
         write(0,*) k,zw(k), ah(k)
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

!     Smoothing algorithm for coordinate surfaces 

      smooth = config_smooth_surfaces
!      smooth = .false.

      if (smooth) then

         dzmin = 0.3

         do k=2,kz-1
            hx(k,:) = hx(k-1,:)
            dzminf = zw(k)-zw(k-1)

!            dzmin = max(.5,1.-.5*zw(k)/hm)

            sm = .05*min(.5*zw(k)/hm,1.)
          
            do i=1,50
               do iCell=1,grid %nCells
                  hs1(iCell) = 0.
                  do j = 1,nEdgesOnCell(iCell)

                     hs1(iCell) = hs1(iCell) + dvEdge(edgesOnCell(j,iCell))    &
                                           / dcEdge(edgesOnCell(j,iCell))    &
                                           *  (hx(k,cellsOnCell(j,iCell))-hx(k,iCell))
                  end do
                  hs1(iCell) = hx(k,iCell) + sm*hs1(iCell)

                  hs(iCell) = 0.
              !    do j = 1,nEdgesOnCell(iCell)
              !       hs(iCell) = hs(iCell) + dvEdge(edgesOnCell(j,iCell))    &
              !                             / dcEdge(edgesOnCell(j,iCell))    &
              !                             *  (hs1(cellsOnCell(j,iCell))-hs1(iCell))
              !    end do
                  hs(iCell) = hs1(iCell) - 0.*hs(iCell)

               end do

               call mpas_dmpar_exch_halo_field1d_real(dminfo, hs(:), &
                                                grid % nCells, &
                                                parinfo % cellsToSend, parinfo % cellsToRecv)

             !  dzmina = minval(hs(:)-hx(k-1,:))
               dzmina = minval(zw(k)+ah(k)*hs(:)-zw(k-1)-ah(k-1)*hx(k-1,:))
               call mpas_dmpar_min_real(dminfo, dzmina, dzmina_global)
             !  write(0,*) ' k,i, dzmina, dzmin, zw(k)-zw(k-1) ', k,i, dzmina, dzmin, zw(k)-zw(k-1)
               if (dzmina_global >= dzmin*(zw(k)-zw(k-1))) then
                  hx(k,:)=hs(:)
                  dzminf = dzmina_global
               else
                  exit
               end if
            end do
            write(0,*) k,i,sm,dzminf/(zw(k)-zw(k-1)),dzmina/(zw(k)-zw(k-1))
         end do

         do k=kz,nz
               hx(k,:) = 0.
         end do
      else

         do k=2,nz1
            dzmina = minval(zw(k)+ah(k)*hx(k,:)-zw(k-1)-ah(k-1)*hx(k-1,:))
            write(0,*) k,dzmina/(zw(k)-zw(k-1))
         end do

      end if

      deallocate(hs )
      deallocate(hs1)

!     Height of coordinate levels (calculation of zgrid)

      do iCell=1,grid % nCells
         do k=1,nz	
            zgrid(k,iCell) = zw(k) + ah(k)*hx(k,iCell)
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
           if (ztemp.gt.zd+.1)  then
               dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
           end if
         end do
      enddo

!      do k=1,nz1
!         write(0,*) ' k, zgrid(k,1),hx(k,1) ',k,zgrid(k,1),hx(k,1)
!      enddo

!      do k=1,nz1
!         write(0,*) ' k, zx(k,1) ',k,zx(k,1)
!      enddo


      ! For z-metric term in omega equation
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
  
                  if (k /= 1) then
                     zf(k,1,iEdge) = ( fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1) )*zb(k,1,iEdge)
                     zf(k,2,iEdge) = ( fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2) )*zb(k,2,iEdge)
                     zf3(k,1,iEdge)= ( fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1) )*zb3(k,1,iEdge)
                     zf3(k,2,iEdge)= ( fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2) )*zb3(k,2,iEdge)
                  end if

            end do

         end if
      end do

      write(0,*) ' grid metrics setup complete '

      end if    ! config_vertical_grid


      if (config_met_interp) then

      !
      ! First, try to locate the LANDSEA field for use as an interpolation mask
      !
      call read_met_init(trim(config_met_prefix), .false., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         write(0,*) 'Error reading initial met data'
         return
      end if

      call read_next_met_field(field, istatus)

      do while (istatus == 0)
         if (index(field % field, 'LANDSEA') /= 0) then

            allocate(maskslab(-3:field % nx+3, field % ny))
            maskslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
            maskslab(0, 1:field % ny)  = field % slab(field % nx, 1:field % ny)
            maskslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
            maskslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
            maskslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
            maskslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
            maskslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)
write(0,*) 'minval, maxval of LANDSEA = ', minval(maskslab), maxval(maskslab)

         end if
   
         deallocate(field % slab)
         call read_next_met_field(field, istatus)
      end do

      call read_met_close()

      edge_mask(:) = 1


      !
      ! Horizontally interpolate meteorological data
      !
      allocate(vert_level(config_nfglevels))
      vert_level(:) = -1.0

      call read_met_init(trim(config_met_prefix), .false., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         write(0,*) 'Error reading initial met data'
         return
      end if

      call read_next_met_field(field, istatus)

      do while (istatus == 0)

         interp_list(1) = FOUR_POINT
         interp_list(2) = SEARCH
         interp_list(3) = 0

         maskval = -1.0
         masked = -1
         fillval = 0.0
         msgval = -1.e30

         mask_array => grid % landmask % array

         if (index(field % field, 'UU') /= 0 .or. &
             index(field % field, 'VV') /= 0 .or. &
             index(field % field, 'TT') /= 0 .or. &
             index(field % field, 'RH') /= 0 .or. &
             index(field % field, 'GHT') /= 0 .or. &
             index(field % field, 'PMSL') /= 0 .or. &
             index(field % field, 'PSFC') /= 0 .or. &
             index(field % field, 'SOILHGT') /= 0 .or. &
             index(field % field, 'SM000010') /= 0 .or. &
             index(field % field, 'SM010040') /= 0 .or. &
             index(field % field, 'SM040100') /= 0 .or. &
             index(field % field, 'SM100200') /= 0 .or. &
             index(field % field, 'SM010200') /= 0 .or. &
             index(field % field, 'ST000010') /= 0 .or. &
             index(field % field, 'ST010040') /= 0 .or. &
             index(field % field, 'ST040100') /= 0 .or. &
             index(field % field, 'ST100200') /= 0 .or. &
             index(field % field, 'ST010200') /= 0 .or. &
             index(field % field, 'PRES') /= 0 .or. &
             index(field % field, 'SNOW') /= 0 .or. &
             index(field % field, 'SEAICE') /= 0 .or. &
             index(field % field, 'SKINTEMP') /= 0) then

            if (index(field % field, 'SM000010') /= 0 .or. &
                index(field % field, 'SM010040') /= 0 .or. &
                index(field % field, 'SM040100') /= 0 .or. &
                index(field % field, 'SM100200') /= 0 .or. &
                index(field % field, 'SM010200') /= 0 .or. &
                index(field % field, 'ST000010') /= 0 .or. &
                index(field % field, 'ST010040') /= 0 .or. &
                index(field % field, 'ST040100') /= 0 .or. &
                index(field % field, 'ST100200') /= 0 .or. &
                index(field % field, 'ST010200') /= 0 .or. &
                index(field % field, 'SNOW') /= 0 .or. &
                index(field % field, 'SEAICE') /= 0 .or. &
                index(field % field, 'SKINTEMP') /= 0) then
               k = 1
            else if (index(field % field, 'PMSL') == 0 .and. &
                     index(field % field, 'PSFC') == 0 .and. &
                     index(field % field, 'SOILHGT') == 0) then
               do k=1,config_nfglevels
                  if (vert_level(k) == field % xlvl .or. vert_level(k) == -1.0) exit
               end do
               if (k > config_nfglevels) write(0,*) 'ERROR: We seem to have more levels than we thought we should!'
               if (vert_level(k) == -1.0) vert_level(k) = field % xlvl
            else
               k = 1
            end if

            !
            ! Set up projection
            !
            call map_init(proj)
          
            if (field % iproj == PROJ_LATLON) then
               call map_set(PROJ_LATLON, proj, &
                            latinc = real(field % deltalat), &
                            loninc = real(field % deltalon), &
                            knowni = 1.0, &
                            knownj = 1.0, &
                            lat1 = real(field % startlat), &
                            lon1 = real(field % startlon))
            else if (field % iproj == PROJ_GAUSS) then
               call map_set(PROJ_GAUSS, proj, &
                            nlat = nint(field % deltalat), &
                            loninc = real(field % deltalon), &
                            lat1 = real(field % startlat), &
                            lon1 = real(field % startlon))
!                            nxmax = nint(360.0 / field % deltalon), &
            end if


            !
            ! Horizontally interpolate the field at level k
            !
            if (index(field % field, 'UU') /= 0) then
write(0,*) 'Interpolating U at ', k, vert_level(k)

               mask_array => edge_mask

               nInterpPoints = grid % nEdges
               latPoints => grid % latEdge % array
               lonPoints => grid % lonEdge % array
               destField2d => fg % u % array
               ndims = 2
            else if (index(field % field, 'VV') /= 0) then
write(0,*) 'Interpolating V at ', k, vert_level(k)

               mask_array => edge_mask

               nInterpPoints = grid % nEdges
               latPoints => grid % latEdge % array
               lonPoints => grid % lonEdge % array
               destField2d => fg % v % array
               ndims = 2
            else if (index(field % field, 'TT') /= 0) then
write(0,*) 'Interpolating T at ', k, vert_level(k)
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % t % array
               ndims = 2
            else if (index(field % field, 'RH') /= 0) then
write(0,*) 'Interpolating RH at ', k, vert_level(k)
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % rh % array
               ndims = 2
            else if (index(field % field, 'GHT') /= 0) then
write(0,*) 'Interpolating GHT at ', k, vert_level(k)
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % z % array
               ndims = 2
            else if (index(field % field, 'PRES') /= 0) then
write(0,*) 'Interpolating PRES at ', k, vert_level(k)
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % p % array
               ndims = 2
            else if (index(field % field, 'PMSL') /= 0) then
write(0,*) 'Interpolating PMSL'
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField1d => fg % pmsl % array
               ndims = 1
            else if (index(field % field, 'PSFC') /= 0) then
write(0,*) 'Interpolating PSFC'
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField1d => fg % psfc % array
               ndims = 1
            else if (index(field % field, 'SOILHGT') /= 0) then
write(0,*) 'Interpolating SOILHGT'
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField1d => fg % soilz % array
               ndims = 1
            else if (index(field % field, 'SM000010') /= 0) then
write(0,*) 'Interpolating SM000010'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % sm_fg % array
               k = 1
               ndims = 2
               fg % dzs_fg % array(k,:) = 10.
               fg % zs_fg  % array(k,:) = 10.
            else if (index(field % field, 'SM010200') /= 0) then
write(0,*) 'Interpolating SM010200'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % sm_fg % array
               k = 2
               ndims = 2
               fg % dzs_fg % array(k,:) = 200.-10.
               fg % zs_fg  % array(k,:) = 200.
            else if (index(field % field, 'SM010040') /= 0) then
write(0,*) 'Interpolating SM010040'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % sm_fg % array
               k = 2
               ndims = 2
               fg % dzs_fg % array(k,:) = 40.-10.
               fg % zs_fg  % array(k,:) = 40.
            else if (index(field % field, 'SM040100') /= 0) then
write(0,*) 'Interpolating SM040100'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % sm_fg % array
               k = 3
               ndims = 2
               fg % dzs_fg % array(k,:) = 100.-40.
               fg % zs_fg  % array(k,:) = 100.
            else if (index(field % field, 'SM100200') /= 0) then
write(0,*) 'Interpolating SM100200'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 1.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % sm_fg % array
               k = 4
               ndims = 2
               fg % dzs_fg % array(k,:) = 200.-100.
               fg % zs_fg  % array(k,:) = 200.
            else if (index(field % field, 'ST000010') /= 0) then
write(0,*) 'Interpolating ST000010'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % st_fg % array
               k = 1
               ndims = 2
               fg % dzs_fg % array(k,:) = 10.
               fg % zs_fg  % array(k,:) = 10.
            else if (index(field % field, 'ST010200') /= 0) then
write(0,*) 'Interpolating ST010200'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % st_fg % array
               k = 2
               ndims = 2
               fg % dzs_fg % array(k,:) = 200.-10.
               fg % zs_fg  % array(k,:) = 200.
            else if (index(field % field, 'ST010040') /= 0) then
write(0,*) 'Interpolating ST010040'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % st_fg % array
               k = 2
               ndims = 2
               fg % dzs_fg % array(k,:) = 40.-10.               
               fg % zs_fg  % array(k,:) = 40.
            else if (index(field % field, 'ST040100') /= 0) then
write(0,*) 'Interpolating ST040100'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % st_fg % array
               k = 3
               ndims = 2
               fg % dzs_fg % array(k,:) = 100.-40.             
               fg % zs_fg  % array(k,:) = 100.
            else if (index(field % field, 'ST100200') /= 0) then
write(0,*) 'Interpolating ST100200'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 0.0
               masked = 0
               fillval = 285.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField2d => fg % st_fg % array
               k = 4
               ndims = 2
               fg % dzs_fg % array(k,:) = 200.-100.
               fg % zs_fg  % array(k,:) = 200.
            else if (index(field % field, 'SNOW') /= 0) then
write(0,*) 'Interpolating SNOW'

               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = 0

               masked = 0
               fillval = 0.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField1d => fg % snow % array
               ndims = 1
            else if (index(field % field, 'SEAICE') /= 0) then
write(0,*) 'Interpolating SEAICE'

               interp_list(1) = SIXTEEN_POINT
               interp_list(2) = FOUR_POINT
               interp_list(3) = W_AVERAGE4
               interp_list(4) = SEARCH
               interp_list(5) = 0

               maskval = 1.0
               masked = 1
               fillval = 0.0

               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField1d => fg % xice % array
               ndims = 1
            else if (index(field % field, 'SKINTEMP') /= 0) then
write(0,*) 'Interpolating SKINTEMP'
               nInterpPoints = grid % nCells
               latPoints => grid % latCell % array
               lonPoints => grid % lonCell % array
               destField1d => fg % skintemp % array
               ndims = 1
            end if

            allocate(rslab(-3:field % nx+3, field % ny))
            rslab(1:field % nx, 1:field % ny) = field % slab(1:field % nx, 1:field % ny)
            rslab(0, 1:field % ny)  = field % slab(field % nx, 1:field % ny)
            rslab(-1, 1:field % ny) = field % slab(field % nx-1, 1:field % ny)
            rslab(-2, 1:field % ny) = field % slab(field % nx-2, 1:field % ny)
            rslab(field % nx+1, 1:field % ny) = field % slab(1, 1:field % ny)
            rslab(field % nx+2, 1:field % ny) = field % slab(2, 1:field % ny)
            rslab(field % nx+3, 1:field % ny) = field % slab(3, 1:field % ny)

            do i=1,nInterpPoints
               if (mask_array(i) /= masked) then
                  lat = latPoints(i)*DEG_PER_RAD
                  lon = lonPoints(i)*DEG_PER_RAD
                  call latlon_to_ij(proj, lat, lon, x, y)
                  if (x < 0.5) then
                     lon = lon + 360.0
                     call latlon_to_ij(proj, lat, lon, x, y)
                  end if
                  if (ndims == 1) then
                     destField1d(i) = interp_sequence(x, y, 1, rslab, -3, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1, maskval=maskval, mask_array=maskslab)
                  else if (ndims == 2) then
                     destField2d(k,i) = interp_sequence(x, y, 1, rslab, -3, field % nx + 3, 1, field % ny, 1, 1, msgval, interp_list, 1, maskval=maskval, mask_array=maskslab)
                  end if
               else
                  if (ndims == 1) then
                     destField1d(i) = fillval
                  else if (ndims == 2) then
                     destField2d(k,i) = fillval
                  end if
               end if
            end do

            deallocate(rslab)
     
         end if
   
         deallocate(field % slab)
         call read_next_met_field(field, istatus)
      end do

      call read_met_close()

      ! Fix for isobaric data
      if (minval(fg % p % array(:,:)) == 0.0) then
         write(0,*) 'Setting pressure field for isobaric data'
         do k=1,config_nfglevels
            if (vert_level(k) /= 200100.0) then
               fg % p % array(k,:) = vert_level(k)
            else
               fg % p % array(k,:) = fg % psfc % array(:)
            end if
         end do
      end if

      ! Set SST based on SKINTEMP field if it wasn't found in input data
      if (minval(fg % sst % array) == 0.0 .and. maxval(fg % sst % array) == 0.0) then
         write(0,*) 'Setting SST from SKINTEMP'
         where (grid % landmask % array == 0) fg % sst % array = fg % skintemp % array
      end if

      ! Set SNOWC (snow-cover flag) based on SNOW
      fg % snowc % array(:) = 0.0
      where (fg % snow % array > 0.0) fg % snowc % array = 1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!MGD CHECK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do iCell=1,grid%nCells
   if (grid % landmask % array(iCell) == 1) then

      do k = 1, config_nfgsoillevels
         if (fg % st_fg % array(k,iCell) <= 0.0) write(0,*) 'Bad st_fg ', k, iCell
      enddo

      do k = 1, config_nfgsoillevels
         if (fg % sm_fg % array(k,iCell) <= 0.0) write(0,*) 'Bad sm_fg ', k, iCell
      enddo
      !LDF end.

   end if
end do
write(0,*) 'Done with soil consistency check'


      !
      ! Get SEAICE from a separate file
      !
      call read_met_init('SEAICE_FRACTIONAL', .true., config_start_time(1:13), istatus)

      if (istatus /= 0) then
         write(0,*) 'Error reading SEAICE_FRACTIONAL data'
      end if

      if (istatus == 0) then
         call read_next_met_field(field, istatus)
         do while (istatus == 0)
            if (index(field % field, 'SEAICE') /= 0) then

write(0,*) 'PROCESSING SEAICE'

               !
               ! Set up projection
               !
               call map_init(proj)
          
               if (field % iproj == PROJ_PS) then
                  call map_set(PROJ_PS, proj, &
                               dx = real(field % dx * 1000.0), &
                               truelat1 = real(field % truelat1), &
                               stdlon = real(field % xlonc), &
                               knowni = real(field % nx / 2.0), &
                               knownj = real(field % ny / 2.0), &
                               lat1 = real(field % startlat), &
                               lon1 = real(field % startlon))
               end if

               if (index(field % field, 'SEAICE') /= 0) then
                  nInterpPoints = grid % nCells
                  latPoints => grid % latCell % array
                  lonPoints => grid % lonCell % array
                  destField1d => fg % xice % array
                  ndims = 1
               end if
   
               interp_list(1) = FOUR_POINT
               interp_list(2) = W_AVERAGE4
               interp_list(3) = 0
   
               masked = 1
               fillval = 0.0
               msgval = 1.01
               mask_array => grid % landmask % array


               allocate(rslab(field % nx, field % ny))
               rslab(:,:) = field % slab(:,:)
               do i=1,nInterpPoints
                  if (mask_array(i) /= masked) then
                     lat = latPoints(i)*DEG_PER_RAD
                     lon = lonPoints(i)*DEG_PER_RAD
                     call latlon_to_ij(proj, lat, lon, x, y)
                     if (x < 0.5) then
                        lon = lon + 360.0
                        call latlon_to_ij(proj, lat, lon, x, y)
                     end if
                     if (ndims == 1) then
                        destField1d(i) = interp_sequence(x, y, 1, rslab, 1, field % nx, 1, field % ny, 1, 1, msgval, interp_list, 1)
                        if (destField1d(i) == msgval) destField1d(i) = fillval
                     else if (ndims == 2) then
                        destField2d(k,i) = interp_sequence(x, y, 1, rslab, 1, field % nx, 1, field % ny, 1, 1, msgval, interp_list, 1)
                        if (destField2d(k,i) == msgval) destField2d(k,i) = fillval
                     end if
                  else
                     if (ndims == 1) then
                        destField1d(i) = fillval
                     else if (ndims == 2) then
                        destField2d(k,i) = fillval
                     end if
                  end if
               end do
               deallocate(rslab)

            end if
      
            deallocate(field % slab)
            call read_next_met_field(field, istatus)
         end do
      end if

      call read_met_close()

      if (allocated(maskslab)) deallocate(maskslab)

      ! Freeze really cold ocean
      where (fg % sst % array < 271.0 .and. grid % landmask % array == 0) fg % xice % array = 1.0

      ! Set SEAICE (0/1 flag) based on XICE (fractional ice coverage)
      fg % seaice % array(:) = 0.0
      where (fg % xice % array >= 0.5) fg % seaice % array = 1.0


      !
      ! For now, hard-wire soil layer depths and thicknesses
      !

      !LDF begin:
      !fg % dzs % array(1,:) = 0.10
      !fg % dzs % array(2,:) = 0.30
      !fg % dzs % array(3,:) = 0.60
      !fg % dzs % array(4,:) = 1.00

      !fg % dz % array(1,:) = 0.05
      !fg % dz % array(2,:) = 0.25
      !fg % dz % array(3,:) = 0.70
      !fg % dz % array(4,:) = 1.50
      !LDF end.


      !  
      ! Compute normal wind component and store in fg%u
      !  
      do iEdge=1,grid%nEdges
         do k=1,config_nfglevels
            fg % u % array(k,iEdge) = cos(grid % angleEdge % array(iEdge)) * fg % u % array(k,iEdge) &
                                    + sin(grid % angleEdge % array(iEdge)) * fg % v % array(k,iEdge)
         end do
      end do


      !  
      ! Vertically interpolate meteorological data
      !  
      allocate(sorted_arr(2,config_nfglevels))

      do iCell=1,grid%nCells

         ! T
         sorted_arr(:,:) = -999.0
         do k=1,config_nfglevels
            sorted_arr(1,k) = fg % z % array(k,iCell)
!NOSFC            if (vert_level(k) == 200100.0) sorted_arr(1,k) = fg % soilz % array(iCell)
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999.0
            sorted_arr(2,k) = fg % t % array(k,iCell)
         end do
         call quicksort(config_nfglevels, sorted_arr)
         do k=1,grid%nVertLevels
            target_z = 0.5 * (grid % zgrid % array(k,iCell) + grid % zgrid % array(k+1,iCell))
            state % theta_m % array(k,iCell) = vertical_interp(target_z, config_nfglevels, sorted_arr, order=1, extrap=1)
         end do


         ! RH
         sorted_arr(:,:) = -999.0
         do k=1,config_nfglevels
            sorted_arr(1,k) = fg % z % array(k,iCell)
!NOSFC            if (vert_level(k) == 200100.0) sorted_arr(1,k) = fg % soilz % array(iCell)
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999.0
            sorted_arr(2,k) = fg % rh % array(k,iCell)
         end do
         call quicksort(config_nfglevels, sorted_arr)
         do k=1,grid%nVertLevels
            target_z = 0.5 * (grid % zgrid % array(k,iCell) + grid % zgrid % array(k+1,iCell))
            state % scalars % array(state % index_qv,k,iCell) = vertical_interp(target_z, config_nfglevels, sorted_arr, order=1, extrap=0)
         end do


         ! GHT
         sorted_arr(:,:) = -999.0
         do k=1,config_nfglevels
            sorted_arr(1,k) = fg % z % array(k,iCell)
!NOSFC            if (vert_level(k) == 200100.0) sorted_arr(1,k) = fg % soilz % array(iCell)
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999.0
            sorted_arr(2,k) = fg % z % array(k,iCell)
         end do
         call quicksort(config_nfglevels, sorted_arr)
         do k=1,grid%nVertLevels
            target_z = 0.5 * (grid % zgrid % array(k,iCell) + grid % zgrid % array(k+1,iCell))
            fg % gfs_z % array(k,iCell) = vertical_interp(target_z, config_nfglevels, sorted_arr, order=1, extrap=1)
         end do


         ! PRESSURE
         sorted_arr(:,:) = -999.0
         do k=1,config_nfglevels
            sorted_arr(1,k) = fg % z % array(k,iCell)
            if (vert_level(k) == 200100.0) then 
!NOSFC               sorted_arr(1,k) = fg % soilz % array(iCell)
               sorted_arr(1,k) = 99999.0
               sfc_k = k
            end if
            sorted_arr(2,k) = log(fg % p % array(k,iCell))
         end do
         call quicksort(config_nfglevels, sorted_arr)
         do k=1,grid%nVertLevels
            target_z = 0.5 * (grid % zgrid % array(k,iCell) + grid % zgrid % array(k+1,iCell))
            diag % pressure % array(k,iCell) = exp(vertical_interp(target_z, config_nfglevels, sorted_arr, order=1, extrap=1))
         end do


         ! PRESSURE
!         sorted_arr(:,:) = -999.0
!         do k=1,config_nfglevels
!            sorted_arr(1,k) = fg % z % array(k,iCell)
!            if (vert_level(k) == 200100.0) then 
!!NOSFC               sorted_arr(1,k) = fg % soilz % array(iCell)
!               sorted_arr(1,k) = 99999.0
!               sfc_k = k
!            end if
!            sorted_arr(2,k) = log(fg % p % array(k,iCell))
!         end do
!         call quicksort(config_nfglevels, sorted_arr)
!         do k=1,grid%nVertLevels+1
!            target_z = grid % zgrid % array(k,iCell)
!            fg % gfs_p % array(k,iCell) = exp(vertical_interp(target_z, config_nfglevels, sorted_arr, order=1, extrap=1))
!         end do

      end do


      do iEdge=1,grid%nEdges

         ! U
         sorted_arr(:,:) = -999.0
         do k=1,config_nfglevels
            sorted_arr(1,k) = 0.5 * (fg % z % array(k,cellsOnEdge(1,iEdge)) + fg % z % array(k,cellsOnEdge(2,iEdge)))
!NOSFC            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 0.5 * (fg % soilz % array(cellsOnEdge(1,iEdge)) + fg % soilz % array(cellsOnEdge(2,iEdge)))
            if (vert_level(k) == 200100.0) sorted_arr(1,k) = 99999.0
            sorted_arr(2,k) = fg % u % array(k,iEdge)
         end do
         call quicksort(config_nfglevels, sorted_arr)
         do k=1,grid%nVertLevels
            target_z = 0.25 * (grid % zgrid % array(k,cellsOnEdge(1,iEdge)) + grid % zgrid % array(k+1,cellsOnEdge(1,iEdge)) + grid % zgrid % array(k,cellsOnEdge(2,iEdge)) + grid % zgrid % array(k+1,cellsOnEdge(2,iEdge)))
            state % u % array(k,iEdge) = vertical_interp(target_z, config_nfglevels, sorted_arr, order=1, extrap=0)
         end do

      end do


      !
      ! Adjust surface pressure for difference in topography
      !
      do sfc_k=1,config_nfglevels
         if (vert_level(sfc_k) == 200100.) exit
      end do 
      do iCell=1,grid%nCells

         ! We need to extrapolate
            sorted_arr(:,:) = -999.0
            do k=1,config_nfglevels
               sorted_arr(1,k) = fg % z % array(k,iCell)
               if (vert_level(k) == 200100.0) then 
!NOSFC                  sorted_arr(1,k) = fg % soilz % array(iCell)
                  sorted_arr(1,k) = 99999.0
               end if
               sorted_arr(2,k) = log(fg % p % array(k,iCell))
            end do
            call quicksort(config_nfglevels, sorted_arr)
            target_z = grid % zgrid % array(1,iCell)
            fg % psfc % array(iCell) = exp(vertical_interp(target_z, config_nfglevels, sorted_arr, order=1, extrap=1))

      end do

      deallocate(sorted_arr)



      !
      ! Diagnose fields needed in initial conditions file (u, w, rho, theta, scalars)
      ! NB: At this point, "rho_zz" is simple dry density, and "theta_m" is regular potential temperature
      !
      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels

            ! QV
            es = 6.112 * exp((17.27*(state % theta_m % array(k,iCell) - 273.16))/(state % theta_m % array(k,iCell) - 35.86))
            rs = 0.622 * es / (diag % pressure % array(k,iCell) - es)
            scalars(state % index_qv,k,iCell) = rs * scalars(state % index_qv,k,iCell)

            ! PI
            p(k,iCell) = (diag % pressure % array(k,iCell) / p0) ** (rgas / cp)

            ! THETA - can compute this using PI instead
!            t(k,iCell) = t(k,iCell) / p(k,iCell)
            t(k,iCell) = t(k,iCell) * (p0 / diag % pressure % array(k,iCell)) ** (rgas / cp)

            ! RHO_ZZ
            rho_zz(k,iCell) = diag % pressure % array(k,iCell) / rgas / (p(k,iCell) * t(k,iCell))
            rho_zz(k,iCell) = rho_zz(k,iCell) / (1.0 + scalars(state % index_qv,k,iCell))
         end do
      end do


      !
      ! Reference state based on a dry isothermal atmosphere
      !
      do iCell=1,grid % nCells
         do k=1,nz1
            ztemp    = 0.5*(zgrid(k+1,iCell)+zgrid(k,iCell))
            ppb(k,iCell) = p0*exp(-gravity*ztemp/(rgas*t0b))      ! pressure_base
            pb (k,iCell) = (ppb(k,iCell)/p0)**(rgas/cp)           ! exner_base
!            rb (k,iCell) = ppb(k,iCell)/(rgas*t0b*zz(k,iCell))    ! rho_base
            rb (k,iCell) = ppb(k,iCell)/(rgas*t0b)                ! rho_base
            tb (k,iCell) = t0b/pb(k,iCell)                        ! theta_base
            rtb(k,iCell) = rb(k,iCell)*tb(k,iCell)                ! rtheta_base
            p  (k,iCell) = pb(k,iCell)                            ! exner
            pp (k,iCell) = 0.                                     ! pressure_p
            rr (k,iCell) = 0.                                     ! rho_p
         end do
      end do

      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            pp(k,iCell) = diag % pressure % array(k,iCell) - ppb(k,iCell) 
            rr(k,iCell) = rho_zz(k,iCell) - rb(k,iCell)
         end do
      end do

      do iCell=1,grid%nCells
         k = 1
         rho_zz(k,iCell) = ((diag % pressure % array(k,iCell) / p0)**(cv / cp)) * (p0 / rgas) / (t(k,iCell)*(1.0 + 1.61*scalars(state % index_qv,k,iCell)))
         rr(k,iCell) = rho_zz(k,iCell) - rb(k,iCell)

         do k=2,grid % nVertLevels
            it = 0
            p_check = 2.0 * 0.0001
            do while ( (it < 30) .and. (p_check > 0.0001) )

               p_check = pp(k,iCell)
               dz = (zgrid(k,iCell) - zgrid(k-1,iCell))
               pp(k,iCell) = pp(k-1,iCell) - 0.5 * (rr(k,iCell) + rr(k-1,iCell))*gravity*dz &
                                           - 0.5 * (rho_zz(k,iCell)*scalars(state % index_qv,k,iCell) + rho_zz(k-1,iCell)*scalars(state % index_qv,k-1,iCell))*gravity*dz
               diag % pressure % array(k,iCell) = pp(k,iCell) + ppb(k,iCell)
               p(k,iCell) = (diag % pressure % array(k,iCell) / p0) ** (rgas / cp)
               rho_zz(k,iCell) = diag % pressure % array(k,iCell) / rgas / (p(k,iCell)*t(k,iCell)*(1.0 + 1.61*scalars(state % index_qv,k,iCell)))
               rr(k,iCell) = rho_zz(k,iCell) - rb(k,iCell)

               p_check = abs(p_check - pp(k,iCell))
                
               it = it + 1
            end do
         end do
      end do

      ! Compute theta_m and rho-tilde
      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            t(k,iCell) = t(k,iCell) * (1.0 + 1.61*scalars(state % index_qv,k,iCell))
            rho_zz(k,iCell) = rho_zz(k,iCell) / zz(k,iCell)
            rb(k,iCell) = rb(k,iCell) / zz(k,iCell)
         end do
      end do

      do iEdge=1,grid%nEdges
         do k=1,grid%nVertLevels
            diag % ru % array(k,iEdge) = state % u % array(k,iEdge) * 0.5*(state % rho_zz % array(k,cellsOnEdge(1,iEdge)) + state % rho_zz % array(k,cellsOnEdge(2,iEdge)))
         end do
      end do


      diag % rw % array = 0.
      state % w % array = 0.
      do iEdge = 1,grid % nEdges

         cell1 = CellsOnEdge(1,iEdge)
         cell2 = CellsOnEdge(2,iEdge)

         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then
         do k = 2, grid%nVertLevels
            flux =  (fzm(k)*diag % ru % array(k,iEdge)+fzp(k)*diag % ru % array(k-1,iEdge))
            diag % rw % array(k,cell2) = diag % rw % array(k,cell2) + zf(k,2,iEdge)*flux
            diag % rw % array(k,cell1) = diag % rw % array(k,cell1) - zf(k,1,iEdge)*flux

            if (config_theta_adv_order ==3) then 
               diag % rw % array(k,cell2) = diag % rw % array(k,cell2)    &
                                            - sign(1.,diag % ru % array(k,iEdge))*config_coef_3rd_order*zf3(k,2,iEdge)*flux
               diag % rw % array(k,cell1) = diag % rw % array(k,cell1)    &
                                            + sign(1.,diag % ru % array(k,iEdge))*config_coef_3rd_order*zf3(k,1,iEdge)*flux
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
   
      deallocate(vert_level)

      end if    ! config_met_interp

     
      ! Calculate surface pressure
      do iCell=1,grid%nCells
         diag % surface_pressure % array(iCell) = 0.5*gravity/rdzw(1)                                              &
                                                * (1.25* rho_zz(1,iCell) * (1. + scalars(state % index_qv, 1, iCell))  &
                                                -  0.25* rho_zz(2,iCell) * (1. + scalars(state % index_qv, 2, iCell)))
         diag % surface_pressure % array(iCell) = diag % surface_pressure % array(iCell) + pp(1,iCell) + ppb(1,iCell)
      end do

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            diag % rho % array(k,iCell) = state % rho_zz % array(k,iCell) * zz(k,iCell)
            diag % theta % array(k,iCell) = state % theta_m % array(k,iCell) / (1.0 + 1.61 * scalars(state % index_qv,k,iCell))
         end do
      end do

   end subroutine init_atm_test_case_gfs

   subroutine init_atm_test_case_sst(domain, dminfo, grid, fg, state, diag, test_case, parinfo)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Real-data test case using SST data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use mpas_dmpar
      use mpas_io_output
      use init_atm_read_met
      use init_atm_llxy
      use init_atm_hinterp

      implicit none

      type (domain_type), intent(inout) :: domain
      type (dm_info), intent(in) :: dminfo
      type (mesh_type), intent(inout) :: grid
      type (fg_type), intent(inout) :: fg
      type (state_type), intent(inout) :: state
      type (diag_type), intent(inout) :: diag
      integer, intent(in) :: test_case
      type (parallel_info), pointer :: parinfo

      integer :: istatus
      integer :: iCell, i, j
      type (met_data) :: field
      type (proj_info) :: proj
      real (kind=RKIND) :: lat, lon, x, y
      integer, dimension(5) :: interp_list
      real (kind=RKIND), allocatable, dimension(:,:) :: slab_r8
      type (io_output_object) :: sfc_update_obj
      type (MPAS_Clock_type) :: fg_clock
      type (MPAS_Time_type) :: start_time, stop_time, curr_time
      type (MPAS_TimeInterval_type) :: fg_interval
      character (len=32) :: timeString


      ! Set interpolation sequence to be used for SST/SKINTEMP field
      interp_list(1) = FOUR_POINT
      interp_list(2) = SEARCH
      interp_list(3) = 0


      ! Set up clock to step through all intermediate file dates to be processed
      call mpas_set_time(start_time, dateTimeString=trim(config_start_time))
      call mpas_set_time(stop_time,  dateTimeString=trim(config_stop_time))
      call mpas_set_timeInterval(fg_interval, S=config_fg_interval)

      call mpas_create_clock(fg_clock, start_time, fg_interval, stopTime=stop_time) 


      ! Initialize the output file
      sfc_update_obj % time = 1
      sfc_update_obj % filename = trim(config_sfc_update_name)

      call mpas_output_state_init(sfc_update_obj, domain, "SFC")

      ! Loop over all times, interpolating the SST/SKINTEMP field from each intermediate file
      curr_time = mpas_get_clock_time(fg_clock, MPAS_NOW) 
      do while (curr_time <= stop_time)
         call mpas_get_time(curr_time, dateTimeString=timeString)
         write(0,*) 'Processing ',trim(config_met_prefix)//':'//timeString(1:13)

         ! Open intermediate file
         call read_met_init(trim(config_met_prefix), .false., timeString(1:13), istatus)
         if (istatus /= 0) then
            write(0,*) 'Error reading ',trim(config_met_prefix)//':'//timeString(1:13)
            exit
         end if

         ! Scan through all fields in the file, looking for the SST or SKINTEMP field
         call read_next_met_field(field, istatus)
         do while (istatus == 0)

            if (index(field % field, 'SKINTEMP') /= 0 .or. index(field % field, 'SST') /= 0) then

               ! Interpolation routines use real(kind=RKIND), so copy from default real array
               allocate(slab_r8(field % nx, field % ny))
               do j=1,field % ny
               do i=1,field % nx
                  slab_r8(i,j) = field % slab(i,j)
               end do
               end do

               !
               ! Set up map projection
               !
               call map_init(proj)
             
               if (field % iproj == PROJ_LATLON) then
                  call map_set(PROJ_LATLON, proj, &
                               latinc = real(field % deltalat), &
                               loninc = real(field % deltalon), &
                               knowni = 1.0, &
                               knownj = 1.0, &
                               lat1 = real(field % startlat), &
                               lon1 = real(field % startlon))
               else if (field % iproj == PROJ_GAUSS) then
                  call map_set(PROJ_GAUSS, proj, &
                               nlat = nint(field % deltalat), &
                               loninc = real(field % deltalon), &
                               lat1 = real(field % startlat), &
                               lon1 = real(field % startlon))
!                               nxmax = nint(360.0 / field % deltalon), &
               end if
   
               ! Interpolate SST/SKINTEMP field to each MPAS grid cell
               do iCell=1,grid % nCells
                  lat = grid % latCell % array(iCell) * DEG_PER_RAD
                  lon = grid % lonCell % array(iCell) * DEG_PER_RAD
                  call latlon_to_ij(proj, lat, lon, x, y)
                  if (y < 0.5) then
                     y = 1.0
                  else if (y >= real(field%ny)+0.5) then
                     y = real(field % ny)
                  end if
                  if (x < 0.5) then
                     lon = lon + 360.0
                     call latlon_to_ij(proj, lat, lon, x, y)
                  else if (x >= real(field%nx)+0.5) then
                     lon = lon - 360.0
                     call latlon_to_ij(proj, lat, lon, x, y)
                  end if
                  fg % sst % array(iCell) = interp_sequence(x, y, 1, slab_r8, 1, field % nx, 1, field % ny, 1, 1, -1.e30, interp_list, 1)
               end do

               deallocate(slab_r8)
               deallocate(field % slab)
               exit
            end if

            if (index(field % field, 'SEAICE') /= 0) then

               ! Interpolation routines use real(kind=RKIND), so copy from default real array
               allocate(slab_r8(field % nx, field % ny))
               do j=1,field % ny
               do i=1,field % nx
                  slab_r8(i,j) = field % slab(i,j)
               end do
               end do

               !
               ! Set up map projection
               !
               call map_init(proj)
             
               if (field % iproj == PROJ_LATLON) then
                  call map_set(PROJ_LATLON, proj, &
                               latinc = real(field % deltalat), &
                               loninc = real(field % deltalon), &
                               knowni = 1.0, &
                               knownj = 1.0, &
                               lat1 = real(field % startlat), &
                               lon1 = real(field % startlon))
               else if (field % iproj == PROJ_GAUSS) then
                  call map_set(PROJ_GAUSS, proj, &
                               nlat = nint(field % deltalat), &
                               loninc = real(field % deltalon), &
                               lat1 = real(field % startlat), &
                               lon1 = real(field % startlon))
!                               nxmax = nint(360.0 / field % deltalon), &
               end if
   
               ! Interpolate SEAICE/SKINTEMP field to each MPAS grid cell
               do iCell=1,grid % nCells
                  lat = grid % latCell % array(iCell) * DEG_PER_RAD
                  lon = grid % lonCell % array(iCell) * DEG_PER_RAD
                  call latlon_to_ij(proj, lat, lon, x, y)
                  if (y < 0.5) then
                     y = 1.0
                  else if (y >= real(field%ny)+0.5) then
                     y = real(field % ny)
                  end if
                  if (x < 0.5) then
                     lon = lon + 360.0
                     call latlon_to_ij(proj, lat, lon, x, y)
                  else if (x >= real(field%nx)+0.5) then
                     lon = lon - 360.0
                     call latlon_to_ij(proj, lat, lon, x, y)
                  end if
                  fg % xice % array(iCell) = interp_sequence(x, y, 1, slab_r8, 1, field % nx, 1, field % ny, 1, 1, -1.e30, interp_list, 1)

               end do

               deallocate(slab_r8)
               deallocate(field % slab)
               exit
            end if

            deallocate(field % slab)
            call read_next_met_field(field, istatus)
         end do

         ! Close intermediate file
         call read_met_close()

         ! Write the interpolated SST/SKINTEMP field as a new time slice in the MPAS output file
         call mpas_output_state_for_domain(sfc_update_obj, domain, sfc_update_obj % time)
         sfc_update_obj % time = sfc_update_obj % time + 1

         call mpas_advance_clock(fg_clock)
         curr_time = mpas_get_clock_time(fg_clock, MPAS_NOW) 

         call mpas_get_time(curr_time, dateTimeString=timeString)
         state % xtime % scalar = timeString

      end do

      call mpas_output_state_finalize(sfc_update_obj, dminfo)
      
   end subroutine init_atm_test_case_sst

#if 0
   real function four_pt(nx, ny, array, xx, yy)

      implicit none

      integer, intent(in) :: nx, ny 
      real (kind=4), dimension(nx, ny), intent(in) :: array
      real (kind=4), intent(in) :: xx, yy

      integer :: min_x, max_x, min_y, max_y

      min_x = floor(xx)
      min_y = floor(yy)
      max_x = ceiling(xx)
      max_y = ceiling(yy)

      if (min_x == 0) min_x = max_x
      if (max_x == nx+1) max_x = min_x
      if (min_y == 0) min_y = max_y
      if (max_y == ny+1) max_y = min_y

      if ((min_x < 1) .or. (max_x > nx) .or. (min_y < 1) .or. (max_y > ny)) then
         write(0,*) '(x,y) location out of bounds'
         four_pt = 0.0
         return
      end if 

      if (min_x == max_x) then
         if (min_y == max_y) then
            four_pt = array(min_x,min_y)
         else
            four_pt = array(min_x,min_y)*(real(max_y)-yy) + &
                      array(min_x,max_y)*(yy-real(min_y))
         end if
      else if (min_y == max_y) then
         if (min_x == max_x) then
            four_pt = array(min_x,min_y)
         else
            four_pt = array(min_x,min_y)*(real(max_x)-xx) + &
                      array(max_x,min_y)*(xx-real(min_x))
         end if
      else
         four_pt = (yy - min_y) * (array(min_x,max_y)*(real(max_x)-xx) + &
                                   array(max_x,max_y)*(xx-real(min_x))) + &
                   (max_y - yy) * (array(min_x,min_y)*(real(max_x)-xx) + &
                                   array(max_x,min_y)*(xx-real(min_x)));
      end if

      return

   end function four_pt
#endif


   real function sphere_distance(lat1, lon1, lat2, lon2, radius)
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


   integer function nearest_cell(target_lat, target_lon, &
                                 start_cell, &
                                 nCells, maxEdges, nEdgesOnCell, cellsOnCell, latCell, lonCell)

      implicit none

      real (kind=RKIND), intent(in) :: target_lat, target_lon
      integer, intent(in) :: start_cell
      integer, intent(in) :: nCells, maxEdges
      integer, dimension(nCells), intent(in) :: nEdgesOnCell
      integer, dimension(maxEdges,nCells), intent(in) :: cellsOnCell
      real (kind=RKIND), dimension(nCells), intent(in) :: latCell, lonCell

      integer :: i
      integer :: iCell
      integer :: current_cell
      real (kind=RKIND) :: current_distance, d
      real (kind=RKIND) :: nearest_distance

      nearest_cell = start_cell
      current_cell = -1

      do while (nearest_cell /= current_cell)
         current_cell = nearest_cell
         current_distance = sphere_distance(latCell(current_cell), lonCell(current_cell), target_lat, target_lon, 1.0)
         nearest_cell = current_cell
         nearest_distance = current_distance
         do i = 1, nEdgesOnCell(current_cell)
            iCell = cellsOnCell(i,current_cell)
            d = sphere_distance(latCell(iCell), lonCell(iCell), target_lat, target_lon, 1.0)
            if (d < nearest_distance) then
               nearest_cell = iCell
               nearest_distance = d
            end if
         end do
      end do

   end function nearest_cell


   integer function nearest_edge(target_lat, target_lon, &
                                 start_edge, &
                                 nCells, nEdges, maxEdges, nEdgesOnCell, edgesOnCell, cellsOnEdge, latCell, lonCell, latEdge, lonEdge)

      implicit none

      real (kind=RKIND), intent(in) :: target_lat, target_lon
      integer, intent(in) :: start_edge
      integer, intent(in) :: nCells, nEdges, maxEdges
      integer, dimension(nCells), intent(in) :: nEdgesOnCell
      integer, dimension(maxEdges,nCells), intent(in) :: edgesOnCell
      integer, dimension(2,nEdges), intent(in) :: cellsOnEdge
      real (kind=RKIND), dimension(nCells), intent(in) :: latCell, lonCell
      real (kind=RKIND), dimension(nEdges), intent(in) :: latEdge, lonEdge

      integer :: i, cell1, cell2, iCell
      integer :: iEdge
      integer :: current_edge
      real (kind=RKIND) :: cell1_dist, cell2_dist
      real (kind=RKIND) :: current_distance, d
      real (kind=RKIND) :: nearest_distance

      nearest_edge = start_edge
      current_edge = -1

      do while (nearest_edge /= current_edge)
         current_edge = nearest_edge
         current_distance = sphere_distance(latEdge(current_edge), lonEdge(current_edge), target_lat, target_lon, 1.0)
         nearest_edge = current_edge
         nearest_distance = current_distance
         cell1 = cellsOnEdge(1,current_edge)
         cell2 = cellsOnEdge(2,current_edge)
         cell1_dist = sphere_distance(latCell(cell1), lonCell(cell1), target_lat, target_lon, 1.0)
         cell2_dist = sphere_distance(latCell(cell2), lonCell(cell2), target_lat, target_lon, 1.0)
         if (cell1_dist < cell2_dist) then
            iCell = cell1
         else
            iCell = cell2
         end if
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            d = sphere_distance(latEdge(iEdge), lonEdge(iEdge), target_lat, target_lon, 1.0)
            if (d < nearest_distance) then
               nearest_edge = iEdge
               nearest_distance = d
            end if
         end do
      end do

   end function nearest_edge


   real function vertical_interp(target_z, nz, zf, order, extrap, surface_val, sealev_val)

      implicit none

      real (kind=RKIND), intent(in) :: target_z
      integer, intent(in) :: nz 
      real (kind=RKIND), dimension(2,nz), intent(in) :: zf      ! zf(1,:) is column of vertical coordinate values, zf(2,:) is column of field values
      integer, intent(in), optional :: order
      integer, intent(in), optional :: extrap
      real (kind=RKIND), intent(in), optional :: surface_val
      real (kind=RKIND), intent(in), optional :: sealev_val

      integer :: k, lm, lp
      real (kind=RKIND) :: wm, wp
      real (kind=RKIND) :: slope

      integer :: interp_order, extrap_type
      real (kind=RKIND) :: surface, sealevel


      if (present(order)) then
         interp_order = order
      else
         interp_order = 2
      end if

      if (present(extrap)) then
         extrap_type = extrap
      else
         interp_order = 1
      end if

      if (present(surface_val)) then
         surface = surface_val
      else
         surface = 200100.0
      end if

      if (present(sealev_val)) then
         sealevel = sealev_val
      else
         sealevel = 201300.0
      end if

      !
      ! Extrapolation required
      !
      if (target_z < zf(1,1)) then
         if (extrap_type == 0) then
            vertical_interp = zf(2,1)
         else if (extrap_type == 1) then
            slope = (zf(2,2) - zf(2,1)) / (zf(1,2) - zf(1,1))
            vertical_interp = zf(2,1) + slope * (target_z - zf(1,1))
         end if
         return
      end if
      if (target_z >= zf(1,nz)) then
         if (extrap_type == 0) then
            vertical_interp = zf(2,nz)
         else if (extrap_type == 1) then
            slope = (zf(2,nz) - zf(2,nz-1)) / (zf(1,nz) - zf(1,nz-1))
            vertical_interp = zf(2,nz) + slope * (target_z - zf(1,nz))
         end if
         return
      end if


      !
      ! No extrapolation required
      !
      do k=1,nz-1
         if (target_z >= zf(1,k) .and. target_z < zf(1,k+1)) then
            lm = k
            lp = k+1
            wm = (zf(1,k+1) - target_z) / (zf(1,k+1) - zf(1,k))
            wp = (target_z - zf(1,k)) / (zf(1,k+1) - zf(1,k))
            exit
         end if
      end do

      vertical_interp = wm*zf(2,lm) + wp*zf(2,lp)

      return

   end function vertical_interp


   subroutine init_atm_check_read_error(istatus, fname, dminfo)

      implicit none

      integer, intent(in) :: istatus
      character (len=*), intent(in) :: fname
      type (dm_info), intent(in) :: dminfo

      if (istatus /= 0) then
         write(0,*) 'ERROR: Could not read file '//trim(fname)
         call mpas_dmpar_abort(dminfo)
      end if

   end subroutine init_atm_check_read_error

end module init_atm_test_cases
