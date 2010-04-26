module test_cases

   use grid_types
   use configure
   use constants


   contains


   subroutine setup_hyd_test_case(domain)
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
         write(0,*) ' need hydrostatic test case configuration, error stop '
         stop

      else if ((config_test_case == 1) .or. (config_test_case == 2) .or. (config_test_case == 3)) then
         write(0,*) ' Jablonowski and Williamson baroclinic wave test case '
         if (config_test_case == 1) write(0,*) ' no initial perturbation '
         if (config_test_case == 2) write(0,*) ' initial perturbation included '
         if (config_test_case == 3) write(0,*) ' normal-mode perturbation included '
         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call hyd_test_case_1(block_ptr % mesh, block_ptr % time_levs(1) % state, config_test_case)
            do i=2,nTimeLevs
               call copy_state(block_ptr % time_levs(1) % state, block_ptr % time_levs(i) % state)
            end do

            block_ptr => block_ptr % next
         end do

      else
         write(0,*) ' Only test case 1 and 2 are currently supported for hydrostatic core '
         stop
      end if

   end subroutine setup_hyd_test_case

!----------------------------------------------------------------------------------------------------------

   subroutine hyd_test_case_1(grid, state, test_case)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup baroclinic wave test case from Jablonowski and Williamson 2008 (QJRMS)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_meta), intent(inout) :: grid
      type (grid_state), intent(inout) :: state
      integer, intent(in) :: test_case

      real (kind=RKIND), parameter :: u0 = 35.0
      real (kind=RKIND), parameter :: alpha_grid = 0.  ! no grid rotation
      real (kind=RKIND), parameter :: omega_e = 7.29212e-05
      real (kind=RKIND), parameter :: t0 = 288., delta_t = 4.8e+05, dtdz = 0.005, eta_t = 0.2
      real (kind=RKIND), parameter :: u_perturbation = 1., pert_radius = 0.1, latitude_pert = 40., longitude_pert = 20.
      real (kind=RKIND), parameter :: theta_c = pii/4.0
      real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
      real (kind=RKIND), parameter :: rh_max = 0.4       ! Maximum relative humidity
      real (kind=RKIND), parameter :: k_x = 9.           ! Normal mode wave number

      real (kind=RKIND), dimension(:), pointer :: rdnu, rdnw, fnm, fnp, dbn, dnu, dnw
      real (kind=RKIND), dimension(:), pointer :: surface_pressure
      real (kind=RKIND), dimension(:,:), pointer :: pressure, theta, alpha, geopotential, h
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      integer :: iCell, iEdge, vtx1, vtx2, ivtx, k, nz, nz1
      real (kind=RKIND) :: u, v, flux, fluxk, lat1, lat2, eta_v, r_pert, u_pert, lat_pert, lon_pert, r

      real (kind=RKIND) :: ptop, p0, phi
      real (kind=RKIND) :: lon_Edge

      real (kind=RKIND), dimension(grid % nVertLevels, grid % nCells) :: rel_hum, temperature
      real (kind=RKIND) :: ptmp, es, qvs
      integer :: iter

!      real (kind=RKIND), dimension(grid % nVertLevelsP1 ) :: hyai, hybi, znu, znw, znwc, znwv, hyam, hybm
!      real (kind=RKIND), dimension(grid % nVertLevelsP1 ) :: znuc, znuv, bn, divh, dpn, teta, phi
      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: hyai, hybi, znu, znw, znwc, znwv, hyam, hybm
      real (kind=RKIND), dimension(grid % nVertLevels + 1 ) :: znuc, znuv, bn, divh, dpn, teta

      real (kind=RKIND) :: HYAI_CAM26(27), HYBI_CAM26(27), HYAM_CAM26(26), HYBM_CAM26(26)

      logical, parameter :: cam26 = .true.

      data hyai_cam26 / 0.002194067, 0.004895209, 0.009882418, 0.018052010,  &
                        0.029837240, 0.044623340, 0.061605870, 0.078512430,  &
                        0.077312710, 0.075901310, 0.074240860, 0.072287440,  &
                        0.069989330, 0.067285740, 0.064105090, 0.060363220,  &
                        0.055961110, 0.050782250, 0.044689600, 0.037521910,  &
                        0.029089490, 0.020847390, 0.013344430, 0.007084990,  &
                        0.002521360, 0.000000000, 0.000000000/,              &
           hybi_cam26 / 0.000000000, 0.000000000, 0.000000000, 0.000000000,  &
                        0.000000000, 0.000000000, 0.000000000, 0.000000000,  &
                        0.015053090, 0.032762280, 0.053596220, 0.078106270,  &
                        0.106941100, 0.140863700, 0.180772000, 0.227722000,  &
                        0.282956200, 0.347936400, 0.424382200, 0.514316800,  &
                        0.620120200, 0.723535500, 0.817676800, 0.896215300,  &
                        0.953476103, 0.985112200, 1.000000000/

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
      
      nz1 = grid % nVertLevels
      nz = nz1 + 1

      rdnu => grid % rdnu % array
      rdnw => grid % rdnw % array
      fnm => grid % fnm % array
      fnp => grid % fnp % array
      dbn => grid % dbn % array
      dnu => grid % dnu % array
      dnw => grid % dnw % array

      surface_pressure => state % surface_pressure % array
      pressure => state % pressure % array
      theta => state % theta % array
      alpha => state % alpha % array
      geopotential => state % geopotential % array
      h => state % h % array
      scalars => state % scalars % array

      scalars(:,:,:) = 0.

      p0      = 100000.
      bn (1) = 1.
      znw(1) = 1.
      znwc(1) = 1.
      !znwv(1) = (znwc(1)-.252)*pii/2.
      znwv(1) = ((znwc(1)-.252)*pii/2.*p0-ptop)/(p0-ptop)
		
      if (cam26) then

        if (grid % nVertLevels /= 26 ) then
          write(0,*) ' init is for 26 levels only, error stop '
          stop
        else
      	  do k=1,nz
            hyai(k) = hyai_cam26(k)
            hybi(k) = hybi_cam26(k)
          end do

          write(0,*) ' initialization using cam 26 levels '

        end if

        ptop    = hyai(1)*p0

        do k=1,nz1
           hyam(k) = .5*(hyai(k)+hyai(k+1))
           hybm(k) = .5*(hybi(k)+hybi(k+1))
           znuc(nz-k)   = hyam(k) + hybm(k)
           znwc(nz-k+1) = hyai(k) + hybi(k)
           znu (nz-k  ) = (znuc(nz-k  )*p0-ptop)/(p0-ptop)
           znw (nz-k+1) = (znwc(nz-k+1)*p0-ptop)/(p0-ptop)
!           znuv(nz-k  ) = (znuc(nz-k  )-.252)*pii/2.
!           znwv(nz-k+1) = (znwc(nz-k+1)-.252)*pii/2.
           bn(k+1) = hybi(nz-k)
        end do

      else ! analytic profile

        ptop = 219.4067
        znw(1) = 1.

        do k=1,nz1

          ! eta profile (constant deta for exp=1,)

          znw(k+1) = (1.-float(k)/float(nz1))**2.

          ! profile for tranisition from sigma to presure coordinate
          ! bn(k)=znw(k) for sigma coord, bn(k)=0 for p coord
          !  bn(1)=1, bn(nz)=0 must be satisfied
				
          bn(k+1) = znw(k+1)*sin(.5*pii*znw(k+1))**2
          !!  bn(k+1) = znw(k+1)
																
          znu (k)   = .5*(znw(k)+znw(k+1))
          znuc(k)   = (znu(k  )*(p0-ptop)+ptop)/p0
          znwc(k+1) = (znw(k+1)*(p0-ptop)+ptop)/p0
        end do

      end if  ! cam or analytic grid-level profile

      !
      !  metrics for vertical stretching
      !

      do k=1,nz1
        !znuv(k  ) = (znuc(k  )-.252)*pii/2.
        !znwv(k+1) = (znwc(k+1)-.252)*pii/2.
        znuv(k  ) = ((znuc(k  )-.252)*pii/2.*p0-ptop)/(p0-ptop)
        znwv(k+1) = ((znwc(k+1)-.252)*pii/2.*p0-ptop)/(p0-ptop)
        dnw (k) = znw(k+1)-znw(k)
        rdnw(k) = 1./dnw(k)
        dbn (k) = rdnw(k)*(bn(k+1)-bn(k))
        dpn (k) = 0.
        divh(k) = 0.
        write (6,*) k,znw(k),dnw(k),bn(k),dbn(k)
      end do

      dpn(nz)=0.
      fnm(1) = 0.
      fnp(1) = 0.
      do k=2,nz1
         dnu (k)  = .5*(dnw(k)+dnw(k-1))
         rdnu(k)  = 1./dnu(k)
         fnp (k)  = .5* dnw(k  )/dnu(k)
         fnm (k)  = .5* dnw(k-1)/dnu(k)
      end do

      !
      ! Initialize wind field
      !

      lat_pert = latitude_pert*pii/180.
      lon_pert = longitude_pert*pii/180.

      do iEdge=1,grid % nEdges

         vtx1 = grid % VerticesOnEdge % array (1,iEdge)
         vtx2 = grid % VerticesOnEdge % array (2,iEdge)
         lat1 = grid%latVertex%array(vtx1)
         lat2 = grid%latVertex%array(vtx2)
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


         do k=1,grid % nVertLevels
           fluxk = u0*flux*(cos(znuv(k))**(1.5))
           state % u % array(k,iEdge) = fluxk + u_pert
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! To get hydrostatic balance with misture -- soln. 2.
! original scheme by Jablonowski
!                            T' = -1./R_d *(p/p_0) * d(phi')/d(eta)
!                               = -1./R_d * eta * d(phi')/d(eta)
! soln. 2 -> derive temperature profile from hydrostatic balance with moisture
! 
!                           T_v = -1/(1+q_v)*(p/R_d)* d(eta)/d(p_d) * d(phi)/d(eta)
!                           phi'(k) = phi(k+1) + d(eta)* alpha_pert * d(eta)/d(p_d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
      do iCell=1,grid % nCells

        phi = grid % latCell % array (iCell)

        surface_pressure(iCell) = p0

        do k=1,nz1
          h(k,iCell) = (1.-dbn(k))*(p0-ptop)+dbn(k)*(surface_pressure(iCell)-ptop)
        end do
			
         pressure(nz,iCell) = ptop
         do k=nz1,1,-1
               pressure(k,iCell) = pressure(k+1,iCell)-dnw(k)*h(k,iCell)
         end do

         do k=1,nz1
            ptmp = 0.5*(pressure(k,iCell)+pressure(k+1,iCell))
            if (znuc(k) >= eta_t) then
               teta(k) = t0*(ptmp/p0)**(rgas*dtdz/gravity)
            else
               teta(k) = t0*(ptmp/p0)**(rgas*dtdz/gravity) + delta_t*(eta_t-ptmp/p0)**5
            end if
            theta (k,iCell) = teta(k)+.75*ptmp/h(k,iCell)*pii*u0/rgas*sin(znuv(k))    &
                              *sqrt(cos(znuv(k)))*                         &
                                ((-2.*sin(phi)**6                          &
                                     *(cos(phi)**2+1./3.)+10./63.)         &
                                     *2.*u0*cos(znuv(k))**1.5              &
                                +(1.6*cos(phi)**3                          &
                                     *(sin(phi)**2+2./3.)-pii/4.)*a*omega)

            theta (k,iCell) = theta(k,iCell)*  &
                      (0.5*(pressure(k,iCell)+pressure(k+1,iCell))/p0)**(-rgas/cp)
            alpha(k,iCell) = ((rgas/p0)*theta(k,iCell)* &
                      (0.5*(pressure(k,iCell)+pressure(k+1,iCell))/p0)**cvpm)

         end do
      end do
!
!     initialization for geopotential
!
      do iCell=1,grid % nCells

         phi = grid % latCell % array (iCell)

         geopotential(1,iCell) = u0*cos(znwv(1))**1.5*                     &
                                 ((-2.*sin(phi)**6                     &
                                      *(cos(phi)**2+1./3.)+10./63.)    &
                                      *(u0)*cos(znwv(1))**1.5          &
                                 +(1.6*cos(phi)**3                     &
                                     *(sin(phi)**2+2./3.)-pii/4.)*a*omega)
         do k=1,nz1
           geopotential(k+1,iCell) = geopotential(k,iCell)-dnw(k)*h(k,iCell)*alpha(k,iCell)
         end do
      end do
		
      write(6,*) 'ptop_dry = ',ptop,'  zt_dry = ',geopotential(nz,1)/gravity

      write(6,*) ' full sounding for dry'
      do k=1,nz1
         write(6,*) k, geopotential(k,1)/gravity, 0.01*pressure(k,1), theta(k,1), &
                    theta(k,1)*(pressure(k,1)/p0)**(rgas/cp)
      end do

!
!     initialization for moisture 
!
      if (config_mp_physics /= 0) then

         do iCell=1,grid % nCells
            do k=1,nz1
               ptmp = 0.5*(pressure(k,iCell) + pressure(k+1,iCell))
               if (ptmp < 50000.) then
                  rel_hum(k,iCell) = 0.0
               else
                  rel_hum(k,iCell) = (1.-((p0-ptmp)/50000.)**1.25)
               end if
               rel_hum(k,iCell) = min(rh_max,rel_hum(k,iCell))
            end do
         end do
      else
         rel_hum(:,:) = 0.
      end if


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! iteration 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      do iter=1,30
         do iCell=1,grid % nCells 
   
            phi = grid % latCell % array (iCell)
            do k=1,nz1
               ptmp = 0.5*(pressure(k+1,iCell)+pressure(k,iCell))
   
               if(znuc(k) >= eta_t)  then
                  teta(k) = t0*(ptmp/p0)**(rgas*dtdz/gravity)
               else
                  teta(k) = t0*(ptmp/p0)**(rgas*dtdz/gravity) + delta_t*(eta_t-ptmp/p0)**5
               end if
   
               temperature (k,iCell) = teta(k)+.75*ptmp/h(k,iCell)*pii*u0/rgas*sin(znuv(k))    &
                                 *sqrt(cos(znuv(k)))*                         &
                                   ((-2.*sin(phi)**6                          &
                                        *(cos(phi)**2+1./3.)+10./63.)         &
                                        *2.*u0*cos(znuv(k))**1.5              &
                                   +(1.6*cos(phi)**3                          &
                                        *(sin(phi)**2+2./3.)-pii/4.)*a*omega)
   
               temperature(k,iCell) = temperature(k,iCell)/(1.+0.61*scalars(index_qv,k,iCell))
   
               theta (k,iCell) = temperature(k,iCell)*  &
                                     (0.5*(pressure(k,iCell)+pressure(k+1,iCell))/p0)**(-rgas/cp)
               alpha (k,iCell) = (rgas/p0)*theta(k,iCell)*(1.+1.61*scalars(index_qv,k,iCell))* &
                                     (0.5*(pressure(k,iCell)+pressure(k+1,iCell))/p0)**cvpm
   
               if (temperature(k,iCell) > 273.15) then
                   es  = 1000.*0.6112*exp(17.67*(temperature(k,iCell)-273.15)/(temperature(k,iCell)-29.65))
               else
                   es  = 1000.*0.6112*exp(21.8745584*(temperature(k,iCell)-273.16)/(temperature(k,iCell)-7.66))
               end if
               qvs = (287.04/461.6)*es/(ptmp-es)
!               qvs =  380.*exp(17.27*(temperature(k,iCell)-273.)/(temperature(k,iCell)-36.))/ptmp
   
               scalars(index_qv,k,iCell) = rel_hum(k,iCell)*qvs
            end do
   
            do k=nz1,1,-1
               pressure(k,iCell) = pressure(k+1,iCell)-dnw(k)*h(k,iCell)*(1.+scalars(index_qv,k,iCell))
               geopotential(k,iCell) = geopotential(k+1,iCell)+dnw(k)*h(k,iCell)*alpha(k,iCell)
            end do
   
         end do
      end do

      write(6,*) 'ptop = ',ptop,'  zt = ',geopotential(nz,1)/gravity

      write(6,*) ' full sounding with moisture'
      do k=1,nz1
         write(6,*) k, geopotential(k,1)/gravity, 0.01*pressure(k,1), theta(k,1), &
                    theta(k,1)*(pressure(k,1)/p0)**(rgas/cp)
      end do

! When initializing a scalar, be sure not to put unreasonably large values
! into indices in the moist class
!      scalars(2,:,:) = 1.  ! transport test
!      scalars(2,:,:) = theta  ! transport test
!      if (num_scalars >= 2) then
!         scalars(2,:,:) = 0.0
!         do iCell=1,grid%nCells
!            r = sphere_distance(theta_c, lambda_c, grid % latCell % array(iCell), grid % lonCell % array(iCell), a)
!            if (r < a/3.0) then
!               do k=1,grid%nVertLevels
!                  scalars(2,k,iCell) = (1.0 / 2.0) * (1.0 + cos(pii*r*3.0/a))
!               end do
!            end if
!         end do
!      end if
!      if (num_scalars >= 3) scalars(3,:,:) = theta + 100.  ! transport test
!      if (num_scalars >= 4) scalars(4,:,:) = theta + 200.  ! transport test
!      if (num_scalars >= 5) scalars(5,:,:) = theta + 300.  ! transport test

   end subroutine hyd_test_case_1


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


   real function AA(theta)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! A, used in height field computation for Rossby-Haurwitz wave
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      real (kind=RKIND), intent(in) :: theta

      AA = 0.5 * w * (2.0 * omega + w) * cos(theta)**2.0 + &
          0.25 * K**2.0 * cos(theta)**(2.0*R) * ((R+1.0)*cos(theta)**2.0 + 2.0*R**2.0 - R - 2.0 - 2.0*(R*cos(theta))**2.0)

   end function AA

   
   real function BB(theta)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! B, used in height field computation for Rossby-Haurwitz wave
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      real (kind=RKIND), intent(in) :: theta

      BB = (2.0*(omega + w)*K / ((R+1.0)*(R+2.0))) * cos(theta)**R * ((R**2.0 + 2.0*R + 2.0) - ((R+1.0)*cos(theta))**2.0)

   end function BB


   real function CC(theta)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C, used in height field computation for Rossby-Haurwitz wave
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      real (kind=RKIND), intent(in) :: theta

      CC = 0.25 * K**2.0 * cos(theta)**(2.0*R) * ((R+1.0)*cos(theta)**2.0 - R - 2.0)

   end function CC

end module test_cases
