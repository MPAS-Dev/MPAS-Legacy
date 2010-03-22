module test_cases

   use grid_types
   use configure
   use constants


   contains


   subroutine setup_sw_test_case(domain)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Configure grid metadata and model state for the shallow water test case 
   !   specified in the namelist
   !
   ! Output: block - a subset (not necessarily proper) of the model domain to be
   !                 initialized
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain

      integer :: i, iCell, iEdge, iVtx, iLevel
      type (block_type), pointer :: block_ptr
      real (kind=RKIND), dimension(:,:), pointer :: h, u, u_src, rho
      real (kind=RKIND) :: delta_rho
      integer :: nCells, nEdges, nVertices, nVertLevels

      if (config_test_case == 0) then
         write(0,*) 'Using initial conditions supplied in input file'

      else if (config_test_case == 1) then
         write(0,*) ' Setting up shallow water test case 1:'
         write(0,*) ' Advection of Cosine Bell over the Pole'

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call sw_test_case_1(block_ptr % mesh, block_ptr % time_levs(1) % state)
            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 2) then
         write(0,*) ' Setup shallow water test case 2: '// &
           'Global Steady State Nonlinear Zonal Geostrophic Flow'

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call sw_test_case_2(block_ptr % mesh, block_ptr % time_levs(1) % state)
            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 5) then
         write(0,*) ' Setup shallow water test case 5:'// &
           ' Zonal Flow over an Isolated Mountain'

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call sw_test_case_5(block_ptr % mesh, block_ptr % time_levs(1) % state)
            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 6) then
         write(0,*) ' Set up shallow water test case 6:'
         write(0,*) ' Rossby-Haurwitz Wave'

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call sw_test_case_6(block_ptr % mesh, block_ptr % time_levs(1) % state)
            block_ptr => block_ptr % next
         end do

      else
         write(0,*) &
           'Only test case 1, 2, 5, and 6 are currently supported.'
         stop
      end if

      ! mrp 100121:
      !
      ! Initialize u_src, rho, alpha
      ! This is a temporary fix until everything is specified correctly 
      ! in an initial conditions file.
      !
      block_ptr => domain % blocklist
      do while (associated(block_ptr))
         h          => block_ptr % time_levs(1) % state % h % array
         u          => block_ptr % time_levs(1) % state % u % array
         rho        => block_ptr % time_levs(1) % state % rho % array

         u_src      => block_ptr % mesh % u_src % array

         nCells      = block_ptr % mesh % nCells
         nEdges      = block_ptr % mesh % nEdges
         nVertices   = block_ptr % mesh % nVertices
         nVertLevels = block_ptr % mesh % nVertLevels


         ! Momentum forcing u_src
         ! for shallow water test cases:
         ! u_src = 0.0
         ! for rectangular basin:
         do iEdge=1,nEdges
            u_src(1:nVertLevels, iEdge) = u_src(1, iEdge) / nVertLevels
         end do

         ! define the density for multiple layers
         delta_rho=0.0
         do iLevel = 1,nVertLevels
           rho(iLevel,1) = delta_rho*(iLevel-1)
         enddo
         rho(:,1) = rho(:,1) - sum(rho(:,1))/nVertLevels + 1000.0
         do iLevel = 1,nVertLevels
           rho(iLevel,:) = rho(iLevel,1)
         enddo

         do i=2,nTimeLevs
            call copy_state(block_ptr % time_levs(1) % state, &
                            block_ptr % time_levs(i) % state)
         end do

         ! print some diagnostics
         print '(10a)', 'ilevel',&
            '  rho      ',&
            '  min u       max u     ',&
            '  min h       max h     ',&
            '  min u_src   max u_src '
         do iLevel = 1,nVertLevels
            print '(i5,20es12.4)', ilevel, rho(ilevel,1), &
              minval(u(iLevel,:)), maxval(u(iLevel,:)), &
              minval(h(iLevel,:)), maxval(h(iLevel,:)), &
              minval(u_src(iLevel,:)), maxval(u_src(iLevel,:))
         enddo

         block_ptr => block_ptr % next
      end do
      ! mrp 100121 end

   end subroutine setup_sw_test_case


   subroutine sw_test_case_1(grid, state)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup shallow water test case 1: Advection of Cosine Bell over the Pole
   !
   ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
   !            Approximations to the Shallow Water Equations in Spherical 
   !            Geometry" J. of Comp. Phys., 102, pp. 211--224
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_meta), intent(inout) :: grid
      type (grid_state), intent(inout) :: state

      real (kind=RKIND), parameter :: u0 = 2.0 * pii * a / (12.0 * 86400.0)
      real (kind=RKIND), parameter :: h0 = 1000.0
      real (kind=RKIND), parameter :: theta_c = 0.0
      real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
      real (kind=RKIND), parameter :: alpha = pii/4.0

      integer :: iCell, iEdge, iVtx
      real (kind=RKIND) :: r, u, v
      real (kind=RKIND), allocatable, dimension(:) :: psiVertex

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

      !
      ! Initialize wind field
      !
      allocate(psiVertex(grid % nVertices))
      do iVtx=1,grid % nVertices
         psiVertex(iVtx) = -a * u0 * ( &
                                       sin(grid%latVertex%array(iVtx)) * cos(alpha) - &
                                       cos(grid%lonVertex%array(iVtx)) * cos(grid%latVertex%array(iVtx)) * sin(alpha) &
                                     )
      end do
      do iEdge=1,grid % nEdges
         state % u % array(1,iEdge) = -1.0 * ( &
                                               psiVertex(grid%verticesOnEdge%array(2,iEdge)) - &
                                               psiVertex(grid%verticesOnEdge%array(1,iEdge)) &
                                             ) / grid%dvEdge%array(iEdge)
      end do
      deallocate(psiVertex)

      !
      ! Initialize cosine bell at (theta_c, lambda_c)
      !
      do iCell=1,grid % nCells
         r = sphere_distance(theta_c, lambda_c, grid % latCell % array(iCell), grid % lonCell % array(iCell), a) 
         if (r < a/3.0) then
            state % h % array(1,iCell) = (h0 / 2.0) * (1.0 + cos(pii*r*3.0/a))
         else
            state % h % array(1,iCell) = 0.0
         end if
      end do

   end subroutine sw_test_case_1


   subroutine sw_test_case_2(grid, state)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup shallow water test case 2: Global Steady State Nonlinear Zonal 
   !                                  Geostrophic Flow
   !
   ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
   !            Approximations to the Shallow Water Equations in Spherical 
   !            Geometry" J. of Comp. Phys., 102, pp. 211--224
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_meta), intent(inout) :: grid
      type (grid_state), intent(inout) :: state

      real (kind=RKIND), parameter :: u0 = 2.0 * pii * a / (12.0 * 86400.0)
      real (kind=RKIND), parameter :: gh0 = 29400.0
      real (kind=RKIND), parameter :: alpha = 0.0

      integer :: iCell, iEdge, iVtx
      real (kind=RKIND) :: u, v
      real (kind=RKIND), allocatable, dimension(:) :: psiVertex


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
      

      !
      ! Initialize wind field
      !
      allocate(psiVertex(grid % nVertices))
      do iVtx=1,grid % nVertices
         psiVertex(iVtx) = -a * u0 * ( &
                                       sin(grid%latVertex%array(iVtx)) * cos(alpha) - &
                                       cos(grid%lonVertex%array(iVtx)) * cos(grid%latVertex%array(iVtx)) * sin(alpha) &
                                     )
      end do
      do iEdge=1,grid % nEdges
         state % u % array(1,iEdge) = -1.0 * ( &
                                               psiVertex(grid%verticesOnEdge%array(2,iEdge)) - &
                                               psiVertex(grid%verticesOnEdge%array(1,iEdge)) &
                                             ) / grid%dvEdge%array(iEdge)
      end do
      deallocate(psiVertex)

      !
      ! Generate rotated Coriolis field
      !
      do iEdge=1,grid % nEdges
         grid % fEdge % array(iEdge) = 2.0 * omega * &
                                       ( -cos(grid%lonEdge%array(iEdge)) * cos(grid%latEdge%array(iEdge)) * sin(alpha) + &
                                         sin(grid%latEdge%array(iEdge)) * cos(alpha) &
                                       )
      end do
      do iVtx=1,grid % nVertices
         grid % fVertex % array(iVtx) = 2.0 * omega * &
                                         (-cos(grid%lonVertex%array(iVtx)) * cos(grid%latVertex%array(iVtx)) * sin(alpha) + &
                                          sin(grid%latVertex%array(iVtx)) * cos(alpha) &
                                         )
      end do

      !
      ! Initialize height field (actually, fluid thickness field)
      !
      do iCell=1,grid % nCells
         state % h % array(1,iCell) = (gh0 - (a * omega * u0 + 0.5 * u0**2.0) * &
                                             (-cos(grid%lonCell%array(iCell)) * cos(grid%latCell%array(iCell)) * sin(alpha) + &
                                              sin(grid%latCell%array(iCell)) * cos(alpha) &
                                             )**2.0 &
                                      ) / &
                                      gravity
      end do

   end subroutine sw_test_case_2


   subroutine sw_test_case_5(grid, state)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup shallow water test case 5: Zonal Flow over an Isolated Mountain
   !
   ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
   !            Approximations to the Shallow Water Equations in Spherical 
   !            Geometry" J. of Comp. Phys., 102, pp. 211--224
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_meta), intent(inout) :: grid
      type (grid_state), intent(inout) :: state

      real (kind=RKIND), parameter :: u0 = 20.
      real (kind=RKIND), parameter :: gh0 = 5960.0*gravity
!      real (kind=RKIND), parameter :: hs0 = 2000. original
      real (kind=RKIND), parameter :: hs0 = 250.  !mrp 100204
      real (kind=RKIND), parameter :: theta_c = pii/6.0
      real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
      real (kind=RKIND), parameter :: rr = pii/9.0
      real (kind=RKIND), parameter :: alpha = 0.0

      integer :: iCell, iEdge, iVtx
      real (kind=RKIND) :: r, u, v
      real (kind=RKIND), allocatable, dimension(:) :: psiVertex


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

      !
      ! Initialize wind field
      !
      allocate(psiVertex(grid % nVertices))
      do iVtx=1,grid % nVertices
         psiVertex(iVtx) = -a * u0 * ( &
                                       sin(grid%latVertex%array(iVtx)) * cos(alpha) - &
                                       cos(grid%lonVertex%array(iVtx)) * cos(grid%latVertex%array(iVtx)) * sin(alpha) &
                                     )
      end do
      do iEdge=1,grid % nEdges
         state % u % array(1,iEdge) = -1.0 * ( &
                                               psiVertex(grid%verticesOnEdge%array(2,iEdge)) - &
                                               psiVertex(grid%verticesOnEdge%array(1,iEdge)) &
                                             ) / grid%dvEdge%array(iEdge)
      end do
      deallocate(psiVertex)

      !
      ! Generate rotated Coriolis field
      !
      do iEdge=1,grid % nEdges
         grid % fEdge % array(iEdge) = 2.0 * omega * &
                                        (-cos(grid%lonEdge%array(iEdge)) * cos(grid%latEdge%array(iEdge)) * sin(alpha) + &
                                         sin(grid%latEdge%array(iEdge)) * cos(alpha) &
                                        )
      end do
      do iVtx=1,grid % nVertices
         grid % fVertex % array(iVtx) = 2.0 * omega * &
                                         (-cos(grid%lonVertex%array(iVtx)) * cos(grid%latVertex%array(iVtx)) * sin(alpha) + &
                                          sin(grid%latVertex%array(iVtx)) * cos(alpha) &
                                         )
      end do

      !
      ! Initialize mountain
      !
      do iCell=1,grid % nCells
         if (grid % lonCell % array(iCell) < 0.0) grid % lonCell % array(iCell) = grid % lonCell % array(iCell) + 2.0 * pii
         r = sqrt(min(rr**2.0, (grid % lonCell % array(iCell) - lambda_c)**2.0 + (grid % latCell % array(iCell) - theta_c)**2.0))
         grid % h_s % array(iCell) = hs0 * (1.0 - r/rr)
      end do
! output about mountain
print *, 'h_s',minval(grid % h_s % array),sum(grid % h_s % array)/grid % nCells, maxval(grid % h_s % array)

      !
      ! Initialize tracer fields
      !
      do iCell=1,grid % nCells
         r = sqrt(min(rr**2.0, (grid % lonCell % array(iCell) - lambda_c)**2.0 + (grid % latCell % array(iCell) - theta_c)**2.0))
         state % tracers % array(1,1,iCell) = 1.0 - r/rr
      end do
      do iCell=1,grid % nCells
         r = sqrt(min(rr**2.0, (grid % lonCell % array(iCell) - lambda_c)**2.0 + &
                      (grid % latCell % array(iCell) - theta_c - pii/6.0)**2.0 &
                     ) &
                 )
         state % tracers % array(2,1,iCell) = 1.0 - r/rr
      end do

      !
      ! Initialize height field (actually, fluid thickness field)
      !
      do iCell=1,grid % nCells
         state % h % array(1,iCell) = (gh0 - (a * omega * u0 + 0.5 * u0**2.0) * &
                                         (-cos(grid%lonCell%array(iCell)) * cos(grid%latCell%array(iCell)) * sin(alpha) + &
                                          sin(grid%latCell%array(iCell)) * cos(alpha) &
                                         )**2.0 &
                                      ) / &
                                      gravity
         state % h % array(1,iCell) = state % h % array(1,iCell) - grid % h_s % array(iCell)
      end do

   end subroutine sw_test_case_5


   subroutine sw_test_case_6(grid, state)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup shallow water test case 6: Rossby-Haurwitz Wave
   !
   ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
   !            Approximations to the Shallow Water Equations in Spherical 
   !            Geometry" J. of Comp. Phys., 102, pp. 211--224
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_meta), intent(inout) :: grid
      type (grid_state), intent(inout) :: state

      real (kind=RKIND), parameter :: h0 = 8000.0
      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      integer :: iCell, iEdge, iVtx
      real (kind=RKIND) :: u, v
      real (kind=RKIND), allocatable, dimension(:) :: psiVertex


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

      !
      ! Initialize wind field
      !
      allocate(psiVertex(grid % nVertices))
      do iVtx=1,grid % nVertices
         psiVertex(iVtx) = -a * a * w * sin(grid%latVertex%array(iVtx)) + &
                            a *a * K * (cos(grid%latVertex%array(iVtx))**R) * &
                            sin(grid%latVertex%array(iVtx)) * cos(R * grid%lonVertex%array(iVtx))
      end do
      do iEdge=1,grid % nEdges
         state % u % array(1,iEdge) = -1.0 * ( &
                                               psiVertex(grid%verticesOnEdge%array(2,iEdge)) - &
                                               psiVertex(grid%verticesOnEdge%array(1,iEdge)) &
                                             ) / grid%dvEdge%array(iEdge)
      end do
      deallocate(psiVertex)

      !
      ! Initialize height field (actually, fluid thickness field)
      !
      do iCell=1,grid % nCells
         state % h % array(1,iCell) = (gravity * h0 + a*a*AA(grid%latCell%array(iCell)) + &
                                                      a*a*BB(grid%latCell%array(iCell)) * cos(R*grid%lonCell%array(iCell)) + &
                                                      a*a*CC(grid%latCell%array(iCell)) * cos(2.0*R*grid%lonCell%array(iCell)) &
                                      ) / gravity
      end do

   end subroutine sw_test_case_6


   real function sphere_distance(lat1, lon1, lat2, lon2, radius)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on a
   !   sphere with given radius.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius

      real (kind=RKIND) :: arg1

      arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
                   cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
      sphere_distance = 2.*radius*asin(arg1)

   end function sphere_distance


   real function AA(theta)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! A, used in height field computation for Rossby-Haurwitz wave
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      real (kind=RKIND), intent(in) :: theta

      AA = 0.5 * w * (2.0 * omega + w) * cos(theta)**2.0 + &
          0.25 * K**2.0 * cos(theta)**(2.0*R) * ((R+1.0)*cos(theta)**2.0 + 2.0*R**2.0 - R - 2.0 - 2.0*R**2*cos(theta)**-2.0)

   end function AA

   
   real function BB(theta)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! B, used in height field computation for Rossby-Haurwitz wave
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      real (kind=RKIND), intent(in) :: theta

      BB = (2.0*(omega + w)*K / ((R+1.0)*(R+2.0))) * cos(theta)**R * ((R**2.0 + 2.0*R + 2.0) - ((R+1.0)*cos(theta))**2.0)

   end function BB


   real function CC(theta)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C, used in height field computation for Rossby-Haurwitz wave
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      real (kind=RKIND), intent(in) :: theta

      CC = 0.25 * K**2.0 * cos(theta)**(2.0*R) * ((R+1.0)*cos(theta)**2.0 - R - 2.0)

   end function CC

end module test_cases
