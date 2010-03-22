module vector_reconstruction

   use grid_types
   use configure
   use constants

   implicit none

   public :: init_reconstruct, reconstruct

   contains

   subroutine init_reconstruct(grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Purpose: reconstruct vector field at cell centers based on radial basis functions 
   !
   ! Input: grid meta data and vector component data residing at cell edges
   !
   ! Output: reconstructed vector field (measured in X,Y,Z) located at cell centers
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (grid_meta), intent(inout) :: grid 

      ! temporary arrays needed in the (to be constructed) init procedure
      integer :: nCells, nEdges, nVertLevels, nCellsSolve
      integer, dimension(:,:), pointer :: cellsOnCell, edgesOnCell, cellsOnEdge
      integer, dimension(:), pointer :: nEdgesOnCell
      integer :: iEdge, iCell, k, cell1, cell2, EdgeMax, j, i, npoints
      integer :: lwork, info
      integer, allocatable, dimension(:) :: ipvt
      real (kind=RKIND), allocatable, dimension(:) :: work
      real (kind=RKIND), dimension(:), pointer :: dcEdge, xCell, yCell, zCell
      real (kind=RKIND) :: r, t, v, X1(3), X2(3), alpha
      real (kind=RKIND), allocatable, dimension(:,:,:) :: xLoc

      real (kind=RKIND), dimension(:,:,:), pointer :: matrix_reconstruct
      real (kind=RKIND), dimension(:,:,:), pointer :: normal
      real (kind=RKIND), dimension(:,:), pointer :: rbf_value
      real (kind=RKIND), dimension(:,:,:), pointer :: coeffs_reconstruct
      real (kind=RKIND), allocatable, dimension(:,:) :: mwork

      !========================================================
      ! arrays filled and saved during init procedure
      !========================================================
      matrix_reconstruct => grid % matrix_reconstruct % array
      normal => grid % normal % array
      rbf_value => grid % rbf_value % array
      coeffs_reconstruct => grid % coeffs_reconstruct % array

      !========================================================
      ! temporary variables needed for init procedure
      !========================================================
      xCell       => grid % xCell % array
      yCell       => grid % yCell % array
      zCell       => grid % zCell % array
      cellsOnCell => grid % cellsOnCell % array
      cellsOnEdge => grid % cellsOnEdge % array
      edgesOnCell => grid % edgesOnCell % array
      nEdgesOnCell=> grid % nEdgesOnCell % array
      dcEdge      => grid % dcEdge % array
      nCells      = grid % nCells
      nCellsSolve = grid % nCellsSolve
      nEdges      = grid % nEdges
      nVertLevels = grid % nVertLevels

      ! allocate arrays
      EdgeMax = maxval(nEdgesOnCell)
      allocate(xLoc(3,EdgeMax,nCells))
      allocate(work(EdgeMax*EdgeMax))
      allocate(ipvt(EdgeMax))
       
       ! init arrays
       matrix_reconstruct = 0.0
       normal = 0.0
       rbf_value = 0.0

       ! loop over all cells to be solved on this block
       do iCell=1,nCellsSolve

         npoints = nEdgesOnCell(iCell)   ! only loop over the number of edges for this cell

         ! fill normal vector and xLoc arrays
         ! X1 is the location of the reconstruction, X2 are the neighboring centers, xLoc is the edge positions
         cell1 = iCell
         X1(1) = xCell(cell1)
         X1(2) = yCell(cell1)
         X1(3) = zCell(cell1)
         do j=1,npoints
           cell2 = cellsOnCell(j,iCell)
           iEdge = edgesOnCell(j,iCell)
           X2(1) = xCell(cell2)
           X2(2) = yCell(cell2)
           X2(3) = zCell(cell2)
           if(cell2.gt.cell1) then
               normal(:,j,iCell) = X2(:) - X1(:)
           else
               normal(:,j,iCell) = X1(:) - X2(:)
           endif
           call unit_vector_in_R3(normal(:,j,iCell))
           xLoc(:,j,iCell) = 0.5*(X2(:) + X1(:))
         enddo

         alpha = 0.0
         do i=1,npoints
           call get_distance(xLoc(:,i,iCell),X1(:),r)
           alpha = alpha + r
         enddo
         alpha = 4*alpha / npoints
         do j=1,npoints
           do i=1,npoints
              call get_distance(xLoc(:,i,iCell),xLoc(:,j,iCell),r)
              r = r / alpha
              call evaluate_rbf(r,t)
              call get_dotproduct(normal(:,i,iCell),normal(:,j,iCell),v)
              matrix_reconstruct(i,j,iCell) = v*t
           enddo
         enddo
 
         ! save value of RBF when evaluated between reconstruction location and edge locations
         do j=1,npoints
           call get_distance(xLoc(:,j,iCell), X1(:), r)
           r = r / alpha
           call evaluate_rbf(r,t)
           rbf_value(j,iCell) = t
         enddo

         ! invert matrix_reconstruct using lapack routines
         ! xsad 10-02-09:
         !!!!! the lapack stuff isn't working on coyote in parallel !!!!!!!
         !allocate(mwork(npoints,npoints))
         !lwork = npoints*npoints
         !mwork(:,:) = matrix_reconstruct(1:npoints,1:npoints,iCell)
         !call dgetrf(npoints, npoints, mwork, npoints, ipvt, info)
         !if(info.ne.0) then
          !write(6,*) info, 'dgetrf'
          !stop
         !endif
         !call dgetri(npoints, mwork, npoints, ipvt, work, lwork, info)
         !matrix_reconstruct(1:npoints,1:npoints,iCell) = mwork(:,:)
         !if(info.ne.0) then
          !write(6,*) info, 'dgetri'
          !stop
         !endif
          !deallocate(mwork)
          ! xsad 10-02-09 end

         ! compute the coefficients for reconstructing uX, uY, uZ at cell centers from u_i normal to edges
         ! uX = sum_j(coeffs(1,j) * u_j) (similarly for Y,Z)
         ! coeffs(:,j) = sum_i(rbf_values(i) * normal(:,i) * matrix(i,j))
         do j=1,npoints
           coeffs_reconstruct(:,j,iCell) = 0.0
           do i=1,npoints
              coeffs_reconstruct(:,j,iCell) = coeffs_reconstruct(:,j,iCell) + rbf_value(i,iCell) * normal(:,i,iCell) &
                * matrix_reconstruct(i,j,iCell)
           enddo
         enddo

      enddo   ! iCell

      deallocate(ipvt)
      deallocate(work)
      deallocate(xLoc)

      end subroutine init_reconstruct


      subroutine reconstruct(s, grid)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Purpose: reconstruct vector field at cell centers based on radial basis functions
      !
      ! Input: grid meta data and vector component data residing at cell edges
      !
      ! Output: reconstructed vector field (measured in X,Y,Z) located at cell centers
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_state), intent(inout) :: s 
      type (grid_meta), intent(in) :: grid

      !   temporary arrays needed in the compute procedure
      integer :: nCells, nEdges, nVertLevels, nCellsSolve
      integer, dimension(:,:), pointer :: edgesOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      integer :: iCell,iEdge, i, j, k, npoints, EdgeMax
      real (kind=RKIND) :: r, t, p!, t1(3), t2(3), t3(3), east(3)
      real (kind=RKIND), dimension(:,:), pointer :: u
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructX, uReconstructY, uReconstructZ
      real (kind=RKIND), allocatable, dimension(:,:) :: rhs, c

      real (kind=RKIND), dimension(:,:,:), pointer :: matrix_reconstruct
      real (kind=RKIND), dimension(:,:,:), pointer :: normal
      real (kind=RKIND), dimension(:,:), pointer :: rbf_value
      real (kind=RKIND), dimension(:,:,:), pointer :: coeffs_reconstruct

      ! stored arrays used during compute procedure
      matrix_reconstruct => grid % matrix_reconstruct % array
      normal => grid % normal % array
      rbf_value => grid % rbf_value % array
      coeffs_reconstruct => grid % coeffs_reconstruct % array

      ! temporary variables
      edgesOnCell => grid % edgesOnCell % array
      nEdgesOnCell=> grid % nEdgesOnCell % array
      nCells      = grid % nCells
      nCellsSolve = grid % nCellsSolve
      nEdges      = grid % nEdges
      nVertLevels = grid % nVertLevels
      xCell => grid % xCell % array
      yCell => grid % yCell % array
      zCell => grid % zCell % array

      u => s % u % array
      uReconstructX => s % uReconstructX % array
      uReconstructY => s % uReconstructY % array
      uReconstructZ => s % uReconstructZ % array

      ! allocate space for temporary arrays
      EdgeMax = maxval(nEdgesOnCell)
!      allocate(c(nVertLevels,EdgeMax))
      allocate(rhs(nVertLevels,EdgeMax))

      ! init the intent(out)
      uReconstructX = 0.0
      uReconstructY = 0.0
      uReconstructZ = 0.0

      ! loop over cell centers
      do iCell=1,nCellsSolve

        npoints = nEdgesOnCell(iCell)

        ! fill the RHS of the matrix equation
        ! for testing purposes, fill rhs with test function
        ! test function assumes uniform flow to the east
        rhs = 0.0
        do j=1,npoints
          iEdge = edgesOnCell(j,iCell)
          do k=1,nVertLevels
            rhs(k,j) = u(k,iEdge)
          enddo
        enddo
  
        ! compute c by multiplying matrix_reconstruct * rhs  (Eq 8 in Tex document)
!        c = 0.0
!        do i=1,npoints
!         do j=1,npoints
!           do k=1,nVertLevels
!             c(k,i) = c(k,i) + matrix_reconstruct(i,j,iCell)*rhs(k,j)
!           enddo
!         enddo
!        enddo
  
        ! reconstruct the velocity field at point X1 (Eq 6 in Tex document)
!        do i=1,npoints
!          do k=1,nVertLevels
!            uReconstructX(k,iCell) = uReconstructX(k,iCell) + c(k,i)*rbf_value(i,iCell)*normal(1,i,iCell)
!            uReconstructY(k,iCell) = uReconstructY(k,iCell) + c(k,i)*rbf_value(i,iCell)*normal(2,i,iCell)
!            uReconstructZ(k,iCell) = uReconstructZ(k,iCell) + c(k,i)*rbf_value(i,iCell)*normal(3,i,iCell)
!          enddo
!        enddo
!
        ! a more efficient reconstruction where rbf_values*matrix_reconstruct has been precomputed
        ! in coeffs_reconstruct
        do i=1,npoints
          uReconstructX(:,iCell) = uReconstructX(:,iCell) + coeffs_reconstruct(1,i,iCell) * rhs(:,i)
          uReconstructY(:,iCell) = uReconstructY(:,iCell) + coeffs_reconstruct(2,i,iCell) * rhs(:,i)
          uReconstructZ(:,iCell) = uReconstructZ(:,iCell) + coeffs_reconstruct(3,i,iCell) * rhs(:,i)
        enddo

        
      enddo   ! iCell

      ! deallocate
      deallocate(rhs)
!      deallocate(c)
 
   end subroutine reconstruct

   subroutine get_distance(x1,x2,xout)
     implicit none
     real(kind=RKIND), intent(in) :: x1(3), x2(3)
     real(kind=RKIND), intent(out) :: xout
     xout = sqrt( (x1(1)-x2(1))**2 + (x1(2)-x2(2))**2 + (x1(3)-x2(3))**2 )
   end subroutine get_distance
 
   subroutine get_dotproduct(x1,x2,xout)
     implicit none
     real(kind=RKIND), intent(in) :: x1(3), x2(3)
     real(kind=RKIND), intent(out) :: xout
     xout = x1(1)*x2(1) + x1(2)*x2(2) + x1(3)*x2(3)
   end subroutine get_dotproduct
 
 
   subroutine unit_vector_in_R3(xin)
     implicit none
     real (kind=RKIND), intent(inout) :: xin(3)
     real (kind=RKIND) :: mag
     mag = sqrt(xin(1)**2+xin(2)**2+xin(3)**2)
     xin(:) = xin(:) / mag
   end subroutine unit_vector_in_R3


   subroutine evaluate_rbf(r,xout)
     real(kind=RKIND), intent(in) ::  r
     real(kind=RKIND), intent(out) :: xout

     ! Gaussian
       xout = exp(-r**2)

     ! multiquadrics
     ! xout = 1.0 / sqrt(1.0**2 + r**2)

     ! other
     ! xout = 1.0 / (1.0**2 + r**2)

   end subroutine evaluate_rbf

!======================================================================
! BEGINNING OF CROSS_PRODUCT_IN_R3
!======================================================================
        subroutine cross_product_in_R3(p_1,p_2,p_out)

!-----------------------------------------------------------------------
! PURPOSE: compute p_1 cross p_2 and place in p_out
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! intent(in)
!-----------------------------------------------------------------------
        real (kind=RKIND), intent(in) ::                            &
                        p_1 (3),                                      &
                        p_2 (3)

!-----------------------------------------------------------------------
! intent(out)
!-----------------------------------------------------------------------
        real (kind=RKIND), intent(out) ::                           &
                        p_out (3)

        p_out(1) = p_1(2)*p_2(3)-p_1(3)*p_2(2)
        p_out(2) = p_1(3)*p_2(1)-p_1(1)*p_2(3)
        p_out(3) = p_1(1)*p_2(2)-p_1(2)*p_2(1)

        end subroutine cross_product_in_R3
!======================================================================
! END OF CROSS_PRODUCT_IN_R3
!======================================================================


end module vector_reconstruction
