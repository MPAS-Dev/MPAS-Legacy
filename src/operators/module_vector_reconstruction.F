module vector_reconstruction

   use grid_types
   use configure
   use constants

   implicit none

   public :: init_reconstruct, reconstruct

   contains

   subroutine init_reconstruct(grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Purpose: pre-compute coefficients used by the reconstruct() routine
   !
   ! Input: grid meta data
   !
   ! Output: grid % coeffs_reconstruct - coefficients used to reconstruct 
   !                                     velocity vectors at cell centers 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (grid_meta), intent(inout) :: grid 

      ! temporary arrays needed in the (to be constructed) init procedure
      integer :: nCells, nEdges, nVertLevels, nCellsSolve
      integer, dimension(:,:), pointer :: cellsOnCell, edgesOnCell, cellsOnEdge
      integer, dimension(:), pointer :: nEdgesOnCell
      integer :: iEdge, iCell, k, cell1, cell2, EdgeMax, j, i, npoints, matrixSize
      integer :: lwork, info
      integer, allocatable, dimension(:) :: pivotingIndices
      real (kind=RKIND), dimension(:), pointer :: dcEdge, xCell, yCell, zCell
      real (kind=RKIND) :: r, rbfValue, v, X1(3), X2(3), alpha, rHat(3), normalDotRHat
      real (kind=RKIND) :: xPlane, yPlane, xNormalPlane, yNormalPlane, xHatPlane(3), yHatPlane(3)
      real (kind=RKIND), allocatable, dimension(:,:,:) :: xLoc

      real (kind=RKIND), dimension(:,:,:), pointer :: coeffs_reconstruct
      real (kind=RKIND), allocatable, dimension(:,:) :: mwork
      real (kind=RKIND), dimension(:,:), pointer :: matrix, invMatrix
      real (kind=RKIND), dimension(:,:), pointer :: normals

      !========================================================
      ! arrays filled and saved during init procedure
      !========================================================
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

      allocate(normals(3,EdgeMax))
       
      ! init arrays
      coeffs_reconstruct = 0.0
      normals = 0

      ! loop over all cells to be solved on this block
      do iCell=1,nCellsSolve

         ! fill normal vector and xLoc arrays
         ! X1 is the location of the reconstruction, X2 are the neighboring centers, 
         !  xLoc is the edge positions
         cell1 = iCell
         X1(1) = xCell(cell1)
         X1(2) = yCell(cell1)
         X1(3) = zCell(cell1)

         rHat = X1
         call unit_vector_in_R3(rHat)

         do j=1,nEdgesOnCell(iCell)
           iEdge = edgesOnCell(j,iCell)
           if (iCell == cellsOnEdge(1,iEdge)) then
               cell2 = cellsOnEdge(2,iEdge)
               X2(1) = xCell(cell2)
               X2(2) = yCell(cell2)
               X2(3) = zCell(cell2)
               normals(:,j) = X2(:) - X1(:)
           else
               cell2 = cellsOnEdge(1,iEdge)
               X2(1) = xCell(cell2)
               X2(2) = yCell(cell2)
               X2(3) = zCell(cell2)
               normals(:,j) = X1(:) - X2(:)
           endif

           call unit_vector_in_R3(normals(:,j))

           xLoc(:,j,iCell) = 0.5*(X2(:) + X1(:))
         enddo

         npoints = nEdgesOnCell(iCell)   ! only loop over the number of edges for this cell
         matrixSize = npoints+2 ! room for 2 vector components for constant flow
                                !  basis functions
         allocate(matrix(matrixSize,matrixSize))
         matrix = 0.0
         alpha = 0.0
         do i=1,npoints
           call get_distance(xLoc(:,i,iCell),X1(:),r)
           alpha = alpha + r
         enddo
         alpha = alpha / npoints
         do j=1,npoints
           do i=1,npoints
              call get_distance(xLoc(:,i,iCell),xLoc(:,j,iCell),r)
              r = r / alpha
              call evaluate_rbf(r,rbfValue)
              call get_dotproduct(normals(:,i),normals(:,j),v)
              matrix(i,j) = v*rbfValue
           enddo
         enddo

         ! add matrix entries to suppoert constant flow
         ! xHat and yHat are a local basis in the plane of the horizontal cell
         ! we arbitrarily choose xHat to point toward the first edge
         call get_dotproduct(normals(:,1),rHat,normalDotRHat)
         xHatPlane = normals(:,1) - normalDotRHat*rHat(:)
         call unit_vector_in_R3(xHatPlane)

         call cross_product_in_R3(rHat, xHatPlane, yHatPlane)
         call unit_vector_in_R3(yHatPlane) ! just to be sure...

         do j=1,npoints
           call get_dotproduct(normals(:,j),xHatPlane, xNormalPlane)
           call get_dotproduct(normals(:,j),yHatPlane, yNormalPlane)
           matrix(j,npoints+1) = xNormalPlane
           matrix(j,npoints+2) = yNormalPlane
           matrix(npoints+1,j) = matrix(j,npoints+1)
           matrix(npoints+2,j) = matrix(j,npoints+2)
         end do
       
 
         ! invert matrix
         allocate(invMatrix(matrixSize,matrixSize))
         allocate(pivotingIndices(matrixSize))
         invMatrix = 0.0
         pivotingIndices = 0
         call migs(matrix,matrixSize,invMatrix,pivotingIndices)

         ! compute the coefficients for reconstructing uX, uY, uZ at cell centers from 
         ! u_i normal to edges
         ! uX = sum_j(coeffs(1,j) * u_j) (similarly for Y,Z)
         ! coeffs(:,j) = sum_i(rbf_values(i) * normal(:,i) * matrix(i,j))
         do i=1,npoints
           ! compute value of RBF when evaluated between reconstruction location and edge locations
           call get_distance(xLoc(:,i,iCell), X1(:), r)
           r = r / alpha
           call evaluate_rbf(r,rbfValue)

           ! project the normals onto tangent plane of the cell
           ! normal = normal - (normal dot r_hat) r_hat
           ! rHat, the unit vector pointing from the domain center to the cell center
           call get_dotproduct(normals(:,i),rHat,normalDotRHat)
           normals(:,i) = normals(:,i) - normalDotRHat*rHat(:)

           do j=1,npoints
              coeffs_reconstruct(:,j,iCell) = coeffs_reconstruct(:,j,iCell) &
                + rbfValue * normals(:,i) * invMatrix(i,j)
           enddo
         enddo
         ! polynomial parts
         do j=1,npoints
            coeffs_reconstruct(:,j,iCell) = coeffs_reconstruct(:,j,iCell) &
              + invMatrix(npoints+1,j)*xHatPlane
            coeffs_reconstruct(:,j,iCell) = coeffs_reconstruct(:,j,iCell) &
              + invMatrix(npoints+2,j)*yHatPlane
         enddo

         deallocate(matrix)
         deallocate(invMatrix)
         deallocate(pivotingIndices)

      enddo   ! iCell

      deallocate(xLoc)
      deallocate(normals)

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
      integer :: nCellsSolve
      integer, dimension(:,:), pointer :: edgesOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      integer :: iCell,iEdge, i
      real (kind=RKIND), dimension(:,:), pointer :: u
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructX, uReconstructY, uReconstructZ

      real (kind=RKIND), dimension(:,:,:), pointer :: coeffs_reconstruct

      ! stored arrays used during compute procedure
      coeffs_reconstruct => grid % coeffs_reconstruct % array

      ! temporary variables
      edgesOnCell => grid % edgesOnCell % array
      nEdgesOnCell=> grid % nEdgesOnCell % array
      nCellsSolve = grid % nCellsSolve
      u => s % u % array
      uReconstructX => s % uReconstructX % array
      uReconstructY => s % uReconstructY % array
      uReconstructZ => s % uReconstructZ % array

      ! init the intent(out)
      uReconstructX = 0.0
      uReconstructY = 0.0
      uReconstructZ = 0.0

      ! loop over cell centers
      do iCell=1,nCellsSolve
        ! a more efficient reconstruction where rbf_values*matrix_reconstruct has been precomputed
        ! in coeffs_reconstruct
        do i=1,nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i,iCell)
          uReconstructX(:,iCell) = uReconstructX(:,iCell) &
            + coeffs_reconstruct(1,i,iCell) * u(:,iEdge)
          uReconstructY(:,iCell) = uReconstructY(:,iCell) &
            + coeffs_reconstruct(2,i,iCell) * u(:,iEdge)
          uReconstructZ(:,iCell) = uReconstructZ(:,iCell) &
            + coeffs_reconstruct(3,i,iCell) * u(:,iEdge)
        enddo
      enddo   ! iCell

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


   subroutine evaluate_rbf(xin,xout)
     real(kind=RKIND), intent(in) ::  xin
     real(kind=RKIND), intent(out) :: xout

     ! Gaussian
     ! xout = exp(-r**2)

     ! multiquadrics
       xout = 1.0 / sqrt(1.0**2 + xin**2)

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

! Updated 10/24/2001.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!   Program 4.4   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! Please Note:                                                          !
!                                                                       !
! (1) This computer program is written by Tao Pang in conjunction with  !
!     his book, "An Introduction to Computational Physics," published   !
!     by Cambridge University Press in 1997.                            !
!                                                                       !
! (2) No warranties, express or implied, are made for this program.     !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE MIGS (A,N,X,INDX)
!
! Subroutine to invert matrix A(N,N) with the inverse stored
! in X(N,N) in the output.  Copyright (c) Tao Pang 2001.
!
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  INTEGER :: I,J,K
  INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
  REAL (kind=RKIND), INTENT (INOUT), DIMENSION (N,N):: A
  REAL (kind=RKIND), INTENT (OUT), DIMENSION (N,N):: X
  REAL (kind=RKIND), DIMENSION (N,N) :: B
!
  DO I = 1, N
    DO J = 1, N
      B(I,J) = 0.0
    END DO
  END DO
  DO I = 1, N
    B(I,I) = 1.0
  END DO
!
  CALL ELGS (A,N,INDX)
!
  DO I = 1, N-1
    DO J = I+1, N
      DO K = 1, N
        B(INDX(J),K) = B(INDX(J),K)-A(INDX(J),I)*B(INDX(I),K)
      END DO
    END DO
  END DO
!
  DO I = 1, N
    X(N,I) = B(INDX(N),I)/A(INDX(N),N)
    DO J = N-1, 1, -1
      X(J,I) = B(INDX(J),I)
      DO K = J+1, N
        X(J,I) = X(J,I)-A(INDX(J),K)*X(K,I)
      END DO
      X(J,I) =  X(J,I)/A(INDX(J),J)
    END DO
  END DO
END SUBROUTINE MIGS


SUBROUTINE ELGS (A,N,INDX)
!
! Subroutine to perform the partial-pivoting Gaussian elimination.
! A(N,N) is the original matrix in the input and transformed matrix
! plus the pivoting element ratios below the diagonal in the output.
! INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.
!
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  INTEGER :: I,J,K,ITMP
  INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
  REAL (kind=RKIND) :: C1,PI,PI1,PJ
  REAL (kind=RKIND), INTENT (INOUT), DIMENSION (N,N) :: A
  REAL (kind=RKIND), DIMENSION (N) :: C
!
! Initialize the index
!
  DO I = 1, N
    INDX(I) = I
  END DO
!
! Find the rescaling factors, one from each row
!
  DO I = 1, N
    C1= 0.0
    DO J = 1, N
      !C1 = AMAX1(C1,ABS(A(I,J)))
      C1 = MAX(C1,ABS(A(I,J)))
    END DO
    C(I) = C1
  END DO
!
! Search the pivoting (largest) element from each column
!
  DO J = 1, N-1
    PI1 = 0.0
    DO I = J, N
      PI = ABS(A(INDX(I),J))/C(INDX(I))
      IF (PI.GT.PI1) THEN
        PI1 = PI
        K   = I
      ENDIF
    END DO
!
! Interchange the rows via INDX(N) to record pivoting order
!
    ITMP    = INDX(J)
    INDX(J) = INDX(K)
    INDX(K) = ITMP
    DO I = J+1, N
      PJ  = A(INDX(I),J)/A(INDX(J),J)
!
! Record pivoting ratios below the diagonal
!
      A(INDX(I),J) = PJ
!
! Modify other elements accordingly
!
      DO K = J+1, N
        A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
      END DO
    END DO
  END DO
!
END SUBROUTINE ELGS

end module vector_reconstruction
