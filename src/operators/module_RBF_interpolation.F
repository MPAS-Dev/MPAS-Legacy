module RBF_interpolation
   use dmpar
   use grid_types

   implicit none
   private
   save

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Purpose: perform interpolation of scalar and vector functions in 2D
!   and 3D using Radial Basis Functions (RBFs).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

   ! Initialize the geometry that will be useful from interpolation
  public :: rbfInterp_initialize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Routines for perofrming interpolation in 2D (including Jacobian and Hessian)
  !  at locations that vary using a function that is fixed.  This is useful
  !  for finding the location on the the RBF reconstruction of a function
  !  (e.g., a height field) that minimizes the distance to a point in 3D space
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  public :: rbfInterp_loc_2D_sca_const_compCoeffs, &
    rbfInterp_loc_2D_sca_lin_compCoeffs, &
    rbfInterp_loc_2D_sca_const_evalWithDerivs, &
    rbfInterp_loc_2D_sca_lin_evalWithDerivs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Routines for computing scalar interpolaiton coefficients in 3D (coplanar points
  !  in 3D) with support for either constant or constant and linear basis
  !  functions in addition to RBFs.  In constrast to the two subroutines
  !  above, these coefficients are valid for computing the value of the
  !  interpolant at a fixe point in space using function values that may
  !  vary (e.g., in time) at each of the interpolation "source" points.
  !  The last 3 routines can be used to compute coefficients for imposing both Neumann
  !  and Dirichlet boundary conditions.
  ! Pseudocode for function reconstruction at the destinationPoint is as follows
  !  Dirichlet:  functionAtDestination = sum(functionAtSources*dirichletCoefficients)
  !  Neumann:    functionAtDestination = sum(functionOrDerivAtSources*neumannCoefficients)
  !    where functionOrDerivAtSource = functionAtSources where .not.(isInterface)
  !                                  = functionNormalDerivAtSources where isInterface
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  public :: rbfInterp_func_3D_sca_const_dir_compCoeffs, &
    rbfInterp_func_3DPlane_sca_lin_dir_compCoeffs, &
    rbfInterp_func_3D_sca_lin_dir_compCoeffs, &
    rbfInterp_func_3D_sca_const_dirNeu_compCoeffs, &
    rbfInterp_func_3DPlane_sca_lin_dirNeu_compCoeffs, &
    rbfInterp_func_3D_sca_lin_dirNeu_compCoeffs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Routines for computing vector interpolaiton coefficients in 3D (coplanar points
  !  in 3D) with support for only constant basis functions in addition to RBFs. 
  !  (Linear basis functions can cause problems: a linear vortex flow u = y xHat - x yHat
  !  cannot be reconstructed from vectors that are only normal to the cell edges in 2D.
  !  Therefore, we don't support them). As with the scalar 3D subroutines
  !  above, these coefficients are valid for computing the value of the
  !  interpolant at a fixe point in space using function values that may
  !  vary (e.g., in time) at each of the interpolation "source" points.
  !  The user supplies to these routines a set of sourcePoints and unitVectors
  !  as well as a destinationPoint and, for the last 2 routines, flags for
  !  which points are tangent to the interface and which of the supplied unitVectors
  !  is the normal at the corresponding point.
  !
  ! The first two routines compute coefficients that can be used to interpolate
  !  a vector function to the destination point given "function dot unitVector" values
  !  at each source point.  These routines are useful, for example, for reconstructing
  !  the full vector velocity at cell centers from the normal components of the velocity
  !  at cell faces (or cell edges in 2D), or for computing the full velocity at an
  !  immersed boundary image point based on the normal velocity at several faces and
  !  the full velocity at boundary points (e.g., a no-slip boundary condition).
  !
  ! The last two routines compute coefficients that can be used to interpolate
  !  a vector function to the destination point given "function dot unitVector" values
  !  at non-tangent source point and "dFunction/dn dot unitVector" values at
  !  tangent source point.  These routines are useful, for example, for computing the 
  !  full velocity at an immersed boundary image point based on the normal velocity at
  !  several faces, the normal velocity at boundary points (e.g., u dot n = 0 for a
  !  no-penetration boundary condition on a fixed boundary), and the normal derivative
  !  of the tangential components of velocity at the boundary points (e.g., a free-slip
  !  boundary condition).
  ! Pseudocode for function reconstruction at the destinationPoint is as follows
  !  dirichlet:  functionAtDestination_i = sum_j(functionDotUnitVectorAtSources_j*coefficients_j,i)
  !    for i = x,y,z
  !  tangentNeumann:    functionAtDestination_i &
  !    = sum_(j where .not isTangent) (functionDotUnitVectorAtSources_j*coefficients_j,i) &
  !    + sum_(j where isTangent) ((dFunction/dn dot UnitVector)_j*coefficients_j,i)
  !    for i = x,y,z
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  public :: rbfInterp_func_3D_vec_const_dir_compCoeffs, &
    rbfInterp_func_3DPlane_vec_const_dir_compCoeffs!, &
    !rbfInterp_func_3D_vec_const_tanNeu_compCoeffs, &
    !rbfInterp_func_3DPlane_vec_const_tanNeu_compCoeffs

  contains

  subroutine rbfInterp_initialize(grid)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: compute geometric fields that will be potentially useful for calling
  !          the interpolation routines
  !
  ! Input: the grid
  !
  ! Output: 
  !  edgeNormalVectors - the unit vector at the center of each edge tangent to the sphere
  !  cellTangentPlane - 2 orthogonal unit vectors in the tangent plane of each cell
  !                     The first unit vector is chosen to point toward the center of the first
  !                     edge on the cell.
  !  localVerticalUnitVectors - the unit normal vector of the tangent plane at the center 
  !                             of each cell
  !       
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    implicit none

    type (grid_meta), intent(inout) :: grid 

    integer :: nCells, nEdges
    integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell
    integer :: iEdge, iCell, cell1, cell2
    real(kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xEdge, yEdge, zEdge
    real(kind=RKIND), dimension(:,:), pointer :: localVerticalUnitVectors, edgeNormalVectors
    real(kind=RKIND), dimension(:,:,:), pointer :: cellTangentPlane
    real(kind=RKIND), dimension(3) :: xHatPlane, yHatPlane, rHat
    real(kind=RKIND) :: normalDotRHat
    logical :: on_a_sphere

    xCell       => grid % xCell % array
    yCell       => grid % yCell % array
    zCell       => grid % zCell % array
    xEdge       => grid % xEdge % array
    yEdge       => grid % yEdge % array
    zEdge       => grid % zEdge % array
    cellsOnEdge => grid % cellsOnEdge % array
    edgesOnCell => grid % edgesOnCell % array
    nCells      = grid % nCells
    nEdges      = grid % nEdges
    on_a_sphere = grid % on_a_sphere

    localVerticalUnitVectors => grid % localVerticalUnitVectors % array
    edgeNormalVectors => grid % edgeNormalVectors % array
    cellTangentPlane => grid % cellTangentPlane % array

    ! init arrays
    edgeNormalVectors = 0
    localVerticalUnitVectors = 0

    ! loop over all cells to be solved on this block
    do iCell=1,nCells
      if(on_a_sphere) then
        localVerticalUnitVectors(1,iCell) = xCell(iCell)
        localVerticalUnitVectors(2,iCell) = yCell(iCell)
        localVerticalUnitVectors(3,iCell) = zCell(iCell)
        call unit_vec_in_R3(localVerticalUnitVectors(:,iCell))
      else ! on a plane
        localVerticalUnitVectors(:,iCell) = (/ 0., 0., 1. /)
      end if
    end do

    do iEdge = 1,nEdges
      iCell = cellsOnEdge(1,iEdge) ! the normal vector points from the first cell toward the edge
      if(iCell == nCells+1) then ! this is a boundary edge
        ! the first cell bordering this edge is not real, use the second cell
        !  The normal should always point outward at boundaries, away from the valid cell center
        iCell = cellsOnEdge(2,iEdge)
      end if
      ! the normal points from the cell location to the edge location
      edgeNormalVectors(1,iEdge) = xEdge(iEdge) - xCell(iCell)
      edgeNormalVectors(2,iEdge) = yEdge(iEdge) - yCell(iCell)
      edgeNormalVectors(3,iEdge) = zEdge(iEdge) - zCell(iCell)
      call unit_vec_in_R3(edgeNormalVectors(:,iEdge))
    end do

    do iCell=1,nCells
      iEdge = edgesOnCell(1,iCell)
      ! xHat and yHat are a local basis in the plane of the horizontal cell
      ! we arbitrarily choose xHat to point toward the first edge
      rHat = localVerticalUnitVectors(:,iCell)
      normalDotRHat = sum(edgeNormalVectors(:,iEdge)*rHat)
      xHatPlane = edgeNormalVectors(:,iEdge) - normalDotRHat*rHat
      call unit_vec_in_R3(xHatPlane)

      call cross_product_in_R3(rHat, xHatPlane, yHatPlane)
      call unit_vec_in_R3(yHatPlane) ! just to be sure...
      cellTangentPlane(:,1,iCell) = xHatPlane
      cellTangentPlane(:,2,iCell) = yHatPlane
    end do

  end subroutine rbfInterp_initialize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: Compute interpolation coefficients in 2D that can be used to
  !  reconstruct a given scalar function at varying locations. This is useful
  !  for finding the location on the the RBF reconstruction of a function
  !  (e.g., a height field) that minimizes the distance to a point in 3D space.
  !  The reconstruction is performed with basis functions that are RBFs and constant 
  ! Input:
  !  pointCount - the number of "source" points and functionValues supplied
  !  coeffCount - the size of coefficients, must be at least pointCount + 1
  !  points - the location of the "source" points in the 2D space where the values of
  !    the function are known
  !  fieldValues - the values of the function of interest at the points
  !  alpha - a constant that give the characteristic length scale of the RBFs,
  !    should be on the order of the distance between points
  ! Output:
  !  coefficients - the coefficients needed to perform interpolation of the funciton
  !    at destination points yet to be specified
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine rbfInterp_loc_2D_sca_const_compCoeffs(pointCount, coeffCount, &
    points, fieldValues, alpha, coefficients)
 
    integer, intent(in) :: pointCount, coeffCount
    real(kind=RKIND), dimension(pointCount,2), intent(in) :: points
    real(kind=RKIND), dimension(pointCount), intent(in) :: fieldValues
    real(kind=RKIND), intent(in) :: alpha
    real(kind=RKIND), dimension(coeffCount), intent(out) :: coefficients

    integer :: i, j, matrixSize
    real(kind=RKIND), dimension(pointCount+1,pointCount+1) :: matrix
    real(kind=RKIND), dimension(pointCount+1) :: rhs
    integer, dimension(pointCount+1) :: pivotIndices
    real(kind=RKIND) :: rSquared

    matrixSize = pointCount+1
    coefficients = 0.0
    matrix = 0.0
    rhs = 0.0

    rhs(1:pointCount) = fieldValues

    do j=1,pointCount
      do i=j,pointCount
        rSquared = sum((points(j,:) - points(i,:))**2)/alpha**2
        matrix(i,j) = evaluateRBF(rSquared)
        matrix(j,i) = matrix(i,j)
      end do
    end do
    do j=1,pointCount
      matrix(pointCount+1,j) = 1.0
      matrix(j,pointCount+1) = 1.0
    end do

    call LEGS(matrix(1:matrixSize,1:matrixSize), matrixSize, rhs(1:matrixSize), &
      coefficients(1:matrixSize), pivotIndices(1:matrixSize))

  end subroutine rbfInterp_loc_2D_sca_const_compCoeffs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: Compute interpolation coefficients in 2D that can be used to
  !  reconstruct a given scalar function at varying locations. This is useful
  !  for finding the location on the the RBF reconstruction of a function
  !  (e.g., a height field) that minimizes the distance to a point in 3D space.
  !  The reconstruction is performed with basis functions that are RBFs plus constant
  !  and linear 
  ! Input:
  !  pointCount - the number of "source" points and functionValues supplied
  !  coeffCount - the size of coefficients, must be at least pointCount + 3
  !  points - the location of the "source" points in the 2D space where the values of
  !    the function are known
  !  fieldValues - the values of the function of interest at the points
  !  alpha - a constant that give the characteristic length scale of the RBFs,
  !    should be on the order of the distance between points
  ! Output:
  !  coefficients - the coefficients needed to perform interpolation of the funciton
  !    at destination points yet to be specified
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine rbfInterp_loc_2D_sca_lin_compCoeffs(pointCount, coeffCount, &
    points, fieldValues, alpha, coefficients)
 
    integer, intent(in) :: pointCount, coeffCount
    real(kind=RKIND), dimension(pointCount,2), intent(in) :: points
    real(kind=RKIND), dimension(pointCount), intent(in) :: fieldValues
    real(kind=RKIND), intent(in) :: alpha
    real(kind=RKIND), dimension(coeffCount), intent(out) :: coefficients

    integer :: i, j, matrixSize
    real(kind=RKIND), dimension(pointCount+3,pointCount+3) :: matrix
    real(kind=RKIND), dimension(pointCount+3) :: rhs
    integer, dimension(pointCount+3) :: pivotIndices
    real(kind=RKIND) :: rSquared

    coefficients = 0.0
    matrix = 0.0
    rhs = 0.0

    rhs(1:pointCount) = fieldValues

    do j=1,pointCount
      do i=j,pointCount
        rSquared = sum((points(j,:) - points(i,:))**2)/alpha**2
        matrix(i,j) = evaluateRBF(rSquared)
        matrix(j,i) = matrix(i,j)
      end do
    end do
    matrixSize = pointCount+3
    do j=1,pointCount
      matrix(pointCount+1,j) = 1.0
      matrix(pointCount+2,j) = points(j,1)
      matrix(pointCount+3,j) = points(j,2)
      matrix(j,pointCount+1:pointCount+3) = matrix(pointCount+1:pointCount+3, j)
    end do
    call LEGS(matrix(1:matrixSize,1:matrixSize), matrixSize, rhs(1:matrixSize), &
      coefficients(1:matrixSize), pivotIndices(1:matrixSize))

  end subroutine rbfInterp_loc_2D_sca_lin_compCoeffs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: Evalute a scalar function in 2D using coefficients computed in
  !  rbfInterp_loc_2D_sca_const_compCoeffs.  This 
  !  function can be called repeatedly with different destination points
  !  to quickly evaluate the interpolating function using the same
  !  coefficients.  This is useful for finding the location on the the 
  !  RBF reconstruction of a function (e.g., a height field) that minimizes
  !  the distance to a point in 3D space. The reconstruction is performed
  !  with basis functions that are RBFs and constant 
  ! Input:
  !  fieldCount - the number fields to be evaluated.  This is useful for reconstructing,
  !    for example, the x-, y- and z-components of a vector field at the same
  !    point in 2D
  !  coeffCount - the size of coefficients, must be at least pointCount + 1
  !  pointCount - the number of "source" points and functionValues supplied
  !  coefficients - the coefficients needed to perform interpolation of the funciton
  !    at the evaluationPoint
  !  evaluationPoint - the point in 2D where the function is to be reconstructed
  !  points - the location of the "source" points in the 2D space where the values of
  !    the function are known
  !  alpha - a constant that give the characteristic length scale of the RBFs,
  !    should be on the order of the distance between points
  ! Output:
  !  derivs - the value of the function, the 2 components of its Jacobian and
  !    the 3 unique components of its Hessian at the evaluationPoint
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine rbfInterp_loc_2D_sca_const_evalWithDerivs(fieldCount, coeffCount, &
    pointCount, coefficients, evaluationPoint, points, alpha, derivs)
    integer, intent(in) :: fieldCount, coeffCount, pointCount
    real(kind=RKIND), dimension(coeffCount, fieldCount), intent(in) :: coefficients
    real(kind=RKIND), dimension(2), intent(in) :: evaluationPoint
    real(kind=RKIND), dimension(pointCount,2), intent(in) :: points
    real(kind=RKIND), intent(in) :: alpha

    real(kind=RKIND), dimension(6,fieldCount), intent(out) :: derivs

    integer :: pointIndex
    real(kind=RKIND) :: x, y, rSquared, rbfValue, rbfDerivOverR, rbfSecondDeriv

    derivs = 0.0
    do pointIndex = 1, pointCount
      x = (evaluationPoint(1) - points(pointIndex,1))
      y = (evaluationPoint(2) - points(pointIndex,2))
      rSquared = x**2 + y**2
      call evaluateRBFAndDerivs(rSquared/alpha**2,  rbfValue, rbfDerivOverR, rbfSecondDeriv)
      rbfDerivOverR = rbfDerivOverR/alpha**2
      rbfSecondDeriv = rbfSecondDeriv/alpha**2
      if(rSquared/alpha**2 < 1e-7) then
        ! terms 2,3 and 5 are zero by radial symmetry of the RBF
        derivs(1,:) = derivs(1,:) + coefficients(pointIndex,:)*rbfValue
        derivs(4,:) = derivs(4,:) + coefficients(pointIndex,:)*rbfSecondDeriv
        derivs(6,:) = derivs(6,:) + coefficients(pointIndex,:)*rbfSecondDeriv
      else
        call evaluateRBFAndDerivs(rSquared,  rbfValue, rbfDerivOverR, rbfSecondDeriv)
        derivs(1,:) = derivs(1,:) + coefficients(pointIndex,:)*rbfValue
        derivs(2,:) = derivs(2,:) + coefficients(pointIndex,:)*rbfDerivOverR*x
        derivs(3,:) = derivs(3,:) + coefficients(pointIndex,:)*rbfDerivOverR*y
        derivs(4,:) = derivs(4,:) + coefficients(pointIndex,:) &
          * (rbfSecondDeriv*x**2 + rbfDerivOverR*y**2)/rSquared
        derivs(5,:) = derivs(5,:) + coefficients(pointIndex,:) &
          * (rbfSecondDeriv - rbfDerivOverR)*x*y/rSquared
        derivs(6,:) = derivs(6,:) + coefficients(pointIndex,:) &
          * (rbfSecondDeriv*y**2 + rbfDerivOverR*x**2)/rSquared
      end if
    end do
    derivs(1,:) = derivs(1,:) + coefficients(pointCount+1,:)
  end subroutine rbfInterp_loc_2D_sca_const_evalWithDerivs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: Evalute a scalar function in 2D using coefficients computed in
  !  rbfInterp_loc_2D_sca_const_compCoeffs.  This 
  !  function can be called repeatedly with different destination points
  !  to quickly evaluate the interpolating function using the same
  !  coefficients.  This is useful for finding the location on the the 
  !  RBF reconstruction of a function (e.g., a height field) that minimizes
  !  the distance to a point in 3D space. The reconstruction is performed
  !  with basis functions that are RBFs, constant and linear
  ! Input:
  !  fieldCount - the number fields to be evaluated.  This is useful for reconstructing,
  !    for example, the x-, y- and z-components of a vector field at the same
  !    point in 2D
  !  coeffCount - the size of coefficients, must be at least pointCount + 1
  !  pointCount - the number of "source" points and functionValues supplied
  !  coefficients - the coefficients needed to perform interpolation of the funciton
  !    at the evaluationPoint
  !  evaluationPoint - the point in 2D where the function is to be reconstructed
  !  points - the location of the "source" points in the 2D space where the values of
  !    the function are known
  !  alpha - a constant that give the characteristic length scale of the RBFs,
  !    should be on the order of the distance between points
  ! Output:
  !  derivs - the value of the function, the 2 components of its Jacobian and
  !    the 3 unique components of its Hessian at the evaluationPoint
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine rbfInterp_loc_2D_sca_lin_evalWithDerivs(fieldCount, coeffCount, &
    pointCount, coefficients, evaluationPoint, points, alpha, derivs)
    integer, intent(in) :: fieldCount, coeffCount, pointCount
    real(kind=RKIND), dimension(coeffCount, fieldCount), intent(in) :: coefficients
    real(kind=RKIND), dimension(2), intent(in) :: evaluationPoint
    real(kind=RKIND), dimension(pointCount,2), intent(in) :: points
    real(kind=RKIND), intent(in) :: alpha

    real(kind=RKIND), dimension(6,fieldCount), intent(out) :: derivs

    integer :: pointIndex
    real(kind=RKIND) :: x, y, rSquared, rbfValue, rbfDerivOverR, rbfSecondDeriv

    derivs = 0.0
    do pointIndex = 1, pointCount
      x = (evaluationPoint(1) - points(pointIndex,1))
      y = (evaluationPoint(2) - points(pointIndex,2))
      rSquared = x**2 + y**2
      call evaluateRBFAndDerivs(rSquared/alpha**2,  rbfValue, rbfDerivOverR, rbfSecondDeriv)
      rbfDerivOverR = rbfDerivOverR/alpha**2
      rbfSecondDeriv = rbfSecondDeriv/alpha**2
      if(rSquared/alpha**2 < 1e-7) then
        ! terms 2,3 and 5 are zero by radial symmetry of the RBF
        derivs(1,:) = derivs(1,:) + coefficients(pointIndex,:)*rbfValue
        derivs(4,:) = derivs(4,:) + coefficients(pointIndex,:)*rbfSecondDeriv
        derivs(6,:) = derivs(6,:) + coefficients(pointIndex,:)*rbfSecondDeriv
      else
        call evaluateRBFAndDerivs(rSquared,  rbfValue, rbfDerivOverR, rbfSecondDeriv)
        derivs(1,:) = derivs(1,:) + coefficients(pointIndex,:)*rbfValue
        derivs(2,:) = derivs(2,:) + coefficients(pointIndex,:)*rbfDerivOverR*x
        derivs(3,:) = derivs(3,:) + coefficients(pointIndex,:)*rbfDerivOverR*y
        derivs(4,:) = derivs(4,:) + coefficients(pointIndex,:) &
          * (rbfSecondDeriv*x**2 + rbfDerivOverR*y**2)/rSquared
        derivs(5,:) = derivs(5,:) + coefficients(pointIndex,:) &
          * (rbfSecondDeriv - rbfDerivOverR)*x*y/rSquared
        derivs(6,:) = derivs(6,:) + coefficients(pointIndex,:) &
          * (rbfSecondDeriv*y**2 + rbfDerivOverR*x**2)/rSquared
      end if
    end do
    derivs(1,:) = derivs(1,:) + coefficients(pointCount+1,:) &
      + coefficients(pointCount+2,:)*evaluationPoint(1) &
      + coefficients(pointCount+3,:)*evaluationPoint(2)
    derivs(2,:) = derivs(2,:) + coefficients(pointCount+2,:)
    derivs(3,:) = derivs(3,:) + coefficients(pointCount+3,:)

  end subroutine rbfInterp_loc_2D_sca_lin_evalWithDerivs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: Compute interpolation coefficients in 3D that can be used to
  !  interpolate a number of scalar functions at a given locations. This is useful
  !  if the interpolation location does not change with time, or if several
  !  fields are to be interpolated at a given time step.  (If both the fields
  !  and the interpolation locations vary with time, there is no clear advantage in
  !  using either this method or the method for 2D interpoaltion above; for simplicity
  !  and because we foresee more uses for the method of this subroutine, we have not
  !  implemented a 3D version of the fixed field, variable interpolation location method
  !  as we have in 2D.) Coefficients are produced for handling Dirichlet boundary
  !  conditions (or no boundaries).  The interpolation is performed with basis functions
  !  that are RBFs plus a constant.
  ! Input:
  !  pointCount - the number of "source" points and functionValues supplied
  !  sourcePoints - the location of the "source" points in the 3D space where the values of
  !    the function are known
  !  destinationPoint - the point where the interpolation will be performed
  !  alpha - a constant that give the characteristic length scale of the RBFs,
  !    should be on the order of the distance between points
  ! Output:
  !  coefficients - the coefficients used to interpolate a function with Dirichlet
  !    boundary conditions to the specified destinationPoint
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine rbfInterp_func_3D_sca_const_dir_compCoeffs( &
    pointCount, sourcePoints, destinationPoint, alpha, coefficients)

    integer, intent(in) :: pointCount
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint
    real(kind=RKIND), intent(in) :: alpha
    real(kind=RKIND), dimension(pointCount), intent(out) :: coefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(:,:), pointer :: dirichletMatrix
    real(kind=RKIND), dimension(:), pointer :: rhs, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+1 !! 1 extra space for constant 

    allocate(dirichletMatrix(matrixSize,matrixSize))  
    allocate(rhs(matrixSize))  
    allocate(coeffs(matrixSize))  
    allocate(pivotIndices(matrixSize))  

    dirichletMatrix = 0.0
    rhs = 0.0
    coeffs = 0.0

    call setUpScalarRBFDirichletMatrixAndRHS(pointCount, sourcePoints, destinationPoint, &
      alpha, dirichletMatrix(1:pointCount,1:pointCount), rhs(1:pointCount))

    do i = 1, pointCount
      dirichletMatrix(i,pointCount+1) = 1.0
      dirichletMatrix(pointCount+1,i) = dirichletMatrix(i,pointCount+1)
    end do

    rhs(pointCount+1) = 1.0

    ! solve each linear system
    call LEGS(dirichletMatrix, matrixSize, rhs, coeffs, pivotIndices)
    coefficients = coeffs(1:pointCount)

    deallocate(dirichletMatrix)  
    deallocate(rhs)  
    deallocate(coeffs)  
    deallocate(pivotIndices) 

  end subroutine rbfInterp_func_3D_sca_const_dir_compCoeffs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: Compute interpolation coefficients in a plane in 3D that can be used to
  !  interpolate a number of scalar functions at a given locations. This is useful
  !  if the interpolation location does not change with time, or if several
  !  fields are to be interpolated at a given time step.  (If both the fields
  !  and the interpolation locations vary with time, there is no clear advantage in
  !  using either this method or the method for 2D interpoaltion above; for simplicity
  !  and because we foresee more uses for the method of this subroutine, we have not
  !  implemented a 3D version of the fixed field, variable interpolation location method
  !  as we have in 2D.) Coefficients are produced for handling both Dirichlet (or no)
  !  boundary conditions.  The interpolation is performed with basis functions that are
  !  RBFs plus constant and linear.  All points are projected into the plane given by the
  !  planeBasisVectors.
  ! Input:
  !  pointCount - the number of "source" points and functionValues supplied
  !  sourcePoints - the location of the "source" points in the 3D space where the values of
  !    the function are known.  The points will be projected into the plane given by 
  !    planeBasisVectors
  !  destinationPoint - the point in 3D where the interpolation will be performed.  The
  !    destinationPoint will be projected into the plane given by planeBasisVectors.
  !  alpha - a constant that give the characteristic length scale of the RBFs,
  !    should be on the order of the distance between points
  !  planeBasisVectors - the basis fectors for the plane where interpolation is performed.
  !    All points are projected into this plane. 
  ! Output:
  !  coefficients - the coefficients used to interpolate a function with Dirichlet
  !    boundary conditions to the specified destinationPoint
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine rbfInterp_func_3DPlane_sca_lin_dir_compCoeffs( &
    pointCount, sourcePoints, destinationPoint, &
    alpha, planeBasisVectors, coefficients)

    integer, intent(in) :: pointCount
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint
    real(kind=RKIND), intent(in) :: alpha
    real(kind=RKIND), dimension(2,3) :: planeBasisVectors
    real(kind=RKIND), dimension(pointCount), intent(out) :: coefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(:,:), pointer :: dirichletMatrix
    real(kind=RKIND), dimension(:), pointer :: rhs, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+3 !! 3 extra space for constant and 2 planar dimensions

    allocate(dirichletMatrix(matrixSize,matrixSize))  
    allocate(rhs(matrixSize))  
    allocate(coeffs(matrixSize))  
    allocate(pivotIndices(matrixSize))  

    dirichletMatrix = 0.0
    rhs = 0.0
    coeffs = 0.0

    call setUpScalarRBFDirichletMatrixAndRHS(pointCount, sourcePoints, destinationPoint, &
      alpha, dirichletMatrix(1:pointCount,1:pointCount), rhs(1:pointCount))

    do i = 1, pointCount
      dirichletMatrix(i,pointCount+1) = 1.0
      dirichletMatrix(i,pointCount+2) = sum(sourcePoints(i,1:3)*planeBasisVectors(1,:))
      dirichletMatrix(i,pointCount+3) = sum(sourcePoints(i,1:3)*planeBasisVectors(2,:))
      dirichletMatrix(pointCount+1:pointCount+3,i) &
        = dirichletMatrix(i,pointCount+1:pointCount+3)
    end do

    rhs(pointCount+1) = 1.0
    rhs(pointCount+2) = sum(destinationPoint(1:3)*planeBasisVectors(1,:))
    rhs(pointCount+3) = sum(destinationPoint(1:3)*planeBasisVectors(2,:))

    ! solve each linear system
    call LEGS(dirichletMatrix, matrixSize, rhs, coeffs, pivotIndices)
    coefficients = coeffs(1:pointCount)

    deallocate(dirichletMatrix)  
    deallocate(rhs)  
    deallocate(coeffs)  
    deallocate(pivotIndices) 

  end subroutine rbfInterp_func_3DPlane_sca_lin_dir_compCoeffs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: Compute interpolation coefficients in 3D that can be used to
  !  interpolate a number of scalar functions at a given locations. This is useful
  !  if the interpolation location does not change with time, or if several
  !  fields are to be interpolated at a given time step.  (If both the fields
  !  and the interpolation locations vary with time, there is no clear advantage in
  !  using either this method or the method for 2D interpoaltion above; for simplicity
  !  and because we foresee more uses for the method of this subroutine, we have not
  !  implemented a 3D version of the fixed field, variable interpolation location method
  !  as we have in 2D.) Coefficients are produced for handling both Dirichlet (or no)
  !  boundary conditions.  The interpolation is performed with basis functions that are
  !  RBFs plus constant and linear.
  ! Input:
  !  pointCount - the number of "source" points and functionValues supplied
  !  sourcePoints - the location of the "source" points in the 3D space where the values of
  !    the function are known
  !  destinationPoint - the point where the interpolation will be performed
  !  alpha - a constant that give the characteristic length scale of the RBFs,
  !    should be on the order of the distance between points
  ! Output:
  !  coefficients - the coefficients used to interpolate a function with Dirichlet
  !    boundary conditions to the specified destinationPoint
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine rbfInterp_func_3D_sca_lin_dir_compCoeffs(pointCount, &
    sourcePoints, destinationPoint, alpha, coefficients)

    integer, intent(in) :: pointCount
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint
    real(kind=RKIND), intent(in) :: alpha
    real(kind=RKIND), dimension(pointCount), intent(out) :: coefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(:,:), pointer :: dirichletMatrix
    real(kind=RKIND), dimension(:), pointer :: rhs, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+4 !! 4 extra space for constant and linear in 3D

    allocate(dirichletMatrix(matrixSize,matrixSize))  
    allocate(rhs(matrixSize))  
    allocate(coeffs(matrixSize))  
    allocate(pivotIndices(matrixSize))  

    dirichletMatrix = 0.0
    rhs = 0.0
    coeffs = 0.0

    call setUpScalarRBFDirichletMatrixAndRHS(pointCount, sourcePoints, destinationPoint, &
      alpha, dirichletMatrix(1:pointCount,1:pointCount), rhs(1:pointCount))

    do i = 1, pointCount
      dirichletMatrix(i,pointCount+1) = 1.0
      dirichletMatrix(i,pointCount+2:pointCount+4) = sourcePoints(i,1:3)
      dirichletMatrix(pointCount+1:pointCount+4,i) &
        = dirichletMatrix(i,pointCount+1:pointCount+4)
    end do

    rhs(pointCount+1) = 1.0
    rhs(pointCount+2:pointCount+4) = destinationPoint(1:3)

    ! solve each linear system
    call LEGS(dirichletMatrix, matrixSize, rhs, coeffs, pivotIndices)
    coefficients = coeffs(1:pointCount)

    deallocate(dirichletMatrix)  
    deallocate(rhs)  
    deallocate(coeffs)  
    deallocate(pivotIndices) 

  end subroutine rbfInterp_func_3D_sca_lin_dir_compCoeffs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: Compute interpolation coefficients in 3D that can be used to
  !  interpolate a number of scalar functions at a given locations. This is useful
  !  if the interpolation location does not change with time, or if several
  !  fields are to be interpolated at a given time step.  (If both the fields
  !  and the interpolation locations vary with time, there is no clear advantage in
  !  using either this method or the method for 2D interpoaltion above; for simplicity
  !  and because we foresee more uses for the method of this subroutine, we have not
  !  implemented a 3D version of the fixed field, variable interpolation location method
  !  as we have in 2D.) Coefficients are produced for handling both Dirichlet and Neumann
  !  boundary conditions.  The interpolation is performed with basis functions that are
  !  RBFs plus a constant.
  ! Input:
  !  pointCount - the number of "source" points and functionValues supplied
  !  sourcePoints - the location of the "source" points in the 3D space where the values of
  !    the function are known
  !  isInterface - a logical array indicating which of the source points (if any) are at
  !    at the domain interface.  These points and their normals will be used to compute the
  !    neumannCoefficients below
  !  interfaceNormals - a 3D normal vector for each sourcePoint.  These vectors are only used
  !    at points where isInterface == .true., and can take arbitrary values elsewehere.  The
  !    normal vector is used to compute coefficients for the normal derivative of the
  !    interpolating function in order to impose the Neumann Boundary condition
  !  destinationPoint - the point where the interpolation will be performed
  !  alpha - a constant that give the characteristic length scale of the RBFs,
  !    should be on the order of the distance between points
  ! Output:
  !  dirichletCoefficients - the coefficients used to interpolate a function with Dirichlet
  !    boundary conditions to the specified destinationPoint
  !  neumannCoefficients - the coefficients used to interpolate a function with Neumann
  !    boundary conditions to the specified destinationPoint
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine rbfInterp_func_3D_sca_const_dirNeu_compCoeffs( &
    pointCount, sourcePoints, isInterface, interfaceNormals, destinationPoint, &
    alpha, dirichletCoefficients, neumannCoefficients)

    integer, intent(in) :: pointCount
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints
    logical, dimension(pointCount), intent(in) :: isInterface
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: interfaceNormals
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint
    real(kind=RKIND), intent(in) :: alpha
    real(kind=RKIND), dimension(pointCount), intent(out) :: &
      dirichletCoefficients, neumannCoefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(:,:), pointer :: dirichletMatrix, neumannMatrix
    real(kind=RKIND), dimension(:), pointer :: rhs, rhsCopy, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+1 !! 1 extra space for constant 

    allocate(dirichletMatrix(matrixSize,matrixSize))  
    allocate(neumannMatrix(matrixSize,matrixSize))  
    allocate(rhs(matrixSize))  
    allocate(rhsCopy(matrixSize))  
    allocate(coeffs(matrixSize))  
    allocate(pivotIndices(matrixSize))  

    dirichletMatrix = 0.0
    neumannMatrix = 0.0
    rhs = 0.0
    rhsCopy = 0.0
    coeffs = 0.0

    call setUpScalarRBFMatrixAndRHS(pointCount, &
      sourcePoints, isInterface, interfaceNormals, destinationPoint, &
      alpha, dirichletMatrix(1:pointCount,1:pointCount), &
      neumannMatrix(1:pointCount,1:pointCount), rhs(1:pointCount))

    do i = 1, pointCount
      dirichletMatrix(i,pointCount+1) = 1.0
      if(isInterface(i)) then
        neumannMatrix(i,pointCount+1) = 0.0
      else
        neumannMatrix(i,pointCount+1) = dirichletMatrix(i,pointCount+1)
      end if
      dirichletMatrix(pointCount+1,i) = dirichletMatrix(i,pointCount+1)
      neumannMatrix(pointCount+1,i) = neumannMatrix(i,pointCount+1)
    end do

    rhs(pointCount+1) = 1.0

    ! solve each linear system
    rhsCopy = rhs
    call LEGS(dirichletMatrix, matrixSize, rhs, coeffs, pivotIndices)
    dirichletCoefficients = coeffs(1:pointCount)

    call LEGS(neumannMatrix, matrixSize, rhsCopy, coeffs, pivotIndices)
    neumannCoefficients = coeffs(1:pointCount)

    deallocate(dirichletMatrix)  
    deallocate(neumannMatrix)  
    deallocate(rhs)  
    deallocate(rhsCopy)  
    deallocate(coeffs)  
    deallocate(pivotIndices) 

  end subroutine rbfInterp_func_3D_sca_const_dirNeu_compCoeffs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: Compute interpolation coefficients in a plane in 3D that can be used to
  !  interpolate a number of scalar functions at a given locations. This is useful
  !  if the interpolation location does not change with time, or if several
  !  fields are to be interpolated at a given time step.  (If both the fields
  !  and the interpolation locations vary with time, there is no clear advantage in
  !  using either this method or the method for 2D interpoaltion above; for simplicity
  !  and because we foresee more uses for the method of this subroutine, we have not
  !  implemented a 3D version of the fixed field, variable interpolation location method
  !  as we have in 2D.) Coefficients are produced for handling both Dirichlet and Neumann
  !  boundary conditions.  The interpolation is performed with basis functions that are
  !  RBFs plus constant and linear.  All points are projected into the plane given by the
  !  planeBasisVectors.
  ! Input:
  !  pointCount - the number of "source" points and functionValues supplied
  !  sourcePoints - the location of the "source" points in the 3D space where the values of
  !    the function are known.  The sourcePoints will be projected into the plane given by
  !    planeBasisVectors
  !  isInterface - a logical array indicating which of the source points (if any) are at
  !    at the domain interface.  These points and their normals will be used to compute the
  !    neumannCoefficients below
  !  interfaceNormals - a 3D normal vector for each sourcePoint.  These vectors are only used
  !    at points where isInterface == .true., and can take arbitrary values elsewehere.  The
  !    normal vector is used to compute coefficients for the normal derivative of the
  !    interpolating function in order to impose the Neumann Boundary condition
  !  destinationPoint - the point in 3D where the interpolation will be performed.  The
  !    destinationPoint will be projected into the plane given by planeBasisVectors.
  !  alpha - a constant that give the characteristic length scale of the RBFs,
  !    should be on the order of the distance between points
  !  planeBasisVectors - the basis fectors for the plane where interpolation is performed.
  !    All points are projected into this plane. 
  ! Output:
  !  dirichletCoefficients - the coefficients used to interpolate a function with Dirichlet
  !    boundary conditions to the specified destinationPoint
  !  neumannCoefficients - the coefficients used to interpolate a function with Neumann
  !    boundary conditions to the specified destinationPoint
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine rbfInterp_func_3DPlane_sca_lin_dirNeu_compCoeffs( &
    pointCount, sourcePoints, isInterface, interfaceNormals, destinationPoint, &
    alpha, planeBasisVectors, dirichletCoefficients, neumannCoefficients)

    integer, intent(in) :: pointCount
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints
    logical, dimension(pointCount), intent(in) :: isInterface
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: interfaceNormals
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint
    real(kind=RKIND), intent(in) :: alpha
    real(kind=RKIND), dimension(2,3) :: planeBasisVectors
    real(kind=RKIND), dimension(pointCount), intent(out) :: &
      dirichletCoefficients, neumannCoefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(:,:), pointer :: dirichletMatrix, neumannMatrix
    real(kind=RKIND), dimension(:), pointer :: rhs, rhsCopy, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+3 !! 3 extra space for constant and 2 planar dimensions

    allocate(dirichletMatrix(matrixSize,matrixSize))  
    allocate(neumannMatrix(matrixSize,matrixSize))  
    allocate(rhs(matrixSize))  
    allocate(rhsCopy(matrixSize))  
    allocate(coeffs(matrixSize))  
    allocate(pivotIndices(matrixSize))  

    dirichletMatrix = 0.0
    neumannMatrix = 0.0
    rhs = 0.0
    rhsCopy = 0.0
    coeffs = 0.0

    call setUpScalarRBFMatrixAndRHS(pointCount, &
      sourcePoints, isInterface, interfaceNormals, destinationPoint, &
      alpha, dirichletMatrix(1:pointCount,1:pointCount), &
      neumannMatrix(1:pointCount,1:pointCount), rhs(1:pointCount))

    do i = 1, pointCount
      dirichletMatrix(i,pointCount+1) = 1.0
      dirichletMatrix(i,pointCount+2) = sum(sourcePoints(i,1:3)*planeBasisVectors(1,:))
      dirichletMatrix(i,pointCount+3) = sum(sourcePoints(i,1:3)*planeBasisVectors(2,:))
      if(isInterface(i)) then
        neumannMatrix(i,pointCount+1) = 0.0
        neumannMatrix(i,pointCount+2) = sum(interfaceNormals(i,1:3)*planeBasisVectors(1,:))
        neumannMatrix(i,pointCount+3) = sum(interfaceNormals(i,1:3)*planeBasisVectors(2,:))
      else
        neumannMatrix(i,pointCount+1:pointCount+3) &
          = dirichletMatrix(i,pointCount+1:pointCount+3)
      end if
      dirichletMatrix(pointCount+1:pointCount+3,i) &
        = dirichletMatrix(i,pointCount+1:pointCount+3)
      neumannMatrix(pointCount+1:pointCount+3,i) &
        = neumannMatrix(i,pointCount+1:pointCount+3)
    end do

    rhs(pointCount+1) = 1.0
    rhs(pointCount+2) = sum(destinationPoint(1:3)*planeBasisVectors(1,:))
    rhs(pointCount+3) = sum(destinationPoint(1:3)*planeBasisVectors(2,:))

    ! solve each linear system
    rhsCopy = rhs
    call LEGS(dirichletMatrix, matrixSize, rhs, coeffs, pivotIndices)
    dirichletCoefficients = coeffs(1:pointCount)

    call LEGS(neumannMatrix, matrixSize, rhsCopy, coeffs, pivotIndices)
    neumannCoefficients = coeffs(1:pointCount)

    deallocate(dirichletMatrix)  
    deallocate(neumannMatrix)  
    deallocate(rhs)  
    deallocate(rhsCopy)  
    deallocate(coeffs)  
    deallocate(pivotIndices) 

  end subroutine rbfInterp_func_3DPlane_sca_lin_dirNeu_compCoeffs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: Compute interpolation coefficients in 3D that can be used to
  !  interpolate a number of scalar functions at a given locations. This is useful
  !  if the interpolation location does not change with time, or if several
  !  fields are to be interpolated at a given time step.  (If both the fields
  !  and the interpolation locations vary with time, there is no clear advantage in
  !  using either this method or the method for 2D interpoaltion above; for simplicity
  !  and because we foresee more uses for the method of this subroutine, we have not
  !  implemented a 3D version of the fixed field, variable interpolation location method
  !  as we have in 2D.) Coefficients are produced for handling both Dirichlet and Neumann
  !  boundary conditions.  The interpolation is performed with basis functions that are
  !  RBFs plus constant and linear.
  ! Input:
  !  pointCount - the number of "source" points and functionValues supplied
  !  sourcePoints - the location of the "source" points in the 3D space where the values of
  !    the function are known
  !  isInterface - a logical array indicating which of the source points (if any) are at
  !    at the domain interface.  These points and their normals will be used to compute the
  !    neumannCoefficients below
  !  interfaceNormals - a 3D normal vector for each sourcePoint.  These vectors are only used
  !    at points where isInterface == .true., and can take arbitrary values elsewehere.  The
  !    normal vector is used to compute coefficients for the normal derivative of the
  !    interpolating function in order to impose the Neumann Boundary condition
  !  destinationPoint - the point where the interpolation will be performed
  !  alpha - a constant that give the characteristic length scale of the RBFs,
  !    should be on the order of the distance between points
  ! Output:
  !  dirichletCoefficients - the coefficients used to interpolate a function with Dirichlet
  !    boundary conditions to the specified destinationPoint
  !  neumannCoefficients - the coefficients used to interpolate a function with Neumann
  !    boundary conditions to the specified destinationPoint
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine rbfInterp_func_3D_sca_lin_dirNeu_compCoeffs(pointCount, &
    sourcePoints, isInterface, interfaceNormals, destinationPoint, &
    alpha, dirichletCoefficients, neumannCoefficients)

    integer, intent(in) :: pointCount
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints
    logical, dimension(pointCount), intent(in) :: isInterface
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: interfaceNormals
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint
    real(kind=RKIND), intent(in) :: alpha
    real(kind=RKIND), dimension(pointCount), intent(out) :: &
      dirichletCoefficients, neumannCoefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(:,:), pointer :: dirichletMatrix, neumannMatrix
    real(kind=RKIND), dimension(:), pointer :: rhs, rhsCopy, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+4 !! 4 extra space for constant and linear in 3D

    allocate(dirichletMatrix(matrixSize,matrixSize))  
    allocate(neumannMatrix(matrixSize,matrixSize))  
    allocate(rhs(matrixSize))  
    allocate(rhsCopy(matrixSize))  
    allocate(coeffs(matrixSize))  
    allocate(pivotIndices(matrixSize))  

    dirichletMatrix = 0.0
    neumannMatrix = 0.0
    rhs = 0.0
    rhsCopy = 0.0
    coeffs = 0.0

    call setUpScalarRBFMatrixAndRHS(pointCount, &
      sourcePoints, isInterface, interfaceNormals, destinationPoint, &
      alpha, dirichletMatrix(1:pointCount,1:pointCount), &
      neumannMatrix(1:pointCount,1:pointCount), rhs(1:pointCount))

    do i = 1, pointCount
      dirichletMatrix(i,pointCount+1) = 1.0
      dirichletMatrix(i,pointCount+2:pointCount+4) = sourcePoints(i,1:3)
      if(isInterface(i)) then
        neumannMatrix(i,pointCount+1) = 0.0
        neumannMatrix(i,pointCount+2:pointCount+4) = interfaceNormals(i,1:3)
      else
        neumannMatrix(i,pointCount+1:pointCount+4) &
          = dirichletMatrix(i,pointCount+1:pointCount+4)
      end if
      dirichletMatrix(pointCount+1:pointCount+4,i) &
        = dirichletMatrix(i,pointCount+1:pointCount+4)
      neumannMatrix(pointCount+1:pointCount+4,i) &
        = neumannMatrix(i,pointCount+1:pointCount+4)
    end do

    rhs(pointCount+1) = 1.0
    rhs(pointCount+2:pointCount+4) = destinationPoint(1:3)

    ! solve each linear system
    rhsCopy = rhs
    call LEGS(dirichletMatrix, matrixSize, rhs, coeffs, pivotIndices)
    dirichletCoefficients = coeffs(1:pointCount)

    call LEGS(neumannMatrix, matrixSize, rhsCopy, coeffs, pivotIndices)
    neumannCoefficients = coeffs(1:pointCount)

    deallocate(dirichletMatrix)  
    deallocate(neumannMatrix)  
    deallocate(rhs)  
    deallocate(rhsCopy)  
    deallocate(coeffs)  
    deallocate(pivotIndices) 

  end subroutine rbfInterp_func_3D_sca_lin_dirNeu_compCoeffs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: Compute interpolation coefficients in 3D that can be used to
  !  interpolate a number of vector functions at a given locations. This is useful
  !  if the interpolation location does not change with time, or if several
  !  fields are to be interpolated at a given time step.  (If both the vector fields
  !  and the interpolation locations vary with time, there is no clear advantage in
  !  using either this method or the method for 2D interpoaltion above; for simplicity
  !  and because we foresee more uses for the method of this subroutine, we have not
  !  implemented a 3D version of the fixed field, variable interpolation location method
  !  as we have in 2D.) Coefficients are produced for handling Dirichlet boundary
  !  conditions (or no boundaries).  The interpolation is performed with basis functions
  !  that are RBFs plus a constant.
  ! Input:
  !  pointCount - the number of "source" points and functionValues supplied
  !  sourcePoints - the location of the "source" points in the 3D space where the values of
  !    the function are known
  !  unitVectors - the unit vectors associated with each of the sourcePoints.  Interpolation
  !    is performed by supplying the value of the vector function dotted into each of these unit
  !    vectors.  If multiple unit vectors are supplied at the same sourcePoint, they *must* be
  !    orthogonal for the interpolation to succeed.
  !  destinationPoint - the point where the interpolation will be performed
  !  alpha - a constant that give the characteristic length scale of the RBFs,
  !    should be on the order of the distance between points
  ! Output:
  !  coefficients - the coefficients used to interpolate a function with Dirichlet
  !    boundary conditions to the specified destinationPoint
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine rbfInterp_func_3D_vec_const_dir_compCoeffs(pointCount, &
    sourcePoints, unitVectors, destinationPoint, &
    alpha, coefficients)

    integer, intent(in) :: pointCount
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: unitVectors
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint
    real(kind=RKIND), intent(in) :: alpha
    real(kind=RKIND), dimension(pointCount, 3), intent(out) :: coefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(:,:), pointer :: matrix, matrixCopy
    real(kind=RKIND), dimension(:,:), pointer :: rhs, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+3 ! extra space for constant vector 

    allocate(matrix(matrixSize,matrixSize))  
    allocate(matrixCopy(matrixSize,matrixSize))  
    allocate(rhs(matrixSize,3))  
    allocate(coeffs(matrixSize,3))  
    allocate(pivotIndices(matrixSize))  

    matrix = 0.0
    rhs = 0.0
    coeffs = 0.0

    call setUpVectorDirichletRBFMatrixAndRHS(pointCount, 3, &
      sourcePoints, unitVectors, destinationPoint, &
      alpha, matrix(1:pointCount,1:pointCount), rhs(1:pointCount,:))

    do i = 1, pointCount
      matrix(i,pointCount+1:pointCount+3) = unitVectors(i,:)
      matrix(pointCount+1:pointCount+3,i) &
        = matrix(i,pointCount+1:pointCount+3)
    end do
    do i = 1, 3
      rhs(pointCount+i,i) = 1.0 ! the unit vector in the ith direction
    end do

    ! solve each linear system
    do i = 1, 3
      matrixCopy = matrix
      call LEGS(matrixCopy, matrixSize, rhs(:,i), coeffs(:,i), pivotIndices)
    end do
    coefficients = coeffs(1:pointCount,:)

    deallocate(matrix)  
    deallocate(matrixCopy)  
    deallocate(rhs)  
    deallocate(coeffs)  
    deallocate(pivotIndices) 

  end subroutine rbfInterp_func_3D_vec_const_dir_compCoeffs 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: Compute interpolation coefficients in 3D that can be used to
  !  interpolate a number of vector functions at a given locations. This is useful
  !  if the interpolation location does not change with time, or if several
  !  fields are to be interpolated at a given time step.  (If both the vector fields
  !  and the interpolation locations vary with time, there is no clear advantage in
  !  using either this method or the method for 2D interpoaltion above; for simplicity
  !  and because we foresee more uses for the method of this subroutine, we have not
  !  implemented a 3D version of the fixed field, variable interpolation location method
  !  as we have in 2D.) Coefficients are produced for handling Dirichlet boundary
  !  conditions (or no boundaries).  The interpolation is performed with basis functions
  !  that are RBFs plus a constant.
  ! Input:
  !  pointCount - the number of "source" points and functionValues supplied
  !  sourcePoints - the location of the "source" points in the 3D space where the values of
  !    the function are known.  The sourcePoints are projected into the plane given by
  !    planeBasisVectors
  !  unitVectors - the unit vectors associated with each of the sourcePoints.  Interpolation
  !    is performed by supplying the value of the vector function dotted into each of these unit
  !    vectors.  If multiple unit vectors are supplied at the same sourcePoint, they *must* be
  !    orthogonal for the interpolation to succeed.  The unitVectors are projected into the
  !    plane given by planeBasisVectors
  !  destinationPoint - the point where the interpolation will be performed.  The destinationPoint
  !    is projected into the plane given by planeBasisVectors
  !  alpha - a constant that give the characteristic length scale of the RBFs,
  !    should be on the order of the distance between points
  !  planeBasisVectors - the basis fectors for the plane where interpolation is performed.
  !    All points are projected into this plane. 
  ! Output:
  !  coefficients - the coefficients used to interpolate a function with Dirichlet
  !    boundary conditions to the specified destinationPoint
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine rbfInterp_func_3DPlane_vec_const_dir_compCoeffs(pointCount, &
    sourcePoints, unitVectors, destinationPoint, &
    alpha, planeBasisVectors, coefficients)

    integer, intent(in) :: pointCount
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: unitVectors
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint
    real(kind=RKIND), intent(in) :: alpha
    real(kind=RKIND), dimension(2,3) :: planeBasisVectors
    real(kind=RKIND), dimension(pointCount, 3), intent(out) :: coefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND) :: rSquared, rbfValue, unitVectorDotProduct

    real(kind=RKIND), dimension(pointCount,2) :: planarSourcePoints
    real(kind=RKIND), dimension(pointCount,2) :: planarUnitVectors
    real(kind=RKIND), dimension(2) :: planarDestinationPoint

    real(kind=RKIND), dimension(:,:), pointer :: matrix, matrixCopy
    real(kind=RKIND), dimension(:,:), pointer :: rhs, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+2 ! space for constant vector in plane

    allocate(matrix(matrixSize,matrixSize))  
    allocate(matrixCopy(matrixSize,matrixSize))  
    allocate(rhs(matrixSize,2))  
    allocate(coeffs(matrixSize,2))  
    allocate(pivotIndices(matrixSize)) 

    matrix = 0.0
    rhs = 0.0
    coeffs = 0.0

    do i = 1, pointCount
      planarSourcePoints(i,1) = sum(sourcePoints(i,:)*planeBasisVectors(1,:)) 
      planarSourcePoints(i,2) = sum(sourcePoints(i,:)*planeBasisVectors(2,:)) 
      planarUnitVectors(i,1) = sum(unitVectors(i,:)*planeBasisVectors(1,:)) 
      planarUnitVectors(i,2) = sum(unitVectors(i,:)*planeBasisVectors(2,:)) 
    end do
    planarDestinationPoint(1) = sum(destinationPoint*planeBasisVectors(1,:)) 
    planarDestinationPoint(2) = sum(destinationPoint*planeBasisVectors(2,:)) 

    call setUpVectorDirichletRBFMatrixAndRHS(pointCount, 2, &
      planarSourcePoints, planarUnitVectors, planarDestinationPoint, &
      alpha, matrix(1:pointCount,1:pointCount), rhs(1:pointCount,:))

    do i = 1, pointCount
      matrix(i,pointCount+1:pointCount+2) = planarUnitVectors(i,:) 
      matrix(pointCount+1:pointCount+2,i) = matrix(i,pointCount+1:pointCount+2)
    end do
    do i = 1,2 
      rhs(pointCount+i,i) = 1.0 ! the unit vector in the ith direction
    end do

    ! solve each linear system
    matrixCopy = matrix
    call LEGS(matrix, matrixSize, rhs(:,1), coeffs(:,1), pivotIndices)
    call LEGS(matrixCopy, matrixSize, rhs(:,2), coeffs(:,2), pivotIndices)


    do i = 1,3 
      coefficients(:,i) = planeBasisVectors(1,i)*coeffs(1:pointCount,1) &
        + planeBasisVectors(2,i)*coeffs(1:pointCount,2) 
    end do

    deallocate(matrix)  
    deallocate(matrixCopy)  
    deallocate(rhs)  
    deallocate(coeffs)  
    deallocate(pivotIndices)

  end subroutine rbfInterp_func_3DPlane_vec_const_dir_compCoeffs 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: Compute interpolation coefficients in 3D that can be used to
  !  interpolate a number of vector functions at a given locations. This is useful
  !  if the interpolation location does not change with time, or if several
  !  fields are to be interpolated at a given time step.  (If both the vector fields
  !  and the interpolation locations vary with time, there is no clear advantage in
  !  using either this method or the method for 2D interpoaltion above; for simplicity
  !  and because we foresee more uses for the method of this subroutine, we have not
  !  implemented a 3D version of the fixed field, variable interpolation location method
  !  as we have in 2D.) Coefficients are produced for handling Dirichlet normal /
  !  Neumann tangential boundary conditions (such as free slip).  The interpolation is 
  !  performed with basis functions that are RBFs plus a constant.
  ! Input:
  !  pointCount - the number of "source" points and functionValues supplied
  !  sourcePoints - the location of the "source" points in the 3D space where the values of
  !    the function are known
  !  isTangentToInterface - a logical array indicating which sourcePoints/unitVectors are
  !    tangent to the interface where the boundary condition will be applied.  A Neumann
  !    boundary condition will be applied at these points in these directions.
  !  normalVectorIndex - where isTangentToInterface == .true., the index into unitVectors that
  !    gives the normal vector at the same sourcePoint.  This information is needed to compute
  !    the Neumann boundary condition at this point.
  !  unitVectors - the unit vectors associated with each of the sourcePoints.  Interpolation
  !    is performed by supplying the value of the vector function dotted into each of these unit
  !    vectors.  If multiple unit vectors are supplied at the same sourcePoint, they *must* be
  !    orthogonal for the interpolation to succeed.  A normal vector and two tangential vectors
  !    are needed at each interface point in order to satisfy the Dirichlet normal boundary
  !    condition and the Neumann tangential boundary conditions at these points.
  !  destinationPoint - the point where the interpolation will be performed
  !  alpha - a constant that give the characteristic length scale of the RBFs,
  !    should be on the order of the distance between points
  ! Output:
  !  coefficients - the coefficients used to interpolate a function with Dirichlet
  !    boundary conditions to the specified destinationPoint
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine rbfInterp_func_3D_vec_const_tanNeu_compCoeffs(pointCount, &
    sourcePoints, isTangentToInterface, normalVectorIndex, unitVectors, destinationPoint, &
    alpha, coefficients)

    integer, intent(in) :: pointCount
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints
    logical, dimension(pointCount), intent(in) :: isTangentToInterface
    integer, dimension(pointCount), intent(in) :: normalVectorIndex
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: unitVectors
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint
    real(kind=RKIND), intent(in) :: alpha
    real(kind=RKIND), dimension(pointCount, 3), intent(out) :: coefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(:,:), pointer :: matrix, matrixCopy
    real(kind=RKIND), dimension(:,:), pointer :: rhs, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+3 ! extra space for constant vector 

    allocate(matrix(matrixSize,matrixSize))  
    allocate(matrixCopy(matrixSize,matrixSize))  
    allocate(rhs(matrixSize,3))  
    allocate(coeffs(matrixSize,3))  
    allocate(pivotIndices(matrixSize))  

    matrix = 0.0
    rhs = 0.0
    coeffs = 0.0

    call setUpVectorFreeSlipRBFMatrixAndRHS(pointCount, 3, &
      sourcePoints, isTangentToInterface, normalVectorIndex, unitVectors, destinationPoint, &
      alpha, matrix(1:pointCount,1:pointCount), rhs(1:pointCount,:))

    do i = 1, pointCount
      matrix(pointCount+1:pointCount+3,i) = unitVectors(i,:)
      if(.not. isTangentToInterface(i)) then
        matrix(i,pointCount+1:pointCount+3) = matrix(pointCount+1:pointCount+3,i)
      end if
    end do
    do i = 1, 3
      rhs(pointCount+i,i) = 1.0 ! the unit vector in the ith direction
    end do

    ! solve each linear system
    do i = 1, 3
      matrixCopy = matrix
      call LEGS(matrixCopy, matrixSize, rhs(:,i), coeffs(:,i), pivotIndices)
    end do
    coefficients = coeffs(1:pointCount,:)

    deallocate(matrix)  
    deallocate(matrixCopy)  
    deallocate(rhs)  
    deallocate(coeffs)  
    deallocate(pivotIndices)

  end subroutine rbfInterp_func_3D_vec_const_tanNeu_compCoeffs 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Purpose: Compute interpolation coefficients in 3D that can be used to
  !  interpolate a number of vector functions at a given locations. This is useful
  !  if the interpolation location does not change with time, or if several
  !  fields are to be interpolated at a given time step.  (If both the vector fields
  !  and the interpolation locations vary with time, there is no clear advantage in
  !  using either this method or the method for 2D interpoaltion above; for simplicity
  !  and because we foresee more uses for the method of this subroutine, we have not
  !  implemented a 3D version of the fixed field, variable interpolation location method
  !  as we have in 2D.) Coefficients are produced for handling Dirichlet normal /
  !  Neumann tangential boundary conditions (such as free slip).  The interpolation is 
  !  performed with basis functions that are RBFs plus a constant.
  ! Input:
  !  pointCount - the number of "source" points and functionValues supplied
  !  sourcePoints - the location of the "source" points in the 3D space where the values of
  !    the function are known.  The sourcePoints are projected into the plane given by
  !    planeBasisVectors
  !  isTangentToInterface - a logical array indicating which sourcePoints/unitVectors are
  !    tangent to the interface where the boundary condition will be applied.  A Neumann
  !    boundary condition will be applied at these points in these directions.
  !  normalVectorIndex - where isTangentToInterface == .true., the index into unitVectors that
  !    gives the normal vector at the same sourcePoint.  This information is needed to compute
  !    the Neumann boundary condition at this point.
  !  unitVectors - the unit vectors associated with each of the sourcePoints.  Interpolation
  !    is performed by supplying the value of the vector function dotted into each of these unit
  !    vectors.  If multiple unit vectors are supplied at the same sourcePoint, they *must* be
  !    orthogonal for the interpolation to succeed.  A normal vector and two tangential vectors
  !    are needed at each interface point in order to satisfy the Dirichlet normal boundary
  !    condition and the Neumann tangential boundary conditions at these points. The unitVectors
  !    are projected into the plane given by planeBasisVectors
  !  destinationPoint - the point where the interpolation will be performed.  The destinationPoint
  !    is projected into the plane given by planeBasisVectors
  !  alpha - a constant that give the characteristic length scale of the RBFs,
  !    should be on the order of the distance between points
  !  planeBasisVectors - the basis fectors for the plane where interpolation is performed.
  !    All points are projected into this plane. 
  ! Output:
  !  coefficients - the coefficients used to interpolate a function with Dirichlet
  !    boundary conditions to the specified destinationPoint
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine rbfInterp_func_3DPlane_vec_const_tanNeu_compCoeffs(&
    pointCount, sourcePoints, isTangentToInterface, normalVectorIndex, unitVectors, &
    destinationPoint, alpha, planeBasisVectors, coefficients)

    integer, intent(in) :: pointCount
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints
    logical, dimension(pointCount), intent(in) :: isTangentToInterface
    integer, dimension(pointCount), intent(in) :: normalVectorIndex
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: unitVectors
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint
    real(kind=RKIND), intent(in) :: alpha
    real(kind=RKIND), dimension(2,3), intent(in) :: planeBasisVectors
    real(kind=RKIND), dimension(pointCount, 3), intent(out) :: coefficients

    integer :: i, j
    integer :: matrixSize

    real(kind=RKIND), dimension(pointCount,2) :: planarSourcePoints
    real(kind=RKIND), dimension(pointCount,2) :: planarUnitVectors
    real(kind=RKIND), dimension(2) :: planarDestinationPoint

    real(kind=RKIND), dimension(:,:), pointer :: matrix, matrixCopy
    real(kind=RKIND), dimension(:,:), pointer :: rhs, coeffs
    integer, dimension(:), pointer :: pivotIndices

    matrixSize = pointCount+2 ! space for constant vector in plane

    allocate(matrix(matrixSize,matrixSize))  
    allocate(matrixCopy(matrixSize,matrixSize))  
    allocate(rhs(matrixSize,2))  
    allocate(coeffs(matrixSize,2))  
    allocate(pivotIndices(matrixSize)) 

    matrix = 0.0
    rhs = 0.0
    coeffs = 0.0

    do i = 1, pointCount
      planarSourcePoints(i,1) = sum(sourcePoints(i,:)*planeBasisVectors(1,:)) 
      planarSourcePoints(i,2) = sum(sourcePoints(i,:)*planeBasisVectors(2,:)) 
      planarUnitVectors(i,1) = sum(unitVectors(i,:)*planeBasisVectors(1,:)) 
      planarUnitVectors(i,2) = sum(unitVectors(i,:)*planeBasisVectors(2,:)) 
    end do
    planarDestinationPoint(1) = sum(destinationPoint*planeBasisVectors(1,:)) 
    planarDestinationPoint(2) = sum(destinationPoint*planeBasisVectors(2,:)) 
    call setUpVectorFreeSlipRBFMatrixAndRHS(pointCount, 2, &
      planarSourcePoints, isTangentToInterface, normalVectorIndex, planarUnitVectors, &
      planarDestinationPoint, alpha, matrix(1:pointCount,1:pointCount), rhs(1:pointCount,:))

    do i = 1, pointCount
      matrix(pointCount+1,i) = sum(unitVectors(i,:)*planeBasisVectors(1,:))
      matrix(pointCount+2,i) = sum(unitVectors(i,:)*planeBasisVectors(2,:))
      if(.not. isTangentToInterface(i)) then
        matrix(i,pointCount+1:pointCount+2) = matrix(pointCount+1:pointCount+2,i)
      end if
    end do
    do i = 1,2 
      rhs(pointCount+i,i) = 1.0 ! the unit vector in the ith direction
    end do

    ! solve each linear system
    matrixCopy = matrix
    call LEGS(matrix, matrixSize, rhs(:,1), coeffs(:,1), pivotIndices)
    call LEGS(matrixCopy, matrixSize, rhs(:,2), coeffs(:,2), pivotIndices)

    coefficients(:,1) = planeBasisVectors(1,1)*coeffs(1:pointCount,1) &
      + planeBasisVectors(2,1)*coeffs(1:pointCount,2) 
    coefficients(:,2) = planeBasisVectors(1,2)*coeffs(1:pointCount,1) &
      + planeBasisVectors(2,2)*coeffs(1:pointCount,2) 
    coefficients(:,3) = planeBasisVectors(1,3)*coeffs(1:pointCount,1) &
      + planeBasisVectors(2,3)*coeffs(1:pointCount,2) 

    deallocate(matrix)  
    deallocate(matrixCopy)  
    deallocate(rhs)  
    deallocate(coeffs)  
    deallocate(pivotIndices)

   end subroutine rbfInterp_func_3DPlane_vec_const_tanNeu_compCoeffs 


!!!!!!!!!!!!!!!!!!!!!
! private subroutines
!!!!!!!!!!!!!!!!!!!!!

  function evaluateRBF(rSquared) result(rbfValue)
    real(kind=RKIND), intent(in) :: rSquared
    real(kind=RKIND) :: rbfValue

    ! inverse multiquadratic
    rbfValue = 1/sqrt(1 + rSquared)

  end function evaluateRBF

  subroutine evaluateRBFAndDeriv(rSquared, rbfValue, rbfDerivOverR)
    real(kind=RKIND), intent(in) :: rSquared
    real(kind=RKIND), intent(out) :: rbfValue, rbfDerivOverR

    ! inverse multiquadratic
    rbfValue = 1/sqrt(1 + rSquared)
    rbfDerivOverR = -rbfValue**3

  end subroutine evaluateRBFAndDeriv

  subroutine evaluateRBFAndDerivs(rSquared, rbfValue, rbfDerivOverR, rbfSecondDeriv)
    real(kind=RKIND), intent(in) :: rSquared
    real(kind=RKIND), intent(out) :: rbfValue, rbfDerivOverR, rbfSecondDeriv

    ! inverse multiquadratic
    rbfValue = 1/sqrt(1 + rSquared)
    rbfDerivOverR = -rbfValue**3
    rbfSecondDeriv = (2*rSquared-1)*rbfValue**5

  end subroutine evaluateRBFAndDerivs

  subroutine setUpScalarRBFDirichletMatrixAndRHS(pointCount, sourcePoints, destinationPoint, &
    alpha, dirichletMatrix, rhs)

    integer, intent(in) :: pointCount
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint
    real(kind=RKIND), intent(in) :: alpha
    real(kind=RKIND), dimension(pointCount,pointCount), intent(out) :: dirichletMatrix
    real(kind=RKIND), dimension(pointCount), intent(out) :: rhs

    integer :: i, j

    real(kind=RKIND) :: rSquared, rbfValue

    do j = 1, pointCount
      do i = 1, pointCount
        rSquared = sum((sourcePoints(i,:)-sourcePoints(j,:))**2)/alpha**2
        rbfValue = evaluateRBF(rSquared)
        dirichletMatrix(i,j) = rbfValue
      end do
    end do

    do j = 1, pointCount
      rSquared = sum((destinationPoint-sourcePoints(j,:))**2)/alpha**2
      rhs(j) = evaluateRBF(rSquared)
    end do

  end subroutine setUpScalarRBFDirichletMatrixAndRHS

  subroutine setUpScalarRBFMatrixAndRHS(pointCount, &
    sourcePoints, isInterface, interfaceNormals, destinationPoint, &
    alpha, dirichletMatrix, neumannMatrix, rhs)

    integer, intent(in) :: pointCount
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints
    logical, dimension(pointCount), intent(in) :: isInterface
    real(kind=RKIND), dimension(pointCount,3), intent(in) :: interfaceNormals
    real(kind=RKIND), dimension(3), intent(in) :: destinationPoint
    real(kind=RKIND), intent(in) :: alpha
    real(kind=RKIND), dimension(pointCount,pointCount), intent(out) :: &
      dirichletMatrix, neumannMatrix
    real(kind=RKIND), dimension(pointCount), intent(out) :: rhs

    integer :: i, j

    real(kind=RKIND) :: rSquared, rbfValue, rbfDerivOverR, normalDotX

    do j = 1, pointCount
      if(isInterface(j)) then
        do i = 1, pointCount 
          rSquared = sum((sourcePoints(i,:)-sourcePoints(j,:))**2)/alpha**2
          normalDotX = sum(interfaceNormals(j,:) &
            * (sourcePoints(j,:)-sourcePoints(i,:)))
          call evaluateRBFAndDeriv(rSquared, rbfValue, rbfDerivOverR)
          rbfDerivOverR = rbfDerivOverR/alpha**2
          dirichletMatrix(i,j) = rbfValue
          neumannMatrix(i,j) = rbfDerivOverR*normalDotX
        end do
      else
        do i = 1, pointCount
          rSquared = sum((sourcePoints(i,:)-sourcePoints(j,:))**2)/alpha**2
          rbfValue = evaluateRBF(rSquared)
          dirichletMatrix(i,j) = rbfValue
          neumannMatrix(i,j) = rbfValue
        end do
      end if
    end do

    do j = 1, pointCount
      rSquared = sum((destinationPoint-sourcePoints(j,:))**2)/alpha**2
      rhs(j) = evaluateRBF(rSquared)
    end do

  end subroutine setUpScalarRBFMatrixAndRHS

  subroutine setUpVectorDirichletRBFMatrixAndRHS(pointCount, dimensions, &
    sourcePoints, unitVectors, destinationPoint, &
    alpha, matrix, rhs)

    integer, intent(in) :: pointCount, dimensions
    real(kind=RKIND), dimension(pointCount,dimensions), intent(in) :: sourcePoints
    real(kind=RKIND), dimension(pointCount,dimensions), intent(in) :: unitVectors
    real(kind=RKIND), dimension(dimensions), intent(in) :: destinationPoint
    real(kind=RKIND), intent(in) :: alpha
    real(kind=RKIND), dimension(pointCount,pointCount), intent(out) :: matrix
    real(kind=RKIND), dimension(pointCount,dimensions), intent(out) :: rhs

    integer :: i, j

    real(kind=RKIND) :: rSquared, rbfValue, unitVectorDotProduct

    do j = 1, pointCount
      do i = j, pointCount
        rSquared = sum((sourcePoints(i,:)-sourcePoints(j,:))**2)/alpha**2
        rbfValue = evaluateRBF(rSquared)
        unitVectorDotProduct = sum(unitVectors(i,:)*unitVectors(j,:))
        matrix(i,j) = rbfValue*unitVectorDotProduct
        matrix(j,i) = matrix(i,j)
      end do
    end do

    do j = 1, pointCount
      rSquared = sum((destinationPoint-sourcePoints(j,:))**2)/alpha**2
      rhs(j,:) = evaluateRBF(rSquared)*unitVectors(j,:)
    end do

  end subroutine setUpVectorDirichletRBFMatrixAndRHS

  subroutine setUpVectorFreeSlipRBFMatrixAndRHS(pointCount, dimensions, &
    sourcePoints, isTangentToInterface, normalVectorIndex, unitVectors, destinationPoint, &
    alpha, matrix, rhs)

    integer, intent(in) :: pointCount, dimensions
    real(kind=RKIND), dimension(pointCount,dimensions), intent(in) :: sourcePoints
    logical, dimension(pointCount), intent(in) :: isTangentToInterface
    integer, dimension(pointCount), intent(in) :: normalVectorIndex
    real(kind=RKIND), dimension(pointCount,dimensions), intent(in) :: unitVectors
    real(kind=RKIND), dimension(dimensions), intent(in) :: destinationPoint
    real(kind=RKIND), intent(in) :: alpha
    real(kind=RKIND), dimension(pointCount,pointCount), intent(out) :: matrix
    real(kind=RKIND), dimension(pointCount,dimensions), intent(out) :: rhs

    integer :: i, j

    real(kind=RKIND) :: rSquared, rbfValue, rbfDerivOverR, normalVector(dimensions), &
      normalDotX, unitVectorDotProduct

    do j = 1, pointCount
      if(isTangentToInterface(j)) then
         do i = 1, pointCount
          rSquared = sum((sourcePoints(i,:)-sourcePoints(j,:))**2)/alpha**2
          normalVector = unitVectors(normalVectorIndex(j),:) 
          normalDotX = sum(normalVector * (sourcePoints(j,:)-sourcePoints(i,:)))
          call evaluateRBFAndDeriv(rSquared, rbfValue, rbfDerivOverR)
          rbfDerivOverR = rbfDerivOverR/alpha**2
          unitVectorDotProduct = sum(unitVectors(i,:)*unitVectors(j,:))
          matrix(i,j) = rbfDerivOverR*normalDotX*unitVectorDotProduct
        end do
      else
        do i = 1, pointCount
          rSquared = sum((sourcePoints(i,:)-sourcePoints(j,:))**2)/alpha**2
          rbfValue = evaluateRBF(rSquared)
          unitVectorDotProduct = sum(unitVectors(i,:)*unitVectors(j,:))
          matrix(i,j) = rbfValue*unitVectorDotProduct
        end do
      end if
    end do

    do j = 1, pointCount
      rSquared = sum((destinationPoint-sourcePoints(j,:))**2)/alpha**2
      rhs(j,:) = evaluateRBF(rSquared)*unitVectors(j,:)
    end do

  end subroutine setUpVectorFreeSlipRBFMatrixAndRHS

  subroutine unit_vec_in_R3(xin)
    implicit none
    real (kind=RKIND), intent(inout) :: xin(3)
    real (kind=RKIND) :: mag
    mag = sqrt(xin(1)**2+xin(2)**2+xin(3)**2)
    xin(:) = xin(:) / mag
  end subroutine unit_vec_in_R3

  subroutine cross_product_in_R3(p_1,p_2,p_out)
    real (kind=RKIND), intent(in) :: p_1 (3), p_2 (3)
    real (kind=RKIND), intent(out) :: p_out (3)

    p_out(1) = p_1(2)*p_2(3)-p_1(3)*p_2(2)
    p_out(2) = p_1(3)*p_2(1)-p_1(1)*p_2(3)
    p_out(3) = p_1(1)*p_2(2)-p_1(2)*p_2(1)
  end subroutine cross_product_in_R3

! Updated 10/24/2001.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!   Program 4.3   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!PROGRAM EX43
!
!
! An example of solving linear equation set A(N,N)*X(N) = B(N)
! with the partial-pivoting Gaussian elimination scheme.  The
! numerical values are for the Wheatstone bridge example discussed
! in Section 4.1 in the book with all resistances being 100 ohms
! and the voltage 200 volts.
!
!  IMPLICIT NONE
!  INTEGER, PARAMETER :: N=3
!  INTEGER :: I,J
!  INTEGER, DIMENSION (N) :: INDX
!  REAL, DIMENSION (N) :: X,B
!  REAL, DIMENSION (N,N) :: A
!  DATA B /200.0,0.0,0.0/, &
!       ((A(I,J), J=1,N),I=1,N) /100.0,100.0,100.0,-100.0, &
!                         300.0,-100.0,-100.0,-100.0, 300.0/
!
!  CALL LEGS (A,N,B,X,INDX)
!
!  WRITE (6, "(F16.8)") (X(I), I=1,N)
!END PROGRAM EX43


SUBROUTINE LEGS (A,N,B,X,INDX)
!
! Subroutine to solve the equation A(N,N)*X(N) = B(N) with the
! partial-pivoting Gaussian elimination scheme.
! Copyright (c) Tao Pang 2001.
!
  IMPLICIT NONE
  integer, INTENT (IN) :: N
  integer :: I,J
  integer, INTENT (OUT), DIMENSION (N) :: INDX
  real(kind=RKIND), INTENT (INOUT), DIMENSION (N,N) :: A
  real(kind=RKIND), INTENT (INOUT), DIMENSION (N) :: B
  real(kind=RKIND), INTENT (OUT), DIMENSION (N) :: X
!
  CALL ELGS (A,N,INDX)
!
  DO I = 1, N-1
    DO J = I+1, N
      B(INDX(J)) = B(INDX(J))-A(INDX(J),I)*B(INDX(I))
    END DO
  END DO
!
  X(N) = B(INDX(N))/A(INDX(N),N)
  DO I = N-1, 1, -1
    X(I) = B(INDX(I))
    DO J = I+1, N
      X(I) = X(I)-A(INDX(I),J)*X(J)
    END DO
    X(I) =  X(I)/A(INDX(I),I)
  END DO
!
END SUBROUTINE LEGS
!



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
  integer, INTENT (IN) :: N
  integer :: I,J,K
  integer, INTENT (OUT), DIMENSION (N) :: INDX
  real(kind=RKIND), INTENT (INOUT), DIMENSION (N,N):: A
  real(kind=RKIND), INTENT (OUT), DIMENSION (N,N):: X
  real(kind=RKIND), DIMENSION (N,N) :: B
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
  integer, INTENT (IN) :: N
  integer :: I,J,K,ITMP
  integer, INTENT (OUT), DIMENSION (N) :: INDX
  real(kind=RKIND) :: C1,PI,PI1,PJ
  real(kind=RKIND), INTENT (INOUT), DIMENSION (N,N) :: A
  real(kind=RKIND), DIMENSION (N) :: C
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

end module RBF_interpolation

