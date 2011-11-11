!=============================================================================================
 module mpas_atmphys_todynamics
 use mpas_configure
 use mpas_grid_types

 use mpas_atmphys_constants, only: R_d,R_v

 implicit none
 private
 public:: physics_addtend

 contains
 
!=============================================================================================
 subroutine physics_addtend(dminfo,cellsToSend,cellsToRecv,mesh,state,diag,tend, &
                            tend_physics,mass,mass_edge)
!=============================================================================================

!input variables:
!----------------
 type(dm_info), intent(in):: dminfo
 type(mesh_type),intent(in):: mesh
 type(exchange_list),pointer:: cellsToSend,cellsToRecv

 type(state_type),intent(in):: state
 type(diag_type),intent(in):: diag
 type(tend_physics_type),intent(in):: tend_physics 
 real(kind=RKIND),dimension(:,:):: mass
 real(kind=RKIND),dimension(:,:):: mass_edge

!inout variables:
!----------------
 type(tend_type),intent(inout):: tend

!local variables:
!----------------
 integer:: i,k,nCells,nCellsSolve,nEdges,nEdgesSolve,nVertLevels

 real(kind=RKIND),dimension(:,:),pointer:: theta,theta_m,qv
 real(kind=RKIND),dimension(:,:),pointer:: rthblten,rqvblten,rqcblten, &
                                           rqiblten,rublten,rvblten
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten, &
                                           rqrcuten,rqicuten,rqscuten
 real(kind=RKIND),dimension(:,:),pointer:: rthratenlw,rthratensw                                    
 
 real(kind=RKIND),dimension(:,:),pointer  :: tend_theta,tend_u
 real(kind=RKIND),dimension(:,:,:),pointer:: tend_scalars

 real(kind=RKIND):: tem
 real(kind=RKIND),dimension(:,:),allocatable:: rublten_Edge

!=============================================================================================
!write(0,*)
!write(0,*) '--- enter subroutine physics_add_tend:'
 
 nCells      = mesh % nCells
 nEdges      = mesh % nEdges
 nCellsSolve = mesh % nCellsSolve
 nEdgesSolve = mesh % nEdgesSolve
 nVertLevels = mesh % nVertLevels

 theta   => diag % theta % array
 theta_m => state % theta_m % array
 qv      => state % scalars % array(state%index_qv,:,:)

 rublten    => tend_physics % rublten    % array
 rvblten    => tend_physics % rvblten    % array
 rthblten   => tend_physics % rthblten   % array
 rqvblten   => tend_physics % rqvblten   % array
 rqcblten   => tend_physics % rqcblten   % array
 rqiblten   => tend_physics % rqiblten   % array

 rthcuten   => tend_physics % rthcuten   % array
 rqvcuten   => tend_physics % rqvcuten   % array
 rqccuten   => tend_physics % rqccuten   % array
 rqrcuten   => tend_physics % rqrcuten   % array
 rqicuten   => tend_physics % rqicuten   % array
 rqscuten   => tend_physics % rqscuten   % array

 rthratenlw => tend_physics % rthratenlw % array
 rthratensw => tend_physics % rthratensw % array

 tend_u       => tend % u % array
 tend_theta   => tend % theta_m % array
 tend_scalars => tend % scalars % array

 tend_scalars = 0.

!add coupled tendencies due to PBL processes:
 if(config_pbl_scheme .ne. 'off') then
    allocate(rublten_Edge(nVertLevels,nEdges))
    rublten_Edge(:,:) = 0.
    call tend_toEdges(dminfo,CellsToSend,CellsToRecv,mesh,rublten,rvblten,rublten_Edge)
    do i = 1, nEdgesSolve
    do k  = 1, nVertLevels
       tend_u(k,i)=tend_u(k,i)+rublten_Edge(k,i)*mass_edge(k,i)
    enddo
    enddo

    deallocate(rublten_Edge)

    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_theta(k,i)=tend_theta(k,i)+rthblten(k,i)*mass(k,i)
       tend_scalars(tend%index_qv,k,i)=tend_scalars(tend%index_qv,k,i)+rqvblten(k,i)*mass(k,i)
       tend_scalars(tend%index_qc,k,i)=tend_scalars(tend%index_qc,k,i)+rqcblten(k,i)*mass(k,i)
       tend_scalars(tend%index_qi,k,i)=tend_scalars(tend%index_qi,k,i)+rqiblten(k,i)*mass(k,i)
    enddo
    enddo
 endif

!add coupled tendencies due to convection:
 if(config_conv_deep_scheme .ne. 'off') then
    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_theta(k,i)=tend_theta(k,i)+rthcuten(k,i)*mass(k,i)
       tend_scalars(tend%index_qv,k,i)=tend_scalars(tend%index_qv,k,i)+rqvcuten(k,i)*mass(k,i)
       tend_scalars(tend%index_qc,k,i)=tend_scalars(tend%index_qc,k,i)+rqccuten(k,i)*mass(k,i)
       tend_scalars(tend%index_qr,k,i)=tend_scalars(tend%index_qr,k,i)+rqrcuten(k,i)*mass(k,i)
       tend_scalars(tend%index_qi,k,i)=tend_scalars(tend%index_qi,k,i)+rqicuten(k,i)*mass(k,i)
       tend_scalars(tend%index_qs,k,i)=tend_scalars(tend%index_qs,k,i)+rqscuten(k,i)*mass(k,i)    
    enddo
    enddo
 endif

!add coupled tendencies due to longwave radiation:
 if(config_radt_lw_scheme .ne. 'off') then
    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_theta(k,i)=tend_theta(k,i)+rthratenlw(k,i)*mass(k,i)
    enddo
    enddo
 endif
    
!add coupled tendencies due to shortwave radiation:
 if(config_radt_sw_scheme .ne. 'off') then
    do i = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_theta(k,i)=tend_theta(k,i)+rthratensw(k,i)*mass(k,i)
    enddo
    enddo
 endif

!if non-hydrostatic core, convert the tendency for the potential temperature to a
!tendency for the modified potential temperature:
#ifdef non_hydrostatic_core
    do i = 1, nCellsSolve
    do k = 1, nVertLevels

       theta(k,i) = theta_m(k,i) / (1. + R_v/R_d * qv(k,i))
       tend_theta(k,i) = (1. + R_v/R_d * qv(k,i)) * tend_theta(k,i) &
                       + R_v/R_d * theta(k,i) * tend_scalars(tend%index_qv,k,i)
    enddo
    enddo 
#endif

!write(0,*) 'max PBL tendencies:'
!write(0,*) 'max rthblten=',maxval(rthblten(:,:))
!write(0,*) 'max rqvblten=',maxval(rqvblten(:,:))
!write(0,*) 'max rqcblten=',maxval(rqcblten(:,:))
!write(0,*) 'max rqiblten=',maxval(rqiblten(:,:))
!write(0,*) 'max rublten =',maxval(rublten(:,:))
!write(0,*) 'max rvblten =',maxval(rvblten(:,:))
!write(0,*)
!write(0,*) 'max CU tendencies:'
!write(0,*) 'max rthcuten=',maxval(rthcuten(:,:))
!write(0,*) 'max rqvcuten=',maxval(rqvcuten(:,:))
!write(0,*) 'max rqccuten=',maxval(rqccuten(:,:))
!write(0,*) 'max rqrcuten=',maxval(rqrcuten(:,:))
!write(0,*) 'max rqicuten=',maxval(rqicuten(:,:))
!write(0,*) 'max rqscuten=',maxval(rqscuten(:,:))
!write(0,*)
!write(0,*) 'max tend_scalars:'
!write(0,*) 'max tend qv=',maxval(tend_scalars(tend%index_qv,:,:))
!write(0,*) 'max tend qc=',maxval(tend_scalars(tend%index_qc,:,:))
!write(0,*) 'max tend qr=',maxval(tend_scalars(tend%index_qr,:,:))
!write(0,*) 'max tend qi=',maxval(tend_scalars(tend%index_qi,:,:))
!write(0,*) 'max tend qs=',maxval(tend_scalars(tend%index_qs,:,:))
!write(0,*)

!formats:
 201 format(2i6,8(1x,e15.8))

 end subroutine  physics_addtend

!=============================================================================================
 subroutine tend_toEdges(dminfo,cellsToSend,cellsToRecv,mesh,Ux_tend,Uy_tend,U_tend)
!=============================================================================================

!input arguments:
!----------------
 type(dm_info),intent(in):: dminfo
 type(mesh_type),intent(in):: mesh
 type(exchange_list),intent(in),pointer:: cellsToSend,cellsToRecv

 real(kind=RKIND),intent(in),dimension(:,:):: Ux_tend,Uy_tend 

!output arguments:
!-----------------
 real(kind=RKIND),intent(out),dimension(:,:):: U_tend

!local variables:
 integer:: iCell,iEdge,k,j,nCells,nCellsSolve,nVertLevels
 integer,dimension(:),pointer  :: nEdgesOnCell
 integer,dimension(:,:),pointer:: edgesOnCell

 real(kind=RKIND),dimension(:,:),pointer:: east,north,edge_normal
 real(kind=RKIND),dimension(:,:),allocatable:: Ux_tend_halo,Uy_tend_halo
 
!---------------------------------------------------------------------------------------------

 nCells = mesh % nCells
 nCellsSolve  = mesh % nCellsSolve
 nVertLevels  = mesh % nVertLevels

 east  => mesh % east  % array
 north => mesh % north % array
 edgesOnCell  => mesh % edgesOnCell  % array
 nEdgesOnCell => mesh % nEdgesOnCell % array
 edge_normal  => mesh % edgeNormalVectors % array

 allocate(Ux_tend_halo(nVertLevels,nCells))
 allocate(Uy_tend_halo(nVertLevels,nCells))

 Ux_tend_halo(:,:) = 0.
 Uy_tend_halo(:,:) = 0.
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       Ux_tend_halo(k,iCell) = Ux_tend(k,iCell)
       Uy_tend_halo(k,iCell) = Uy_tend(k,iCell)
    enddo
 enddo

 call mpas_dmpar_exch_halo_field2d_real( &
            dminfo,Ux_tend_halo,nVertLevels,nCells,cellsToSend,cellsToRecv)
 call mpas_dmpar_exch_halo_field2d_real( &
            dminfo,Uy_tend_halo,nVertLevels,nCells,cellsToSend,cellsToRecv)

 U_tend(:,:) = 0.0
 do iCell = 1, nCells
 do j = 1, nEdgesOnCell(iCell)
    iEdge = edgesOnCell(j,iCell)
    do k = 1, nVertLevels
       U_tend(k,iEdge) = U_tend(k,iEdge) & 
                + 0.5 * Ux_tend_halo(k,iCell) * (edge_normal(1,iEdge) * east(1,iCell)   &
                                              +  edge_normal(2,iEdge) * east(2,iCell)   &
                                              +  edge_normal(3,iEdge) * east(3,iCell))  &
                + 0.5 * Uy_tend_halo(k,iCell) * (edge_normal(1,iEdge) * north(1,iCell)  &
                                              +  edge_normal(2,iEdge) * north(2,iCell)  &
                                              +  edge_normal(3,iEdge) * north(3,iCell))
    enddo
 enddo
 enddo

 deallocate(Ux_tend_halo)
 deallocate(Uy_tend_halo)
 
 end subroutine tend_toEdges

!=============================================================================================
 end module mpas_atmphys_todynamics
!=============================================================================================
