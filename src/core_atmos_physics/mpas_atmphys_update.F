!=============================================================================================
 module mpas_atmphys_update
 use mpas_configure
 use mpas_grid_types

 use mpas_atmphys_driver_convection_deep

 implicit none
 private
 public:: physics_update

 contains
 
!=============================================================================================
 subroutine physics_update(domain,dt)
!=============================================================================================

!input arguments:
 type(domain_type),intent(in):: domain
 real(kind=RKIND),intent(in) :: dt 

!local variables:
 type(block_type),pointer:: block
 integer:: i,j,k

!---------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- begin physics_update:'
 block => domain % blocklist
 do while(associated(block))

    !parameterization of convection: update accumulated precipitation.
    call update_convection_deep(dt,block%mesh,block%diag_physics)

    block => block % next
 end do
 write(0,*) '--- end physics_update:'

 end subroutine physics_update

!=============================================================================================
 end module mpas_atmphys_update
!=============================================================================================
