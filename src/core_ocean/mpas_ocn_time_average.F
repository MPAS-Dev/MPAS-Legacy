module ocn_time_average

    use mpas_grid_types

    implicit none
    save
    public

    contains 

    subroutine ocn_time_average_init(state)!{{{
        type (state_type), intent(inout) :: state

        real (kind=RKIND), pointer :: nAccumulate

        real (kind=RKIND), dimension(:), pointer :: acc_ssh, acc_sshVar
        real (kind=RKIND), dimension(:,:), pointer :: acc_uReconstructZonal, acc_uReconstructMeridional, acc_uReconstructZonalVar, acc_uReconstructMeridionalVar
        real (kind=RKIND), dimension(:,:), pointer :: acc_u, acc_uVar

        nAccumulate => state % nAccumulate % scalar

        acc_ssh => state % acc_ssh % array
        acc_sshVar => state % acc_sshVar % array
        acc_uReconstructZonal => state % acc_uReconstructZonal % array
        acc_uReconstructMeridional => state % acc_uReconstructMeridional % array
        acc_uReconstructZonalVar => state % acc_uReconstructZonalVar % array
        acc_uReconstructMeridionalVar => state % acc_uReconstructMeridionalVar % array
        acc_u => state % acc_u % array
        acc_uVar => state % acc_uVar % array

        nAccumulate = 0

        acc_ssh = 0.0
        acc_sshVar = 0.0
        acc_uReconstructZonal = 0.0
        acc_uReconstructMeridional = 0.0
        acc_uReconstructZonalVar = 0.0
        acc_uReconstructMeridionalVar = 0.0
        acc_u = 0.0
        acc_uVar = 0.0

    end subroutine ocn_time_average_init!}}}

    subroutine ocn_time_average_accumulate(state, old_state)!{{{
        type (state_type), intent(inout) :: state
        type (state_type), intent(in) :: old_state

        real (kind=RKIND), pointer :: nAccumulate, old_nAccumulate

        real (kind=RKIND), dimension(:), pointer :: ssh
        real (kind=RKIND), dimension(:,:), pointer :: uReconstructZonal, uReconstructMeridional, u

        real (kind=RKIND), dimension(:,:), pointer :: acc_u, acc_uVar
        real (kind=RKIND), dimension(:,:), pointer :: acc_uReconstructZonal, acc_uReconstructMeridional, acc_uReconstructZonalVar, acc_uReconstructMeridionalVar
        real (kind=RKIND), dimension(:), pointer :: acc_ssh, acc_sshVar

        real (kind=RKIND), dimension(:,:), pointer :: old_acc_u, old_acc_uVar
        real (kind=RKIND), dimension(:,:), pointer :: old_acc_uReconstructZonal, old_acc_uReconstructMeridional, old_acc_uReconstructZonalVar, old_acc_uReconstructMeridionalVar
        real (kind=RKIND), dimension(:), pointer :: old_acc_ssh, old_acc_sshVar

        old_nAccumulate => old_state % nAccumulate  % scalar
        nAccumulate => state % nAccumulate  % scalar

        ssh => state % ssh % array
        uReconstructZonal => state % uReconstructZonal % array
        uReconstructMeridional => state % uReconstructMeridional % array
        u => state % u % array

        acc_ssh => state % acc_ssh % array
        acc_sshVar => state % acc_sshVar % array
        acc_uReconstructZonal => state % acc_uReconstructZonal % array
        acc_uReconstructMeridional => state % acc_uReconstructMeridional % array
        acc_uReconstructZonalVar => state % acc_uReconstructZonalVar % array
        acc_uReconstructMeridionalVar => state % acc_uReconstructMeridionalVar % array
        acc_u => state % acc_u % array
        acc_uVar => state % acc_uVar % array

        old_acc_ssh => old_state % acc_ssh % array
        old_acc_sshVar => old_state % acc_sshVar % array
        old_acc_uReconstructZonal => old_state % acc_uReconstructZonal % array
        old_acc_uReconstructMeridional => old_state % acc_uReconstructMeridional % array
        old_acc_uReconstructZonalVar => old_state % acc_uReconstructZonalVar % array
        old_acc_uReconstructMeridionalVar => old_state % acc_uReconstructMeridionalVar % array
        old_acc_u => old_state % acc_u % array
        old_acc_uVar => old_state % acc_uVar % array

        acc_ssh = old_acc_ssh + ssh
        acc_sshVar = old_acc_sshVar + ssh**2
        acc_uReconstructZonal = old_acc_uReconstructZonal + uReconstructZonal
        acc_uReconstructMeridional = old_acc_uReconstructMeridional + uReconstructMeridional
        acc_uReconstructZonalVar = old_acc_uReconstructZonalVar + uReconstructZonal**2
        acc_uReconstructMeridionalVar = old_acc_uReconstructMeridionalVar + uReconstructMeridional**2
        acc_u = old_acc_u + u
        acc_uVar = old_acc_uVar + u**2

        nAccumulate = old_nAccumulate + 1
    end subroutine ocn_time_average_accumulate!}}}

    subroutine ocn_time_average_normalize(state)!{{{
        type (state_type), intent(inout) :: state

        real (kind=RKIND), pointer :: nAccumulate

        real (kind=RKIND), dimension(:), pointer :: acc_ssh, acc_sshVar
        real (kind=RKIND), dimension(:,:), pointer :: acc_uReconstructZonal, acc_uReconstructMeridional, acc_uReconstructZonalVar, acc_uReconstructMeridionalVar
        real (kind=RKIND), dimension(:,:), pointer :: acc_u, acc_uVar

        nAccumulate => state % nAccumulate  % scalar

        acc_ssh => state % acc_ssh % array
        acc_sshVar => state % acc_sshVar % array
        acc_uReconstructZonal => state % acc_uReconstructZonal % array
        acc_uReconstructMeridional => state % acc_uReconstructMeridional % array
        acc_uReconstructZonalVar => state % acc_uReconstructZonalVar % array
        acc_uReconstructMeridionalVar => state % acc_uReconstructMeridionalVar % array
        acc_u => state % acc_u % array
        acc_uVar => state % acc_uVar % array

        if(nAccumulate > 0) then
          acc_ssh = acc_ssh / nAccumulate
          acc_sshVar = acc_sshVar / nAccumulate
          acc_uReconstructZonal = acc_uReconstructZonal / nAccumulate
          acc_uReconstructMeridional = acc_uReconstructMeridional / nAccumulate
          acc_uReconstructZonalVar = acc_uReconstructZonalVar / nAccumulate
          acc_uReconstructMeridionalVar = acc_uReconstructMeridionalVar / nAccumulate
          acc_u = acc_u / nAccumulate
          acc_uVar = acc_uVar / nAccumulate
        end if
    end subroutine ocn_time_average_normalize!}}}

end module ocn_time_average
