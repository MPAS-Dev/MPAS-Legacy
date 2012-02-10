module ocn_gm

   use mpas_grid_types
   use mpas_configure
   use mpas_timer
   
   implicit none
   private 
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_gm_compute_uBolus

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

contains

   subroutine ocn_gm_compute_uBolus(s, grid)
      implicit none
      type(state_type), intent(inout)        :: s
      type(mesh_type), intent(in)            :: grid

      real(kind=RKIND), dimension(:,:), pointer :: uBolusGM, hEddyFlux, h_edge

      integer, dimension(:), pointer   :: maxLevelEdgeTop
      integer                          :: k, iEdge, nEdges

      uBolusGM         => s % uBolusGM % array
      h_edge         => s % h_edge % array
      hEddyFlux      => s % hEddyFlux % array

      maxLevelEdgeTop => grid % maxLevelEdgeTop % array

      nEdges = grid % nEdges

      call ocn_gm_compute_hEddyFlux(s, grid)

      if (config_vert_grid_type .EQ. 'z-level') then
         ! Nothing for now
         uBolusGM(:,:) = 0.0

      else if (config_vert_grid_type .EQ. 'isopycnal') then
         do iEdge = 1, nEdges
            do k = 1, maxLevelEdgeTop(iEdge)
               uBolusGM(k,iEdge) = hEddyFlux(k,iEdge)/h_edge(k,iEdge)
            end do
         end do
      end if

   end subroutine ocn_gm_compute_uBolus


   subroutine ocn_gm_compute_hEddyFlux(s, grid)
      implicit none
      type(state_type), intent(inout)     :: s
      type(mesh_type), intent(in)         :: grid

      real(kind=RKIND), dimension(:,:), pointer  :: hEddyFlux, h
      real(kind=RKIND), dimension(:), pointer    :: dcEdge
      integer, dimension(:,:), pointer           :: cellsOnEdge
      integer, dimension(:), pointer             :: maxLevelEdgeTop
      integer                                    :: k, cell1, cell2, iEdge, nEdges

      hEddyFlux      => s % hEddyFlux % array
      h              => s % h % array

      dcEdge         => grid % dcEdge % array
      cellsOnEdge    => grid % cellsOnEdge % array
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array

      nEdges = grid % nEdges

      hEddyFlux(:,:) = 0.0

      if (config_vert_grid_type .EQ. 'z-level') then

         !Nothing for now

      else if (config_vert_grid_type .EQ. 'isopycnal') then
            do iEdge = 1,nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               do k=1,maxLevelEdgeTop(iEdge)
                  hEddyFlux(k,iEdge) = -config_h_kappa * (h(k,cell2) - h(k,cell1)) / dcEdge(iEdge)
               end do
            end do
       end if
                  
   end subroutine ocn_gm_compute_hEddyFlux



   subroutine ocn_get_h_kappa(s, grid)

      type (state_type), intent(inout) :: s
      type (mesh_type), intent(in) :: grid

      real(kind=RKIND), dimension(:,:), pointer    :: h_kappa


      h_kappa  => s % h_kappa % array

      h_kappa(:,:) = config_h_kappa


   end subroutine ocn_get_h_kappa


   subroutine ocn_get_h_kappa_q(s, grid)

      type (state_type), intent(inout) :: s
      type (mesh_type), intent(in) :: grid

      real(kind=RKIND), dimension(:,:), pointer    :: h_kappa_q


      h_kappa_q  => s % h_kappa_q % array

      h_kappa_q(:,:) = config_h_kappa_q


   end subroutine ocn_get_h_kappa_q

end module ocn_gm
