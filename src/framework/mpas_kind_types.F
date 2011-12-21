module mpas_kind_types

#ifdef SINGLE_PRECISION
   integer, parameter :: RKIND  = selected_real_kind(6)
#else
   integer, parameter :: RKIND  = selected_real_kind(12)
#endif

   contains

   subroutine dummy_kinds()

   end subroutine dummy_kinds

end module mpas_kind_types
