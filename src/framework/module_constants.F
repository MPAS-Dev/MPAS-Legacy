module constants

   real (kind=RKIND), parameter :: pii     = 3.141592653589793   
   real (kind=RKIND), parameter :: a       = 6371229.0
   real (kind=RKIND), parameter :: omega   = 7.29212e-5
   real (kind=RKIND), parameter :: gravity = 9.80616
   real (kind=RKIND), parameter :: rgas = 287.
   real (kind=RKIND), parameter :: cp = 1003.
   real (kind=RKIND), parameter :: cv = 716.  ! cp - rgas
   real (kind=RKIND), parameter :: cvpm = -.71385842 ! -cv/cp
   real (kind=RKIND), parameter :: prandtl = 1.0


   contains

   subroutine dummy()

   end subroutine dummy

end module constants
