! subroutine readHDF5(filename)
!   use h5fortran
!   implicit none
!   character (len = 100), intent(in) :: filename
!   real, allocatable                 :: particleParams(:)


!   call h5f%open('h5fortran_example2.h5', action='w')
!   call h5f%n('/x', 123)
!   call h5f%close()

  
! end subroutine readHDF5


module SubRoutineModule
  use TypeModule
  implicit none
  
  contains

  
  
end module SubRoutineModule