module SubRoutineModule
  use TypeModule
  use PrecMod
  implicit none
  
  contains

  subroutine EulerStep(ParamAM, IndexPart)
    implicit none
    type(ParamType), intent(inout) :: ParamAM
    integer(wpi), intent(in) :: IndexPart

    !Allocate
    integer(wpi) :: iNeighbours, Currentneighbour
    
    do iNeighbours = 1, ParamAM%MaxNeighbour
      Currentneighbour = ParamAM%NeighbourMatrix(IndexPart,iNeighbours)

    end do

  end subroutine EulerStep

  function LinearElasticForce(ParamAM, DynVarAM, Neighbours) result(ElForce)
    implicit none
    !Input
    type(ParamType), intent(inout) :: ParamAM
    type(DynamicVars), intent(inout) :: DynVarAM
    type(NeighbourType), dimension(:), intent(inout) :: Neighbours

    !Output
    real(wpf) :: ElForce

    

  end function
  
  function PolarizationForce(ParamAM) result(PolForce)
    implicit none
    !Input
    type(ParamType), intent(inout) :: ParamAM

    !Output
    real(wpf) :: PolForce

  end function 

  function PolarizationTorque(ParamAM) result(PolTorque)
    implicit none
    !Input
    type(ParamType), intent(inout) ::ParamAM

    !Output
    real(wpf) :: PolTorque


  end function 

  
  
end module SubRoutineModule









! subroutine readHDF5(filename)
!   use h5fortran
!   implicit none
!   character (len = 100), intent(in) :: filename
!   real, allocatable                 :: particleParams(:)


!   call h5f%open('h5fortran_example2.h5', action='w')
!   call h5f%n('/x', 123)
!   call h5f%close()

  
! end subroutine readHDF5
