module SubRoutineModule
  use TypeModule
  use PrecMod
  implicit none
  
  contains

  subroutine TimeStep(ParamAM, DynVarAM)
    implicit none
    type(ParamType), intent(in) :: ParamAM
    type(DynamicVars), intent(inout) :: DynVarAM


  end subroutine TimeStep


  subroutine EulerStep(ParamAM, DynVarAM, IndexPart)
    implicit none
    type(ParamType), intent(inout) :: ParamAM
    type(DynamicVars), intent(inout) :: DynVarAM
    integer(wpi), intent(in) :: IndexPart

    !Allocate
    integer(wpi) :: iNeighbours, Currentneighbour
    real(wpf), dimension(2) :: TotalForce

    TotalForce = (0._wpf,0._wpf)

    do iNeighbours = 1, ParamAM%MaxNeighbour
      
      TotalForce = TotalForce + LinearElasticForce(ParamAM,DynVarAM,Neighbours(iNeighbours))

    end do

  end subroutine EulerStep

  function LinearElasticForce(ParamAM, DynVarAM, Neighbours) result(ElForce)
    implicit none
    !Input
    type(ParamType), intent(in) :: ParamAM
    type(DynamicVars), intent(in) :: DynVarAM
    type(NeighbourType), dimension(:), intent(in) :: Neighbours

    !Output
    real(wpf), intent(out) :: ElForce

    

  end function
  
  function PolarizationForce(ParamAM) result(PolForce)
    implicit none
    !Input
    type(ParamType), intent(in) :: ParamAM

    !Output
    real(wpf), intent(out) :: PolForce

  end function 

  function PolarizationTorque(ParamAM) result(PolTorque)
    implicit none
    !Input
    type(ParamType), intent(in) ::ParamAM

    !Output
    real(wpf), intent(out) :: PolTorque


  end function 

  function EucledianNormVec(Coords) result(length)
    implicit none
    !Input
    real(wpf), dimension(2), intent(in) :: Coords

    !Output
    real(wpf), intent(out) :: length

    length = sqrt(Coords(1)**2 + Coords(2)**2)

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
