module SubRoutineModule
  use TypeModule
  use PrecMod
  implicit none
  
  contains

  subroutine TimeStep(ParamAM, DynVarAM, Neighbours)
    implicit none
    !Input
    type(ParamType), intent(in) :: ParamAM
    type(DynamicVars), intent(inout) :: DynVarAM
    type(NeighbourType), dimension(:), intent(in) :: Neighbours
    
    !Allocate
    integer(wpi) :: iParticle

    !Loop over particles, this loop has been deliberately made to be very easy to parallelize
    do iParticle = 1,ParamAM%NumPart
      call EulerStep(ParamAM, DynVarAM, Neighbours, iParticle)
    end do

  end subroutine TimeStep


  subroutine EulerStep(ParamAM, DynVarAM, Neighbours, iParticle)
    implicit none
    !Input
    type(ParamType), intent(in) :: ParamAM
    type(DynamicVars), intent(inout) :: DynVarAM
    type(NeighbourType), dimension(:), intent(in) :: Neighbours
    integer(wpi), intent(in) :: iParticle
    
    !Allocate
    integer(wpi) :: iNeighbour
    real(wpf), dimension(2) :: TotalForce, PolVec
    real(wpf) :: TotalTorque

    TotalForce = (/0._wpf,0._wpf/)
    PolVec = (/cos(DynVarAM%PolAng(iParticle,DynVarAM%OldIndex)),sin(DynVarAM%PolAng(iParticle,DynVarAM%OldIndex))/)
    TotalTorque = 0_wpf

    do iNeighbour = 1, Neighbours(iParticle)%NumNeighbours
      
      !Adding all elastic contributions
      TotalForce = TotalForce + LinearElasticForce(ParamAM,DynVarAM,iParticle,Neighbours(iParticle),iNeighbour)


    end do

    !Calculating the torque using the TotalForce, as it per now only contains elastic force, and is thus the total elastic force
    TotalTorque = TotalTorque + PolarizationTorque(ParamAM,TotalForce,PolVec)

    !Adding on the polarization force, to the elastic energy to produce the total force and dividing by zeta to pruce the "velocity"
    TotalForce = (TotalForce + PolarizationForce(ParamAM,PolVec))/ParamAM%zeta

    !TODO: Setting some kind of boundary thing.

    !Setting the new timestep using the previous, with an added velocity-Euler-step
    DynVarAM%Coords(:,iParticle,DynVarAM%NewIndex) = DynVarAM%Coords(:,iParticle,DynVarAM%OldIndex) + TotalForce*ParamAM%deltaT
    DynVarAM%PolAng(iParticle,DynVarAM%NewIndex) = DynVarAM%PolAng(iParticle,DynVarAM%OldIndex) + TotalTorque*ParamAM%deltaT



  end subroutine EulerStep

  function LinearElasticForce(ParamAM, DynVarAM, iPart, Neighbour, iNeigh) result(ElForce)
    implicit none
    !Input
    type(ParamType), intent(in) :: ParamAM
    type(DynamicVars), intent(in) :: DynVarAM
    integer(wpi), intent(in) :: iPart, iNeigh
    type(NeighbourType), intent(in) :: Neighbour
    !Allocate
    real(wpf), dimension(2) :: DeltaCoords
    real(wpf) :: length
    !Output
    real(wpf), dimension(2) :: ElForce

    !Calculating the difference in coordinates between the particle and the neighbour
    DeltaCoords = DynVarAM%Coords(:,iPart,DynVarAM%OldIndex) &
    - DynVarAM%Coords(:,Neighbour%SpecificNeighbours(iNeigh),DynVarAM%OldIndex)

    !Calculating the length
    length = EuclideanNormVec(DeltaCoords)

    !Force is calculated as -k(l-l_0)\vec(d), where d is in the direction of the displacement, as a unitary vector
    Elforce = -ParamAM%k*(length-Neighbour%EquilLength(iNeigh))*DeltaCoords/length

  end function
  
  function PolarizationForce(ParamAM, PolVec) result(PolForce)
    implicit none
    !Input
    type(ParamType), intent(in) :: ParamAM
    real(wpf), dimension(2), intent(in) :: PolVec
    !Output
    real(wpf), dimension(2) :: PolForce

    PolForce = ParamAM%Fa*PolVec

  end function 

  function PolarizationTorque(ParamAM,ElasticForce,PolVec) result(PolTorque)
    implicit none
    !Input
    type(ParamType), intent(in) :: ParamAM
    real(wpf), dimension(2), intent(in) :: ElasticForce, PolVec
    !Allocate
    real(wpf), dimension(2) :: NormPolVec
    !Output
    real(wpf) :: PolTorque

    NormPolVec = (/-PolVec(2),PolVec(1)/)

    PolTorque = ParamAM%xi * dot_product(ElasticForce,NormPolVec)

  end function

  subroutine SwitchIndex(DynVarAM)
    implicit none
    !Input/Output
    type(DynamicVars), intent(inout) :: DynVarAM

    !Switching, or flip-floping, such that the index that previously was 1 becomes 2, while 2 becomes 1
    !The new index becomes the old index - to be read in the next timestep, and old becomes new index - to we written over
    DynVarAM%OldIndex = 3_wpi - DynVarAM%OldIndex
    DynVarAM%NewIndex = 3_wpi - DynVarAM%NewIndex

  end subroutine SwitchIndex

  
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
