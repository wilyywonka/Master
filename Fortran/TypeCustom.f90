module TypeModule
  use PrecMod
  implicit none
  
  type :: ParamType

    integer(wpi) :: NumPart, MaxNeighbour

    real(wpf) :: zeta, xi, Fa, k, deltaT, endT

  end type ParamType

  type :: DynamicVars

    real(wpf), dimension(:,:,:), allocatable :: Coords, PolVec
    real(wpf), dimension(:,:), allocatable :: PolAng
    integer(wpi) :: NewIndex, OldIndex

  end type DynamicVars

  type :: NeighbourType

    integer(wpi), dimension(:), allocatable :: SpecificNeighbours
    real(wpf), dimension(:), allocatable :: EquilLength
    integer(wpi) :: NumNeighbours

  end type NeighbourType

  type :: InitValues

    integer(wpi), dimension(:,:), allocatable :: NeighbourMatrix
    real(wpf), dimension(:,:), allocatable :: Coords, PolVec
    real(wpf), dimension(:), allocatable :: PolAng

  end type InitValues

  contains

  function EucledianNormVec(Coords) result(length)
    implicit none
    !Input
    real(wpf), dimension(2), intent(in) :: Coords

    !Output
    real(wpf) :: length

    length = sqrt(Coords(1)**2 + Coords(2)**2)

  end function

  subroutine Initialize(ParamAM, DynVarAM, InitSetup)
    implicit none
    type(ParamType), intent(inout) :: ParamAM
    type(DynamicVars), intent(inout) :: DynVarAM
    type(InitValues), intent(inout) :: InitSetup

    ! READ NAMELIST

    DynVarAM%OldIndex = 1_wpi
    DynVarAM%NewIndex = 2_wpi

    
    ! Saved as (2,NumPart), as it is usual to need both x and y when considering lengths
    allocate(InitSetup%NeighbourMatrix(ParamAM%MaxNeighbour, ParamAM%NumPart))
    allocate(InitSetup%Coords(2,ParamAM%NumPart))
    allocate(InitSetup%PolVec(2,ParamAM%NumPart))
    allocate(InitSetup%PolAng(ParamAM%NumPart))

    ! READ INITFILES

    ! Saved as (2,NumPart,2) = (x/y,iPart,new/old), as it is usual to need both x and y when considering lengths
    allocate(DynVarAM%Coords(2,ParamAM%NumPart,2)) 
    allocate(DynVarAM%PolVec(2,ParamAM%NumPart,2))
    allocate(DynVarAM%PolAng(ParamAM%NumPart,2))
    
    
    DynVarAM%Coords(:,:,DynVarAM%OldIndex) = InitSetup%Coords
    DynVarAM%PolVec(:,:,DynVarAM%OldIndex) = InitSetup%PolVec
    DynVarAM%PolAng(:,DynVarAM%OldIndex) = InitSetup%PolAng

    
  end subroutine Initialize
  
  subroutine InitializeNeighbours(ParamAM,InitSetup,Neighbours)
    !In/Out
    type(ParamType), intent(in) :: ParamAM
    type(InitValues), intent(in) :: InitSetup
    type(NeighbourType), dimension(:), allocatable, intent(inout) :: Neighbours
    

    !Allocate
    integer(wpi) :: iParticle, iNeighbour, SpecNumNeighbours

    allocate(Neighbours(ParamAM%NumPart))

    do iParticle = 1,ParamAM%NumPart
      SpecNumNeighbours = 0
      do iNeighbour = 1,ParamAM%MaxNeighbour
        if ( InitSetup%NeighbourMatrix(iNeighbour,iParticle) /= 0 ) then
          SpecNumNeighbours = SpecNumNeighbours + 1_wpi
        end if
      end do
      Neighbours(iParticle)%NumNeighbours = SpecNumNeighbours
      allocate(Neighbours(iParticle)%SpecificNeighbours(SpecNumNeighbours))
      allocate(Neighbours(iParticle)%EquilLength(SpecNumNeighbours))
      do iNeighbour = 1, SpecNumNeighbours
        Neighbours(iParticle)%SpecificNeighbours(iNeighbour) = InitSetup%NeighbourMatrix(iNeighbour,iParticle)
        Neighbours(iParticle)%EquilLength(iNeighbour) = EucledianNormVec(InitSetup%Coords(:,iParticle) &
        - InitSetup%Coords(:,iNeighbour))
      end do
    end do

  end subroutine InitializeNeighbours

  


end module TypeModule