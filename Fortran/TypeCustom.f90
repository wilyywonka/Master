module TypeModule
  use PrecMod
  implicit none
  
  type :: ParamType

    integer(wpi) :: NumPart, MaxNeighbour

    real(wpf) :: zeta, xi, Fa, k, deltaT, endT

  end type ParamType

  type :: DynamicVars

    real(wpf), dimension(:,:), allocatable :: Coords, PolVec
    real(wpf), dimension(:), allocatable :: PolAng

  end type DynamicVars

  type :: NeighbourType

    integer(wpi), dimension(:), allocatable :: SpecificNeighbours

  end type NeighbourType

  type :: InitValues

    integer(wpi), dimension(:,:), allocatable :: NeighbourMatrix
    real(wpf), dimension(:,:), allocatable :: Coords, PolVec
    real(wpf), dimension(:), allocatable :: PolAng

  end type InitValues

  contains

  subroutine Initialize(ParamAM, DynVarAM, InitSetup)
    implicit none
    type(ParamType), intent(inout) :: ParamAM
    type(DynamicVars), intent(inout) :: DynVarAM
    type(InitValues), intent(in) :: InitSetup

    ! READ NAMELIST

    allocate(DynVarAM%Coords(2,ParamAM%NumPart)) ! Saved as (2,NumPart), as it is usual to need both x and y when considering lengths
    allocate(DynVarAM%)
    


    
  end subroutine Initialize
  
  subroutine InitializeNeighbours(ParamAM,InitSetup,Neighbours)
    !In/Out
    type(ParamType), intent(in) :: ParamAM
    type(InitValues), intent(in) :: InitSetup
    type(NeighbourType), dimension(:), allocatable, intent(inout) :: Neighbours
    

    !Allocate
    integer(wpi) :: iParticle, iNeighbour, NumNeighbours

    allocate(Neighbours(ParamAM%NumPart))

    do iParticle = 1,ParamAM%NumPart
      NumNeighbours = 0
      do iNeighbour = 1,ParamAM%MaxNeighbour
        if ( InitSetup%NeighbourMatrix(iNeighbour,iParticle) /= 0 ) then
          NumNeighbours = NumNeighbours + 1
        end if
      end do
      allocate(Neighbours(iParticle)%SpecificNeighbours(NumNeighbours))
      do iNeighbour = 1, NumNeighbours
        Neighbours(iParticle)%SpecificNeighbours(iNeighbour) = InitSetup%NeighbourMatrix(iNeighbour,iParticle)
      end do
    end do

  end subroutine InitializeNeighbours


end module TypeModule