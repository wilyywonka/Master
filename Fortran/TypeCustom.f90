module TypeModule
  use PrecMod
  implicit none
  
  type :: ParamType

    integer(wpi) :: NumPart, MaxNeighbour, NumTimeSteps

    real(wpf) :: zeta, xi, Fa, k, deltaT

  end type ParamType

  type :: DynamicVars

    real(wpf), dimension(:,:,:), allocatable :: Coords
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
    real(wpf), dimension(:,:), allocatable :: Coords
    real(wpf), dimension(:), allocatable :: PolAng
    character(len = :), allocatable :: InitSetupFileName, NeighbourMatrixFileName

  end type InitValues

  type :: AnalysisParameters

    real(wpf) :: tester

  end type AnalysisParameters
  contains

  function EuclideanNormVec(Coords) result(length)
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
    allocate(InitSetup%PolAng(ParamAM%NumPart))

    ! READ INITFILES

    ! Saved as (2,NumPart,2) = (x/y,iPart,new/old), as it is usual to need both x and y when considering lengths
    allocate(DynVarAM%Coords(2,ParamAM%NumPart,2)) 
    allocate(DynVarAM%PolAng(ParamAM%NumPart,2))
    
    
    DynVarAM%Coords(:,:,DynVarAM%OldIndex) = InitSetup%Coords
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
        Neighbours(iParticle)%EquilLength(iNeighbour) = EuclideanNormVec(InitSetup%Coords(:,iParticle) &
        - InitSetup%Coords(:,iNeighbour))
      end do
    end do

  end subroutine InitializeNeighbours

  subroutine DeallocateInit(InitSetup)
    !Input
    type(InitValues), intent(inout) :: InitSetup

    deallocate(InitSetup%NeighbourMatrix)
    deallocate(InitSetup%Coords)
    deallocate(InitSetup%PolAng)

  end subroutine DeallocateInit

  subroutine ReadNamelist(ParameterFileName, ParamAM, InitSetup)
    ! This is a modified version of the code available at https://cyber.dabamos.de/programming/modernfortran/namelists.html
    implicit none
    character(len = *), intent(in) :: ParameterFileName
    type(ParamType), intent(inout) :: ParamAM
    type(InitValues), intent(inout) :: InitSetup

  
    ! Allocate
    integer(wpi) :: NumPart, MaxNeighbour, NumTimeSteps
    real(wpf) :: zeta, xi, Fa, k, deltaT
    character(len=100) :: InitFileName, InitNeighbour
    integer(wpi) :: fu, rc
  
    ! Namelist definition.
    namelist /Parameters/ NumPart, MaxNeighbour, NumTimeSteps, zeta, xi, Fa, k, deltaT, InitFileName, InitNeighbour
  
    ! Open and read Namelist file.
    ! We have decided to not do error handling, as all parameters are printed in the end.
    ! This is also a simple implementation with minimal risks.
    open (action='read', file=ParameterFileName, iostat=rc, newunit=fu)
    read (nml=Parameters, iostat=rc, unit=fu)
  
    close (fu)
  
    ! The filenames that are read in have big whitespaces, we trim these.
    InitSetup%InitSetupFileName = trim(InitFileName)
    InitSetup%NeighbourMatrixFileName = trim(InitNeighbour)

    ParamAM%NumPart = NumPart
    ParamAM%MaxNeighbour = MaxNeighbour
    ParamAM%NumTimeSteps = NumTimeSteps
    ParamAM%zeta = zeta
    ParamAM%xi = xi
    ParamAM%Fa = Fa
    ParamAM%k = k
    ParamAM%deltaT = deltaT

    ! Print parameters
    write(*,"(a)") "----------------------------------------------------------------------"
    write(*,"(a)") "Parameters have been read, the following parameters have been set."
    write(*,"(a,g0)") "Number of particles: - - - - - - - - - - - ", NumPart
    write(*,"(a,g0)") "Maximum number of neighbours:- - - - - - - ", MaxNeighbour
    write(*,"(a,g0)") "Number of timesteps: - - - - - - - - - - - ", NumTimeSteps
    write(*,"(a,g0)") "Zeta:- - - - - - - - - - - - - - - - - - - ", zeta
    write(*,"(a,g0)") "Xi:- - - - - - - - - - - - - - - - - - - - ", xi
    write(*,"(a,g0)") "Fa:- - - - - - - - - - - - - - - - - - - - ", Fa
    write(*,"(a,g0)") "k: - - - - - - - - - - - - - - - - - - - - ", k
    write(*,"(a,g0)") "DeltaT:- - - - - - - - - - - - - - - - - - ", deltaT
    write(*,"(a,g0)") "Filename of initial system configuration:- ", InitSetup%InitSetupFileName
    write(*,"(a,g0)") "Filename of neighbour matrix:- - - - - - - ", InitSetup%NeighbourMatrixFileName
    write(*,"(a)") "----------------------------------------------------------------------"
    
  end subroutine ReadNamelist


end module TypeModule