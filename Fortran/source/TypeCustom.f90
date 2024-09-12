module TypeModule
  use PrecMod
  use h5fortran
  implicit none
  
  type :: ParamType

    integer(wpi) :: NumPart, MaxNeighbour, NumTimeSteps, SaveEvery

    real(wpf) :: zeta, xi, Fa, k, deltaT, A, R, b, kBoundary

    character(len = :), allocatable :: SaveFileName, IterMethod, BaseDirName, BoundaryMethod

  end type ParamType

  type :: DynamicVars

    real(wpf), dimension(:,:,:), allocatable :: Coords
    real(wpf), dimension(:,:), allocatable :: PolAng
    integer(wpi) :: NewIndex, OldIndex

  end type DynamicVars

  type :: IterationVars

  real(wpf), dimension(:,:,:), allocatable :: kCoords
  real(wpf), dimension(:,:), allocatable :: kPolAng, TmpCoords
  real(wpf), dimension(:), allocatable :: TmpPolAng

  end type IterationVars

  type :: NeighbourType

    integer(wpi), dimension(:), allocatable :: SpecificNeighbours
    real(wpf), dimension(:), allocatable :: EquilLength
    integer(wpi) :: NumNeighbours

  end type NeighbourType

  type :: InitValues

    integer(wpi), dimension(:,:), allocatable :: NeighbourMatrix
    real(wpf), dimension(:,:), allocatable :: Coords
    real(wpf), dimension(:), allocatable :: PolAng
    character(len = :), allocatable :: InitSetupFileName

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

  subroutine Initialize(ParameterFileName,ParamAM, DynVarAM, InitSetup, IterVarAM)
    implicit none
    character(len = *), intent(in) :: ParameterFileName
    type(ParamType), intent(inout) :: ParamAM
    type(DynamicVars), intent(inout) :: DynVarAM
    type(InitValues), intent(inout) :: InitSetup
    type(IterationVars), intent(inout) :: IterVarAM

    call ReadNamelist(ParameterFileName, ParamAM, InitSetup)

    DynVarAM%OldIndex = 1_wpi
    DynVarAM%NewIndex = 2_wpi

    ! Saved as (2,NumPart), as it is usual to need both x and y when considering lengths
    allocate(InitSetup%NeighbourMatrix(ParamAM%MaxNeighbour, ParamAM%NumPart))
    allocate(InitSetup%Coords(2,ParamAM%NumPart))
    allocate(InitSetup%PolAng(ParamAM%NumPart))

    call ReadHDF5(InitSetup, ParamAM)

    ! Saved as (2,NumPart,2) = (x/y,iPart,new/old), as it is usual to need both x and y when considering lengths
    allocate(DynVarAM%Coords(2,ParamAM%NumPart,2)) 
    allocate(DynVarAM%PolAng(ParamAM%NumPart,2))

    ! Only two possibilities, RK4 is the default if RK2 is not chosen.
    if (ParamAM%IterMethod == "RK2") then
      allocate(IterVarAM%kCoords(2,ParamAM%NumPart,2))
      allocate(IterVarAM%kPolAng(ParamAM%NumPart,2))
    else
      allocate(IterVarAM%kCoords(2,ParamAM%NumPart,4))
      allocate(IterVarAM%kPolAng(ParamAM%NumPart,4))
    end if

    allocate(IterVarAM%TmpCoords(2,ParamAM%NumPart))
    allocate(IterVarAM%TmpPolAng(ParamAM%NumPart))
    
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
        - InitSetup%Coords(:,InitSetup%NeighbourMatrix(iNeighbour,iParticle)))
      end do
    end do

  end subroutine InitializeNeighbours

  subroutine DeallocateInit(InitSetup)
    implicit none
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
    integer(wpi) :: NumPart, MaxNeighbour, NumTimeSteps, SaveEvery
    real(wpf) :: zeta, xi, Fa, k, deltaT, b, kBoundary
    character(len=100) :: InitFileName, SaveFileName, IterMethod, BaseDirName, BoundaryMethod
    integer(wpi) :: fu, rc
  
    ! Namelist definition.
    namelist /Parameters/ NumPart, NumTimeSteps, MaxNeighbour, SaveEvery, zeta, xi, Fa, b, k, kBoundary, deltaT, IterMethod, &
      InitFileName, SaveFileName, BaseDirName, BoundaryMethod
    
    
    ! Open and read Namelist file.
    ! We have decided to not do error handling, as all parameters are printed in the end.
    ! This is also a simple implementation with minimal risks.
    open (action='read', file=ParameterFileName, iostat=rc, newunit=fu)
    read (nml=Parameters, iostat=rc, unit=fu)
  
    close (fu)
  
    ! The filenames that are read in have big whitespaces, we trim these.
    InitSetup%InitSetupFileName = trim(InitFileName)

    ParamAM%NumPart = NumPart
    ParamAM%MaxNeighbour = MaxNeighbour
    ParamAM%NumTimeSteps = NumTimeSteps
    ParamAM%zeta = zeta
    ParamAM%xi = xi
    ParamAM%Fa = Fa
    ParamAM%b = b
    ParamAM%k = k
    ParamAM%kBoundary = kBoundary
    ParamAM%deltaT = deltaT
    ParamAM%SaveEvery = SaveEvery
    ParamAM%IterMethod = trim(IterMethod)
    ParamAM%SaveFileName = trim(SaveFileName)
    ParamAM%BaseDirName = trim(BaseDirName)
    ParamAM%BoundaryMethod = trim(BoundaryMethod)

    ! Print parameters
    write(*,"(a)") "----------------------------------------------------------------------"
    write(*,"(a)") "Parameters have been read, the following parameters have been set."
    write(*,"(a,g0)") "Number of particles: - - - - - - - - - - - ", NumPart
    write(*,"(a,g0)") "Maximum number of neighbours:- - - - - - - ", MaxNeighbour
    write(*,"(a,g0)") "Number of timesteps: - - - - - - - - - - - ", NumTimeSteps
    write(*,"(a,g0)") "Timesteps between save:- - - - - - - - - - ", SaveEvery
    write(*,"(a,g0)") "Zeta:- - - - - - - - - - - - - - - - - - - ", zeta
    write(*,"(a,g0)") "Xi:- - - - - - - - - - - - - - - - - - - - ", xi
    write(*,"(a,g0)") "Fa:- - - - - - - - - - - - - - - - - - - - ", Fa
    write(*,"(a,g0)") "b, boundary parameter: - - - - - - - - - - ", b
    write(*,"(a,g0)") "k - Inter particle:- - - - - - - - - - - - ", k
    write(*,"(a,g0)") "k - Boundary:- - - - - - - - - - - - - - - ", kBoundary
    write(*,"(a,g0)") "DeltaT:- - - - - - - - - - - - - - - - - - ", deltaT
    write(*,"(a,g0)") "Iteration method:- - - - - - - - - - - - - ", IterMethod
    write(*,"(a,g0)") "Filename of initial system configuration:- ", InitSetup%InitSetupFileName
    write(*,"(a,g0)") "Filename of savefile:- - - - - - - - - - - ", ParamAM%SaveFileName
    write(*,"(a,g0)") "HDF5 base directory: - - - - - - - - - - - ", ParamAM%BaseDirName
    write(*,"(a,g0)") "Boundary method: - - - - - - - - - - - - - ", ParamAM%BoundaryMethod
    write(*,"(a)") "----------------------------------------------------------------------"
    
  end subroutine ReadNamelist

  subroutine ReadHDF5(InitSetup, ParamAM)
    implicit none
    ! Input
    type(InitValues), intent(inout) :: InitSetup
    type(ParamType), intent(inout) :: ParamAM

    call h5read(InitSetup%InitSetupFileName, '/InitGroup/Coords', InitSetup%Coords)
    call h5read(InitSetup%InitSetupFileName, '/InitGroup/PolAng', InitSetup%PolAng)
    call h5read(InitSetup%InitSetupFileName, '/InitGroup/NeighbourMatrix', InitSetup%NeighbourMatrix)
    call h5read(InitSetup%InitSetupFileName, '/InitGroup/R', ParamAM%R)

  end subroutine ReadHDF5

  subroutine WriteHDF5(DynVarAM, ParamAM, nSave)
    implicit none

    ! Input
    type(DynamicVars), intent(in) :: DynVarAM
    type(ParamType), intent(in) :: ParamAM
    integer(wpi), intent(in) :: nSave

    ! Allocate
    character(len = 64) :: TmpDirNameCoord, TmpDirNamePol, NumStr, DirNameCoord, DirNamePol


    ! Concatenate the base directory and the subdirectory
    TmpDirNameCoord = trim(ParamAM%BaseDirName) // trim("/Coords/")
    TmpDirNamePol = trim(ParamAM%BaseDirName) // trim("/PolAng/")

    ! Convert number to string
    write (NumStr, *) nSave

    ! Adjust length of string
    NumStr = adjustl(NumStr)

    ! Concatenate the number onto the strings
    DirNameCoord = trim(TmpDirNameCoord) // trim(NumStr)
    DirNamePol = trim(TmpDirNamePol) // trim(NumStr)


    ! Writing DynVarAM%Coords(:,:,DynVarAM%OldIndex) to file ParamAM%SaveFileName in HDF5 directory DirName
    call h5write(ParamAM%SaveFileName,trim(DirNameCoord), DynVarAM%Coords(:,:,DynVarAM%OldIndex))

    ! Writing DynVarAM%Coords(:,:,DynVarAM%OldIndex) to file ParamAM%SaveFileName in HDF5 directory DirName
    call h5write(ParamAM%SaveFileName,trim(DirNamePol), DynVarAM%PolAng(:,DynVarAM%OldIndex))

  end subroutine WriteHDF5


end module TypeModule