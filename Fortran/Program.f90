program ActiveSolidProgram
  use PrecMod
  use TypeModule
  use SubRoutineModule
  implicit none

  !Allocate
  type(ParamType) :: ParamAM
  type(DynamicVars) :: DynVarAM
  type(InitValues) :: InitSetup
  type(NeighbourType), dimension(:), allocatable :: Neighbours
  integer(wpi) :: iTimeStep
  character (len = :), allocatable :: ParameterFileName

  ParameterFileName = "Parameters.nml"

  !Initialize coordinates, and polarity direction
  call Initialize(ParamAM, DynVarAM, InitSetup)

  !Initialize neighbour type
  call InitializeNeighbours(ParamAM,InitSetup,Neighbours)

  call DeallocateInit(InitSetup)

  !Starting the timestep loop
  do iTimeStep = 1, ParamAM%NumTimeSteps

    !During this time  in the timestep the old index is old, and the new is being written to
    call TimeStep(ParamAM, DynVarAM, Neighbours)
    !The new index is to be analyzed
    call SystemAnalysis(ParamAM, DynVarAM, Neighbours)
    !The old and new indices are switched. The new becomes old, and the old is to be written over.
    call SwitchIndex(DynVarAM)
  
  end do

  call WriteToFile(DynVarAM)

end program ActiveSolidProgram