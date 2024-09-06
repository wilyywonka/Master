program ActiveSolidProgram
  use SubRoutineModule
  implicit none

  ! Allocate
  type(ParamType) :: ParamAM
  type(DynamicVars) :: DynVarAM
  type(InitValues) :: InitSetup
  type(NeighbourType), dimension(:), allocatable :: Neighbours
  type(IterationVars) :: IterVarAM
  character (len = :), allocatable :: ParameterFileName
  integer(wpi) :: iTimeStep

  ParameterFileName = "Parameters.nml"

  ! Initialize coordinates, and polarity direction
  call Initialize(ParameterFileName, ParamAM, DynVarAM, InitSetup, IterVarAM)

  ! Initialize neighbour type
  call InitializeNeighbours(ParamAM,InitSetup,Neighbours)

  ! The initialization is now done, and the InitSetup can be deallocated
  call DeallocateInit(InitSetup)

  !Starting the timestep loop
  do iTimeStep = 1, ParamAM%NumTimeSteps

    ! During this time  in the timestep the old index is old, and the new is being written to
    call TimeStep(ParamAM, DynVarAM, IterVarAM, Neighbours)
 
    ! The new index is to be analyzed
!    call SystemAnalysis(ParamAM, DynVarAM, Neighbours)
 
    ! The old and new indices are switched. The new becomes old, and the old is to be written over.
    call SwitchIndex(DynVarAM)
  
    ! Save the system if the index is 
    if (modulo(iTimeStep, ParamAM%SaveEvery) == 0) then
      call WriteToFile(DynVarAM)
    end if

  end do
  

end program ActiveSolidProgram