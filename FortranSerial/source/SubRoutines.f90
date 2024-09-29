module SubRoutineModule
  use TypeModule
  implicit none
  
  contains

  subroutine TimeStep(ParamAM, DynVarAM, IterVarAM, Neighbours)
    implicit none
    ! Input
    type(ParamType), intent(in) :: ParamAM
    type(DynamicVars), intent(inout) :: DynVarAM
    type(IterationVars), intent(inout) :: IterVarAM
    type(NeighbourType), dimension(:), intent(in) :: Neighbours
    
    ! Allocate
    integer(wpi) :: iParticle


    if (ParamAM%IterMethod == "RK2") then

      ! ----------------------------
      ! Runge-Kutta 2
      ! ----------------------------

      ! Loop over particles, this loop has been deliberately made to be very easy to parallelize
      do iParticle = 1,ParamAM%NumPart
        ! Calling the RK-Step
        call RKStep(ParamAM, IterVarAM, Neighbours, iParticle, DynVarAM%Coords(:,:,DynVarAM%OldIndex), &
          DynVarAM%PolAng(:,DynVarAM%OldIndex), 1_wpi)

        IterVarAM%TmpCoords(:,iParticle) = DynVarAM%Coords(:,iParticle,DynVarAM%OldIndex) + &
          ParamAM%deltaT*IterVarAM%kCoords(:,iParticle,1)
  
        IterVarAM%TmpPolAng(iParticle) = DynVarAM%PolAng(iParticle,DynVarAM%OldIndex) + &
          ParamAM%deltaT*IterVarAM%kPolAng(iParticle,1)
      end do

      ! Loop over particles, this loop has been deliberately made to be very easy to parallelize
      do iParticle = 1,ParamAM%NumPart
        ! Calling the RK-Step
        call RKStep(ParamAM, IterVarAM, Neighbours, iParticle, IterVarAM%TmpCoords, IterVarAM%TmpPolAng, 2_wpi)

        ! All RK-Steps are done and now the new iteration is written
        DynVarAM%Coords(:,iParticle,DynVarAM%NewIndex) = DynVarAM%Coords(:,iParticle,DynVarAM%OldIndex) + &
          (ParamAM%deltaT/2_wpi)*(IterVarAM%kCoords(:,iParticle,1_wpi) + IterVarAM%kCoords(:,iParticle,2_wpi))
  
        DynVarAM%PolAng(iParticle,DynVarAM%NewIndex) = DynVarAM%PolAng(iParticle,DynVarAM%OldIndex) + &
        (ParamAM%deltaT/2_wpi)*(IterVarAM%kPolAng(iParticle,1_wpi) + IterVarAM%kPolAng(iParticle,2_wpi))
      end do
    

    else

      ! ----------------------------
      ! Runge-Kutta 4
      ! ----------------------------

      ! Loop over particles, this loop has been deliberately made to be very easy to parallelize
      do iParticle = 1,ParamAM%NumPart
        ! Calling the RK-Step
        call RKStep(ParamAM, IterVarAM, Neighbours, iParticle, DynVarAM%Coords(:,:,DynVarAM%OldIndex), &
        DynVarAM%PolAng(:,DynVarAM%OldIndex), 1_wpi)

        IterVarAM%TmpCoords(:,iParticle) = DynVarAM%Coords(:,iParticle,DynVarAM%OldIndex) + &
          (ParamAM%deltaT/2)*IterVarAM%kCoords(:,iParticle,1)

        IterVarAM%TmpPolAng(iParticle) = DynVarAM%PolAng(iParticle,DynVarAM%OldIndex) + &
          (ParamAM%deltaT/2)*IterVarAM%kPolAng(iParticle,1)
      end do

      ! Loop over particles, this loop has been deliberately made to be very easy to parallelize
      do iParticle = 1,ParamAM%NumPart
        ! Calling the RK-Step
        call RKStep(ParamAM, IterVarAM, Neighbours, iParticle, IterVarAM%TmpCoords, IterVarAM%TmpPolAng, 2_wpi)

        IterVarAM%TmpCoords(:,iParticle) = DynVarAM%Coords(:,iParticle,DynVarAM%OldIndex) + &
          (ParamAM%deltaT/2)*IterVarAM%kCoords(:,iParticle,2)
        IterVarAM%TmpPolAng(iParticle) = DynVarAM%PolAng(iParticle,DynVarAM%OldIndex) + &
          (ParamAM%deltaT/2)*IterVarAM%kPolAng(iParticle,2)
      end do

      ! Loop over particles, this loop has been deliberately made to be very easy to parallelize
      do iParticle = 1,ParamAM%NumPart
        ! Calling the RK-Step
        call RKStep(ParamAM, IterVarAM, Neighbours, iParticle, IterVarAM%TmpCoords, IterVarAM%TmpPolAng, 3_wpi)

        IterVarAM%TmpCoords(:,iParticle) = DynVarAM%Coords(:,iParticle,DynVarAM%OldIndex) + &
          ParamAM%deltaT*IterVarAM%kCoords(:,iParticle,3)

        IterVarAM%TmpPolAng(iParticle) = DynVarAM%PolAng(iParticle,DynVarAM%OldIndex) + &
          ParamAM%deltaT*IterVarAM%kPolAng(iParticle,3)
      end do
      
      ! Loop over particles, this loop has been deliberately made to be very easy to parallelize
      do iParticle = 1,ParamAM%NumPart
        ! Calling the RK-Step
        call RKStep(ParamAM, IterVarAM, Neighbours, iParticle, IterVarAM%TmpCoords, IterVarAM%TmpPolAng, 4_wpi)

        ! All RK-Steps are done and now the new iteration is written
        DynVarAM%Coords(:,iParticle,DynVarAM%NewIndex) = DynVarAM%Coords(:,iParticle,DynVarAM%OldIndex) + &
          (ParamAM%deltaT/6_wpi)*(IterVarAM%kCoords(:,iParticle,1_wpi) + 2_wpi*IterVarAM%kCoords(:,iParticle,2_wpi) + &
            2_wpi*IterVarAM%kCoords(:,iParticle,3_wpi) + IterVarAM%kCoords(:,iParticle,4_wpi))

        DynVarAM%PolAng(iParticle,DynVarAM%NewIndex) = DynVarAM%PolAng(iParticle,DynVarAM%OldIndex) + &
          (ParamAM%deltaT/6_wpi)*(IterVarAM%kPolAng(iParticle,1_wpi) + 2_wpi*IterVarAM%kPolAng(iParticle,2_wpi) + &
            2_wpi*IterVarAM%kPolAng(iParticle,3_wpi) + IterVarAM%kPolAng(iParticle,4_wpi))

      end do
    end if

  end subroutine TimeStep


  subroutine RKStep(ParamAM, IterVarAM, Neighbours, iParticle, TmpCoords, TmpPolAng, Iteration)
    implicit none
    !Input
    type(ParamType), intent(in) :: ParamAM
    type(IterationVars), intent(inout) :: IterVarAM
    type(NeighbourType), dimension(:), intent(in) :: Neighbours
    integer(wpi), intent(in) :: iParticle
    real(wpf), dimension(:,:), intent(in) :: TmpCoords
    real(wpf), dimension(:), intent(in) :: TmpPolAng
    integer(wpi), intent(in) :: Iteration
    !Allocate
    integer(wpi) :: iNeighbour
    real(wpf), dimension(2) :: TotalForce, PolVec
    real(wpf) :: TotalTorque, PartCentreRad

    TotalForce = (/0._wpf, 0._wpf/)
    PolVec = (/cos(TmpPolAng(iParticle)), sin(TmpPolAng(iParticle))/)

    ! Elasticity method
    select case (ParamAM%ElasticityMethod)
      case ("LinearElasticity")
        ! Looping over all the neighbours of this particle to add up the linear elastic forces
        do iNeighbour = 1, Neighbours(iParticle)%NumNeighbours
          TotalForce = TotalForce + LinearElasticForce(ParamAM,TmpCoords,iParticle,Neighbours(iParticle),iNeighbour)
        end do
      case ("FENEElasticity")
        ! Looping over all the neighbours of this particle to add up the nonlinear elastic forces
        do iNeighbour = 1, Neighbours(iParticle)%NumNeighbours
          TotalForce = TotalForce + FENEElasticForce(ParamAM,TmpCoords,iParticle,Neighbours(iParticle),iNeighbour)
        end do
      case ("NonLinearElasticity")
        ! Looping over all the neighbours of this particle to add up the nonlinear elastic forces
        do iNeighbour = 1, Neighbours(iParticle)%NumNeighbours
          TotalForce = TotalForce + NonLinearElasticForce(ParamAM,TmpCoords,iParticle,Neighbours(iParticle),iNeighbour)
        end do
      case default
        print*, "Invalid elasticity method!"
        stop
    end select

    PartCentreRad = EuclideanNormVec(TmpCoords(:,iParticle))

    ! Boundary condition
    select case (ParamAM%BoundaryMethod)
      case ("AttractiveRepellingBoundary")
        if (ParamAM%R-PartCentreRad < 2*ParamAM%b) then
          ! Adding the attractive-repelling boundary force
          TotalForce = TotalForce + BoundaryAttrRep(ParamAM,TmpCoords(:,iParticle), PartCentreRad)
        end if
      case ("RepellingBoundary")
        if (PartCentreRad>ParamAM%R) then
          ! Adding the repelling boundary force
          TotalForce = TotalForce + BoundaryRep(ParamAM,TmpCoords(:,iParticle), PartCentreRad)
        end if 
      case default
        print*, "Invalid boundary method!"
        stop
    end select

    ! Calculating the torque using the TotalForce, as it per now only contains elastic force, and is thus the total elastic force
    TotalTorque = PolarizationTorque(ParamAM,TotalForce,PolVec)

    ! Adding on the polarization force, to the elastic forces to produce the total force and dividing by zeta to get the "velocity"
    TotalForce = (TotalForce + PolarizationForce(ParamAM,PolVec))/ParamAM%zeta

    ! Setting the k-vectors using the "velocity" values
    IterVarAM%kCoords(:,iParticle, Iteration) = TotalForce
    IterVarAM%kPolAng(iParticle, Iteration) = TotalTorque

  end subroutine RKStep

  function BoundaryAttrRep(ParamAM, TmpPartCoords, PartCentreRad) result(BoundForce)
    ! Input
    type(ParamType), intent(in) :: ParamAM
    real(wpf), dimension(2), intent(in) :: TmpPartCoords
    real(wpf), intent(in) :: PartCentreRad
    ! Allocate
    real(wpf) :: ProxLength
    ! Output
    real(wpf), dimension(2) :: BoundForce
  
    ProxLength = ParamAM%R-PartCentreRad

    BoundForce = -((ProxLength**2 - 2*ProxLength*ParamAM%b)/PartCentreRad)*TmpPartCoords*ParamAM%kBoundary

  end function BoundaryAttrRep

  function BoundaryRep(ParamAM, TmpPartCoords, PartCentreRad) result(BoundForce)
     ! Input
    type(ParamType), intent(in) :: ParamAM
    real(wpf), dimension(2), intent(in) :: TmpPartCoords
    real(wpf), intent(in) :: PartCentreRad
    ! Allocate
    real(wpf) :: ProxLength
    ! Output
    real(wpf), dimension(2) :: BoundForce
  
    ProxLength = PartCentreRad-ParamAM%R

    BoundForce = -(ProxLength/PartCentreRad)*TmpPartCoords*ParamAM%kBoundary
   
  end function BoundaryRep


  function LinearElasticForce(ParamAM, TmpCoords, iPart, Neighbour, iNeigh) result(ElForce)
    implicit none
    ! Input
    type(ParamType), intent(in) :: ParamAM
    real(wpf), dimension(:,:), intent(in) :: TmpCoords
    integer(wpi), intent(in) :: iPart, iNeigh
    type(NeighbourType), intent(in) :: Neighbour
    ! Allocate
    real(wpf), dimension(2) :: DeltaCoords
    real(wpf) :: length
    ! Output
    real(wpf), dimension(2) :: ElForce

    !Calculating the difference in coordinates between the particle and the neighbour
    DeltaCoords = TmpCoords(:,iPart) - TmpCoords(:,Neighbour%SpecificNeighbours(iNeigh))

    !Calculating the length
    length = EuclideanNormVec(DeltaCoords)

    !Force is calculated as -k(l-l_0)\vec(d), where d is in the direction of the displacement, as a unitary vector
    Elforce = -ParamAM%k*(length-Neighbour%EquilLength(iNeigh))*DeltaCoords/length

  end function


  function NonLinearElasticForce(ParamAM, TmpCoords, iPart, Neighbour, iNeigh) result(ElForce)
    implicit none
    ! Input
    type(ParamType), intent(in) :: ParamAM
    real(wpf), dimension(:,:), intent(in) :: TmpCoords
    integer(wpi), intent(in) :: iPart, iNeigh
    type(NeighbourType), intent(in) :: Neighbour
    ! Allocate
    real(wpf), dimension(2) :: DeltaCoords
    real(wpf) :: length
    ! Output
    real(wpf), dimension(2) :: ElForce

    !Calculating the difference in coordinates between the particle and the neighbour
    DeltaCoords = TmpCoords(:,iPart) - TmpCoords(:,Neighbour%SpecificNeighbours(iNeigh))

    !Calculating the length
    length = EuclideanNormVec(DeltaCoords)

    !Force is calculated as -k(l-l_0)\vec(d), where d is in the direction of the displacement, as a unitary vector
    Elforce = -(ParamAM%k*(length-Neighbour%EquilLength(iNeigh)) + ParamAM%kNonLin*(length-Neighbour%EquilLength(iNeigh))**3) &
      *DeltaCoords/length

  end function


  function FENEElasticForce(ParamAM, TmpCoords, iPart, Neighbour, iNeigh) result(ElForce)
    implicit none
    ! Input
    type(ParamType), intent(in) :: ParamAM
    real(wpf), dimension(:,:), intent(in) :: TmpCoords
    integer(wpi), intent(in) :: iPart, iNeigh
    type(NeighbourType), intent(in) :: Neighbour
    ! Allocate
    real(wpf), dimension(2) :: DeltaCoords
    real(wpf) :: length
    ! Output
    real(wpf), dimension(2) :: ElForce

    !Calculating the difference in coordinates between the particle and the neighbour
    DeltaCoords = TmpCoords(:,iPart) - TmpCoords(:,Neighbour%SpecificNeighbours(iNeigh))

    !Calculating the length
    length = EuclideanNormVec(DeltaCoords)

    !Force is calculated as -k(l-l_0)\vec(d), where d is in the direction of the displacement, as a unitary vector
    Elforce = -(ParamAM%k*(length-Neighbour%EquilLength(iNeigh))/(1-(length/Neighbour%EquilLength(iNeigh)-1)**2)) &
      *DeltaCoords/length

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

  subroutine TruncateAngles(ParamAM,DynVarAM)
    implicit none
    ! Input/Output
    type(ParamType), intent(in) :: ParamAM
    type(DynamicVars), intent(inout) :: DynVarAM

    DynVarAM%PolAng = modulo(DynVarAM%PolAng, 2*ParamAM%pi)

  end subroutine TruncateAngles

  subroutine CalculateDisplacement(ParamAM,DynVarAM)
    implicit none
    type(ParamType), intent(in) :: ParamAM
    type(DynamicVars), intent(inout) :: DynVarAM
    integer(wpi) :: iParticle
    real(wpf), dimension(2) :: Displacement
    do iParticle = 1,ParamAM%NumPart
      Displacement = (DynVarAM%Coords(:,iParticle,DynVarAM%OldIndex) - DynVarAM%Coords(:,iParticle,DynVarAM%NewIndex))
      DynVarAM%DisplacementVectors(:,iParticle) = Displacement/EuclideanNormVec(Displacement)
    end do

  end subroutine CalculateDisplacement
  
end module SubRoutineModule
