using Plots
using Delaunay
using HDF5


# ---------------------------------------------------------------
# SET PARAMETERS - SIMULATION -
NPart::Int64 = 50000
NSimulationIterations::Int64 = 5000
SaveEverySimulation::Int64 = 10
kSimulation::Float64 = 60.0
kNonLin::Float64 = 60.0
kBoundarySimulation::Float64 = 20.0
xiSimulation::Float64 = 5.0
zetaSimulation::Float64 = 1.0
FaSimulation::Float64 = 1.0
b::Float64 = 0.3
deltaTSimulation::Float64 = 0.001
IterMethodSimulation::String = "RK4"
BaseDirName::String = "/SaveParams"
FileName::String = "InitState.h5"
InitFileName::String = "../../../HDF5Files/"*FileName
SaveFileName::String = "../../../HDF5Files/SaveFiles.h5"
BoundaryMethodSimulation::String = "RepellingBoundary" #RepellingBoundary : AttractiveRepellingBoundary
ElasticityMethod::String = "LinearElasticity" #LinearElasticity : NonLinearElasticity : FENEElasticity


# SET PARAMETERS - SAVE -
PathName::String = "./HDF5Files/"
NameListPath::String = "./Parameters/Parameters.nml"


# SET PARAMETERS - INITSTATE -
ρ::Float64 = 1.0
spread::Vector{Float64} = [0.85,1.15]
finalMeanRad::Float64 = 1.0
growSize::Float64 = 0.3
D::Float64 = 0.01
NTimeSteps::Int64 = 10000
Extrasteps::Int64 = 20000 # May need to be even longer
A::Float64 = 5.0
dt::Float64 = 0.02
saveEvery::Int64 = 1
BrownianMethod::String = "BrownianTemperatureCartesian" # BrownianTemperaturePolar : BrownianTemperatureCartesian
BoundaryMethod::String = "AttractiveRepellingBoundary" #RepellingBoundary : AttractiveRepellingBoundary
ExtraSlots::Int64 = 5
MoveMultiplier::Float64 = 1.2

# SET PARAMETERS - PLOTTING -
PlotEnd::Bool = false
CirclePlot::Bool = false
DelaunayPlot::Bool = false
NeighbourPlot::Bool = false

# ---------------------------------------------------------------

function InitCoords!(Coords::Array{Float64}, NPart::Int64, R::Float64, iOld::Int64, randAngle::Vector{Float64}, randRadius::Vector{Float64})
    # Define N random angles
    randAngle[:] = rand(NPart)*2*pi
    # Define N random, uniformally distributed values. Multiply these with R², and take the square root.
    # This will ensure a uniformally distributed points, even when using polar coordinates.
    randRadius[:] = sqrt.(rand(NPart)*(R^2))    
    # Set the coordinates using the cosine and sine of the angles, times the radius
    Coords[1,:,iOld] = cos.(randAngle).*randRadius
    Coords[2,:,iOld] = sin.(randAngle).*randRadius
    # Return the coordinates
    return nothing
end

#MAY NEED CHANGES; NO DIVISION/MULTIPLICATION WITH E.G. √3
function BrownianTemperaturePolar!(RandDisp::Matrix{Float64}, randAngle::Vector{Float64}, randRadius::Vector{Float64}, NPart::Int64, D::Float64, dt::Float64)
    # Define N random angles
    randAngle[:] = rand(NPart)*2*pi
    # Define N random, uniformally distributed values and take the square root.
    # This will ensure a uniformally distributed points, even when using polar coordinates.
    randRadius[:] = sqrt.(rand(NPart))
    # Set the coordinates using the cosine and sine of the angles, times the radius
    RandDisp[1,:] = cos.(randAngle).*randRadius.*sqrt(D*dt)
    RandDisp[2,:] = sin.(randAngle).*randRadius.*sqrt(D*dt)
    return nothing
end

function BrownianTemperatureCartesian!(RandDisp::Matrix{Float64}, randAngle::Vector{Float64}, randRadius::Vector{Float64}, NPart::Int64, D::Float64, dt::Float64)
    # Set the coordinates using the cosine and sine of the angles, times the radius
    @view(RandDisp[1,:]) .= (2 .*rand(NPart).-1).*sqrt(D*dt*6)
    @view(RandDisp[2,:]) .= (2 .*rand(NPart).-1).*sqrt(D*dt*6)
    #RandDisp[2,:] = (2*rand(NPart)-1).*sqrt(D*dt*3)
    return nothing
end

function RepellingBoundary!(Coords::Matrix{Float64}, BoundaryForce::Matrix{Float64}, R::Float64, radCenter::Vector{Float64}, b::Float64, iParticle::Int64)
    # Initialize a force vector CoordArray, BoundaryForce, R, radCenter, b, iParticle
    @view(BoundaryForce[:,iParticle]) .= 0.0
    # Check if the CENTER of the particles is outside the circle
    if radCenter[iParticle] > R
        # Add it to the force vector
        @view(BoundaryForce[:,iParticle]) .+= -(((radCenter[iParticle]-R))/radCenter[iParticle]).*[Coords[1,iParticle],Coords[2,iParticle]]
    end
    return nothing
end

function AttractiveRepellingBoundary!(Coords::Matrix{Float64}, BoundaryForce::Matrix{Float64}, R::Float64, radCenter::Vector{Float64}, b::Float64, iParticle::Int64)
    # Initialize a force vector
    BoundaryForce[:,iParticle] .= 0.0
    # Check if the CENTER of the particles is within 2b, the attractive/repulsive zone
    if (R-radCenter[iParticle]) < 2*b
        # Add it to the force vector
        BoundaryForce[:,iParticle] += -(((R-radCenter[iParticle])^2 -2*(R-radCenter[iParticle])*b)/radCenter[iParticle]).*[Coords[1,iParticle],Coords[2,iParticle]]
    end
    return nothing
end

# function InitCloseParticles!(CloseParticlesArray::Vector{Vector{Int64}}, NPart::Int64)
#     for i in 1:NPart
#         push!(CloseParticlesArray,zeros(Int64,1))
#     end
#     return nothing
# end

# function CloseParticles!(CloseParticlesArray::Vector{Vector{Int64}}, IndexVec::Vector{Int64}, CoordArray::Array{Float64}, iOld::Int64, radArray::Vector{Float64}, iParticle::Int64)
#     CloseParticlesArray[iParticle] = IndexVec[1:end .!= iParticle][findall((CoordArray[1,1:end .!= iParticle,iOld].-CoordArray[1,iParticle,iOld]).^2 .+ (CoordArray[2,1:end .!= iParticle,iOld].-CoordArray[2,iParticle,iOld]).^2 .< (radArray[1:end .!= iParticle].+radArray[iParticle]).^2)]
#     return nothing
# end

function InitNeighbours(NeighbourCount::Matrix{Int64}, NCells::Int64, CoordArrayShift::Matrix{Float64},l::Float64,DivVec::Vector{Float64},ExtraSlots::Int64)
    iCell::Int64 = 0
    jCell::Int64 = 0
    IndexPlacement::Int64 = 0
    MaxPart::Int64 = 0
    NeighbourCount .= 0
    for Vec in eachslice(CoordArrayShift,dims=2)

        DivVec .= Vec./l
      
        iCell = ceil(Int64,DivVec[1])

        jCell = ceil(Int64,DivVec[2])

        @inbounds NeighbourCount[iCell,jCell] += 1
    end

    MaxPart = maximum(NeighbourCount) + ExtraSlots
    NeighbourIndices::Array{Int64} = zeros(Int64,MaxPart,NCells,NCells)
    NeighbourCount .= 0

    for (iParticle,Vec) in enumerate(eachslice(CoordArrayShift,dims=2))

        DivVec .= Vec./l

        iCell = ceil(Int64,DivVec[1])

        jCell = ceil(Int64,DivVec[2])

        @inbounds NeighbourCount[iCell,jCell] += 1
        IndexPlacement = NeighbourCount[iCell,jCell]
        NeighbourIndices[IndexPlacement,iCell,jCell] = iParticle

    end
    return NeighbourIndices, MaxPart
end

function UpdateNeighbours(NeighbourIndices::Array{Int64},NeighbourCount::Matrix{Int64}, NCells::Int64, CoordArrayShift::Matrix{Float64},l::Float64,DivVec::Vector{Float64}, MaxPart::Int64, ExtraSlots::Int64)
    iCell::Int64 = 0
    jCell::Int64 = 0
    IndexPlacement::Int64 = 0
    NeighbourCount .= 0
    for (iParticle,Vec) in enumerate(eachslice(CoordArrayShift,dims=2))

        DivVec .= Vec./l

        iCell = ceil(Int64,DivVec[1])

        jCell = ceil(Int64,DivVec[2])

        @inbounds @view(NeighbourCount[iCell,jCell]) .+= 1::Int64

        @inbounds IndexPlacement = @view(NeighbourCount[iCell,jCell])[1]
        if IndexPlacement <= MaxPart
            @inbounds NeighbourIndices[IndexPlacement,iCell,jCell] = iParticle
        else
            println("The number of particles in the cell is greater than MaxNeighbour! Reshaping NeighbourMatrix and recursively calling function.")
            MaxPart += 5
            NeighbourIndices = zeros(Int64,MaxPart+ExtraSlots,NCells,NCells)
            MaxPart, NeighbourIndices = UpdateNeighbours(NeighbourIndices,NeighbourCount, NCells, CoordArrayShift, l, DivVec, MaxPart, ExtraSlots)
            return MaxPart, NeighbourIndices
        end

    end
    return MaxPart, NeighbourIndices
end

function CalcElasticRepulsion!(NeighbourCount::Matrix{Int64}, NeighbourIndices::Array{Int64}, forceVec::Matrix{Float64}, radArray::Vector{Float64}, CountCoords::CartesianIndices{2, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}}, DeltaVec::Matrix{Float64}, TempVec::Matrix{Float64}, DeltaDist::Vector{Float64}, CoordArray::Matrix{Float64}, NCells::Int64)
    forceVec .= 0.0
    #CoordPairs = eachslice(CoordArray[:,:], dims= 2)
    #Loop over
    Threads.@threads for IdxCart in CountCoords
        NumInCell::Int64 = @view(NeighbourCount[IdxCart])[1]
        # nested task error: MethodError: no method matching firstindex(::Base.Pairs{CartesianIndex{2}, Int64, CartesianIndices{2, Tuple{…}}, Matrix{Int64}})
        if NumInCell > 0
            # Defining vectors and local variables
            IdxCurrent::Vector{Int64} = [IdxCart[1],IdxCart[2]]
            IdxCells::Vector{Int64} = [IdxCart[1],IdxCart[2]]
            PartInCell::Vector{Int64} = NeighbourIndices[1:NumInCell,IdxCells[1],IdxCells[2]]
            NumInNeighCell::Int64 = 0

            # Defining all of the relative cells to be looped over
            iCells::Vector{Int64} = [-1, 0 ,1]
            jCells::Vector{Int64} = [-1, 0 ,1]
            iNeighbour::Int64 = 0
           # TestIndices::Vector{Int64} = []

            # In the edge cases, the extra indices will not be included
            if IdxCells[1] == 1 
                iCells = [ 0, 1]
            elseif IdxCells[1] == NCells
                iCells = [-1, 0]
            end
            if IdxCells[2] == 1 
                jCells = [ 0, 1]
            elseif IdxCells[2] == NCells
                jCells = [-1, 0]
            end

            # Loop over neighbouring cells
            for iCellIdx in iCells, jCellIdx in jCells
                IdxCells .= IdxCurrent
                @view(IdxCells[1]) .+= iCellIdx
                @view(IdxCells[2]) .+= jCellIdx
                # If the current cell to check is the cell we are in
                if IdxCells == IdxCurrent && NumInCell > 1
                    # Loop over particles and neighbours
                    for iParticle in PartInCell, iNeighbour in PartInCell
                        # Check if the indexes are the same
                        if iNeighbour != iParticle
                            # Save distance between particles
                            @inbounds @view(DeltaVec[:,iParticle]) .= @view(CoordArray[:,iParticle]) .- @view(CoordArray[:,iNeighbour]) #CoordPairs[iParticle] - CoordPairs[iNeighbour]
                            # Calculate the force between the particles
                            LinearElasticity!(DeltaVec, TempVec, DeltaDist, forceVec, radArray, iNeighbour, iParticle)
                        end
                    end
                # Cell is not the same as the current cell
                elseif IdxCells != IdxCurrent
                    # Number of particles in the neighbouring cell
                    NumInNeighCell = NeighbourCount[IdxCells[1],IdxCells[2]]
                    # Check if there are particles in the neighbouring cell
                    if NumInNeighCell > 0
                       # TestIndices = NeighbourIndices[1:NumInNeighCell,IdxCells[1],IdxCells[2]]
                        # Loop over particles and neighbours
                        for iParticle in PartInCell, iNeighbourIdx in 1:NumInNeighCell
                            # Save distance between particles
                            iNeighbour = @view(NeighbourIndices[iNeighbourIdx,IdxCells[1],IdxCells[2]])[1]
                            @inbounds @view(DeltaVec[:,iParticle]) .= @view(CoordArray[:,iParticle]) .- @view(CoordArray[:,iNeighbour]) #CoordPairs[iParticle] .- CoordPairs[iNeighbour]
                            # Calculate the force between the particles
                            LinearElasticity!(DeltaVec, TempVec, DeltaDist, forceVec, radArray, iNeighbour, iParticle)
                        end
                    end
                end

            end
        end
    end
    return nothing
end


function LinearElasticity!(DeltaVec::Matrix{Float64}, TempVec::Matrix{Float64}, DeltaDist::Vector{Float64}, forceVec::Matrix{Float64}, radArray::Vector{Float64}, iNeighbour::Int64, iParticle::Int64)
    # Setting Tempvec as squared of Deltavec
    @inbounds @view(TempVec[:,iParticle]) .= @view(DeltaVec[:,iParticle]).^2
    # Saving the distance between the particles
    @inbounds @view(DeltaDist[iParticle]) .= sqrt(TempVec[1,iParticle]+TempVec[2,iParticle])
    # Checking if the particles are overlapping
    if @view(DeltaDist[iParticle])[1] < (radArray[iNeighbour] + radArray[iParticle])
        # Calculating the force of the overlap and adding the force to the force vector
        @inbounds @view(forceVec[:,iParticle]) .+= (abs.(@view(radArray[iParticle]) .+ @view(radArray[iNeighbour]) .- @view(DeltaDist[iParticle])).^(3/2)) .* @view(DeltaVec[:,iParticle])./@view(DeltaDist[iParticle])
    end
    return nothing
end


function CalcBoundaryRepulsion!(NPart::Int64,CoordArray::Matrix{Float64},BoundaryMethod!::Function,forceVec::Matrix{Float64},CoordArraySliced::Vector{Tuple{Int64, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}},TempVec::Matrix{Float64},R::Float64,radCenter::Vector{Float64},b::Float64,BoundaryForce::Matrix{Float64})
    #CoordArraySliced .= collect(enumerate(eachslice(CoordArray,dims=2)))
    Threads.@threads for iParticle in 1:NPart
        # Assign Vec squared to TempVec
        @inbounds @view(TempVec[:,iParticle]) .= @view(CoordArray[:,iParticle]).^2
        # Calculate the distance from the center of the circle
        @inbounds @view(radCenter[iParticle]) .= sqrt(TempVec[1,iParticle].+TempVec[2,iParticle])
        # Calculate the boundary force and add to the forceVec array
        BoundaryMethod!(CoordArray, BoundaryForce, R, radCenter, b, iParticle)
        # Add the boundary force to the force vector
        @inbounds @view(forceVec[:,iParticle]) .+= @view(BoundaryForce[:,iParticle])
    end
end

function CalcNextStep!(CoordArray::Array{Float64},iNew::Int64,iOld::Int64,forceVec::Matrix{Float64},A::Float64,dt::Float64,RandDisp::Matrix{Float64},NPart::Int64)
    Threads.@threads for iParticle in 1:NPart
        @view(CoordArray[:,iParticle,iNew]) .= @view(CoordArray[:,iParticle,iOld]) .+ @view(forceVec[:,iParticle]).*A.*dt .+ @view(RandDisp[:,iParticle])
    end
end

# NPart = 1000000
# iOld = 1
# finalMeanRad = 1
# ρ = 1
# CoordArray::Array{Float64} = zeros(2,NPart,2)
# R::Float64 = sqrt((NPart*(finalMeanRad)^2)/(ρ))
# randAngle::Vector{Float64} = zeros(NPart)
# randRadius::Vector{Float64} = zeros(NPart)

# InitCoords!(CoordArray,NPart,R,iOld,randAngle,randRadius)

# using Plots
# using BenchmarkTools

# #scatter(CoordArray[1,:,1],CoordArray[2,:,1],size=(900,900))

# MoveMultiplier = 1.2

# CoordArrayShift::Matrix{Float64} = zeros(2,NPart)

# CoordArrayShift = view(CoordArray,:,:,iOld) .+ MoveMultiplier*R


# GridSpan = [0,2*MoveMultiplier*R]
# NCells = Int(floor((GridSpan[2]-GridSpan[1])/(2*1.15)))

# xVec = LinRange(GridSpan[1],GridSpan[2],NCells + 1)
# yVec = LinRange(GridSpan[1],GridSpan[2],NCells + 1)

# l::Float64 = xVec[2]-xVec[1]

# #scatter(CoordArrayShift[1,:,1],CoordArrayShift[2,:,1],size=(900,900))
# #hline!(yVec)
# #vline!(xVec)

# testmat = rand(3,4)
# Testmat2 = reshape(1:24,(3,4,2))

# for (i,v) in pairs(testmat)
#     a = [getindex(i,1),getindex(i,2)]
# end

# @benchmark eachslice(CoordArray[:,:,iOld], dims= 2)


# function Tester1(CoordArray,iOld)
#     sliced = eachslice(CoordArray[:,:,iOld], dims= 2)
#     DivVec::Matrix{Float64} = zeros(2,1000000)
#     DivVecTester = deepcopy(DivVec)
#     a::Vector{Float64} = zeros(1000000)
#     for (i,vec) in enumerate(sliced)
#         @inbounds @view(DivVec[:,i]) .= vec.^2
#         @inbounds @view(a[i]) .= sqrt(sum(@view(DivVec[:,i])[1]))
#     end
#     return a
# end


# @benchmark Tester1(CoordArray,iOld)
# @profview Tester1(CoordArray,iOld)


# eachslice(testmat; dims=1)

# NeighbourCount::Matrix{Int64} = zeros(Int64,NCells,NCells)
# NeighbourCount2::Matrix{Int64} = zeros(Int64,NCells,NCells)

# #@benchmark NeighbourCount .= 0

# Vec::Vector{Float64} = zeros(2)
# DivVec::Vector{Float64} = zeros(2)



# #@time InitNeighbours!(zeros(Int64,2,2),NeighbourCount,zeros(Int64,2),CoordArrayShift,NCells,l)
# #@profview InitNeighbours2!(zeros(Int64,2,2),NeighbourCount2,zeros(Int64,2),CoordArrayShift[:,:,iOld],l,DivVec)
# NeighbourIdx, MaxPart = InitNeighbours(NeighbourCount2,NCells,CoordArrayShift[:,:,iOld],l,DivVec,5)
# CoordArrayShift .-= 50

# UpdateNeighbours!(NeighbourIdx,NeighbourCount2,NCells,CoordArrayShift,l,DivVec,MaxPart,ExtraSlots)
# typeof(CoordArrayShift[:,:,iOld])
# heatmap(NeighbourIdx[7,:,:],size = (1000, 900))
# heatmap(NeighbourCount2,size = (1000, 900))

# ceil((l/0.85)^2)
# @benchmark maximum(NeighbourCount2)

a= abs.(rand(2,5))
bT = collect(enumerate(eachslice(a,dims=2)))


function InitializeSystem(NPart::Int64,ρ::Float64,D::Float64,spread::Vector{Float64},finalMeanRad::Float64,growSize::Float64,NTimeSteps::Int64,Extrasteps::Int64, ExtraSlots::Int64, MoveMultiplier::Float64, 
                            A::Float64,dt::Float64,BrownianMethod!::Function,BoundaryMethod!::Function,b::Float64,save::Bool)

    # ------- Allocate -------

    # Initialize a temporary coordinate array
    CoordArray::Array{Float64} = zeros(2,NPart,2)
    CoordArrayShift::Matrix{Float64} = zeros(2,NPart)

    # Initialize the radius of the particles
    radArray::Vector{Float64} = rand(NPart).*(spread[2]-spread[1]) .+ spread[1] .- growSize

    TempVec::Matrix{Float64} = zeros(2,NPart)

    DeltaDist::Vector{Float64} = zeros(NPart)

    # Defining the radius of the boudary using the packing fraction, mean radius and the number of particles
    R::Float64 = sqrt((NPart*(finalMeanRad)^2)/(ρ))

    # Initializing random angle and random radius arrays
    randAngle::Vector{Float64} = zeros(NPart)
    randRadius::Vector{Float64} = zeros(NPart)

    # Initializing 
    RandDisp::Matrix{Float64} = zeros(2,NPart)

    # Define an old and a new index, these will flip every iteration
    iNew::Int64 = 2
    iOld::Int64 = 1

    # Set the incremental increase in radius equal to the growSize/ NTimeSteps
    IncrementRadIncrease::Float64 = growSize/NTimeSteps
    IncrementDReduce::Float64 = 2*D/(Extrasteps+1)
    
    # Initializing a vector to be filled with radii to each particle
    radCenter::Vector{Float64} = zeros(NPart)
    
    #Initializing the deltacoords matrix
    DeltaVec::Matrix{Float64} = zeros(2,NPart)

    # Initializing force vectors
    forceVec::Matrix{Float64} = zeros(2,NPart)
    BoundaryForce::Matrix{Float64} = zeros(2,NPart)


    GridSpan::Vector{Float64} = [0,2*MoveMultiplier*R]
    NCells::Int64 = Int(floor((GridSpan[2]-GridSpan[1])/(2*spread[2])))

    xVec = LinRange(GridSpan[1],GridSpan[2],NCells + 1)

    l::Float64 = xVec[2]-xVec[1]

    NeighbourCount::Matrix{Int64} = zeros(Int64,NCells,NCells)

    DivVec::Vector{Float64} = zeros(2)

    # Initialize save arrays
    if save
        # Initialize zero arrays
        SaveRadArray = zeros(NPart,floor(Int,(NTimeSteps+Extrasteps)/saveEvery)+1)
        SaveCoordArray = zeros(2,NPart,floor(Int,(NTimeSteps+Extrasteps)/saveEvery)+1)

        # Save the first config
        SaveRadArray[:,1] = radArray
        SaveCoordArray[:,:,1] = CoordArray[:,:,iOld]

        # Set the next saveindex to be 2
        saveCounter = 2
    end
    
    # Fill the coordinate array with initial coordinates given by the InitCoords function.
    InitCoords!(CoordArray,NPart,R,iOld,randAngle,randRadius)

    CoordArrayShift .= @view(CoordArray[:,:,iOld]) .+ MoveMultiplier*R

    NeighbourIndices, MaxPart = InitNeighbours(NeighbourCount,NCells,CoordArrayShift[:,:,iOld],l,DivVec,ExtraSlots)
    CountCoords::CartesianIndices{2, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}} = CartesianIndices(NeighbourCount)
    CountLength::Int64 = length(NeighbourCount)
    CoordArraySliced::Vector{Tuple{Int64, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}} = collect(enumerate(eachslice(CoordArrayShift,dims=2)))
    # Loop over all timesteps AND extrasteps
    for iTimeStep in 1:(NTimeSteps + Extrasteps)
        # Calculate the thermal displacements
        BrownianMethod!(RandDisp, randAngle, randRadius, NPart, D, dt)
        
        # if iTimeStep%300 == 0


        #     CircleArray = zeros(2,100)
        #     PhiCirc = LinRange(0,2*pi,100)
        #     CircleArray[1,:] = R*cos.(PhiCirc)
        #     CircleArray[2,:] = R*sin.(PhiCirc)


        #     function circle(x, y, r=1; n=30)
        #         θ = 0:360÷n:360
        #         Plots.Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
        #     end


        #     # TimeHop = 400
        #     # i = 1
        #     circles = circle.(CoordArray[1,:,iOld],CoordArray[2,:,iOld],radArray[:].*0.5)


        #     plot_kwargs = (aspect_ratio=:equal, fontfamily="Helvetica", legend=false, line="red",
        #         color=:white, grid=false, xlims=(-R-1,R+1), ylims=(-R-1,R+1))


        #     aPlot =  Plots.plot(circles; plot_kwargs...)
        #     Plots.plot!(CircleArray[1,:],CircleArray[2,:],size=(800,800))
        #     display(aPlot)
        # end

        # Sort particles in boxes
        MaxPart, NeighbourIndices = UpdateNeighbours(NeighbourIndices,NeighbourCount,NCells,CoordArrayShift,l,DivVec,MaxPart,ExtraSlots)
        if length(NeighbourCount) != CountLength
            CountCoords = CartesianIndices(NeighbourCount)
            #CountCoordPairs .= collect(pairs(NeighbourCount))
        end
        # Parallel loop over all cells
        CalcElasticRepulsion!(NeighbourCount,NeighbourIndices,forceVec,radArray,CountCoords,DeltaVec,TempVec,DeltaDist,CoordArray[:,:,iOld],NCells)
        # Parallel loop over all particles
        CalcBoundaryRepulsion!(NPart,CoordArray[:,:,iOld],BoundaryMethod!,forceVec,CoordArraySliced,TempVec,R,radCenter,b,BoundaryForce)
        # Multiply this force with A and dt, and add this onto the cordinates of the particle
        CalcNextStep!(CoordArray,iNew,iOld,forceVec,A,dt,RandDisp,NPart)
        # Calculate the shifted coordinates
        CoordArrayShift .= @view(CoordArray[:,:,iNew]) .+ MoveMultiplier*R
        # Increase radius if iTimeStep <= NTimeSteps
        if iTimeStep <= NTimeSteps
            radArray .+= IncrementRadIncrease
        end

        # Decrease brownian noise if iTimeStep > NTimeSteps && iTimeStep <= NTimeSteps+Extrasteps/2
        if iTimeStep > NTimeSteps && iTimeStep <= NTimeSteps+Extrasteps/2
            D -= IncrementDReduce
        end

        # Index Swap
        iOld = 3::Int64 - iOld
        iNew = 3::Int64 - iNew
        
        # Print percentage
        if iTimeStep%1000 == 0
            println("Percentage finished: ", iTimeStep*100/(NTimeSteps+Extrasteps),"%")
        end
        # Save snapshot
        if save && iTimeStep%saveEvery == 0
            @view(SaveRadArray[:,saveCounter]) .= radArray
            @view(SaveCoordArray[:,:,saveCounter]) .= @view(CoordArray[:,:,iOld])
            saveCounter += 1
        end
    end



    # CircleArray = zeros(2,100)
    # PhiCirc = LinRange(0,2*pi,100)
    # CircleArray[1,:] = R*cos.(PhiCirc)
    # CircleArray[2,:] = R*sin.(PhiCirc)


    # function circle(x, y, r=1; n=30)
    #     θ = 0:360÷n:360
    #     Plots.Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
    # end


    # # TimeHop = 400
    # # i = 1
    # circles = circle.(CoordArray[1,:,iOld],CoordArray[2,:,iOld],radArray[:])


    # plot_kwargs = (aspect_ratio=:equal, fontfamily="Helvetica", legend=false, line="red",
    #     color=:black, grid=false)


    # aPlot =  Plots.plot(circles; plot_kwargs...)
    # Plots.plot!(CircleArray[1,:],CircleArray[2,:],size=(800,800))
    # display(aPlot)

    # Return
    if save
        return SaveCoordArray, SaveRadArray, R
    else
        return CoordArray[:,:,iOld], radArray, R
    end
end

# Setting the function variables equal to the functions
if BoundaryMethod == "AttractiveRepellingBoundary"
    BoundaryMethodFunc!::Function = AttractiveRepellingBoundary!
elseif BoundaryMethod =="RepellingBoundary"
    BoundaryMethodFunc!::Function = RepellingBoundary!
end
if BrownianMethod == "BrownianTemperaturePolar"
    BrownianMethodFunc!::Function = BrownianTemperaturePolar!
elseif BrownianMethod =="BrownianTemperatureCartesian"
    BrownianMethodFunc!::Function = BrownianTemperatureCartesian!
end

using BenchmarkTools

# ----------------------------------------------------------------------------------------------------------------------------
# CALC INITIAL VALS #
# ----------------------------------------------------------------------------------------------------------------------------

@time SaveCoordArray, SaveRadArray,R = InitializeSystem(NPart,ρ,D,spread,finalMeanRad,growSize,NTimeSteps,Extrasteps,ExtraSlots,MoveMultiplier,A,dt,BrownianMethodFunc!,BoundaryMethodFunc!,b, false)


# NumberParticles = [10,30,60,100,300,600,1000,3000,6000,10000,20000,30000,40000,50000,60000]
# timeTaken = zeros(length(NumberParticles))
# for iP in 1:length(timeTaken)
#     time = @benchmark InitializeSystem(NumberParticles[$iP],ρ,D,spread,finalMeanRad,growSize,NTimeSteps,Extrasteps,ExtraSlots,MoveMultiplier,A,dt,BrownianMethodFunc!,BoundaryMethodFunc!,b, true)
#     timeTaken[iP] = mean(time).time /10^9
# end
# # ----------------------------------------------------------------------------------------------------------------------------
# aPlot = plot(NumberParticles, timeTaken,xaxis=:lin,yaxis=:lin,label="Time Elapsed",size = (900,800))
# scatter!(NumberParticles, timeTaken,label=false)
# plot!(LinRange(10,60000,1000), 0.0031.*LinRange(10,60000,1000),label="f(x) = 0.0031x")
# ylabel!("Time Elapsed [s]")
# xlabel!("Number of particles [1]")
# timeTaken
# TimePartnumArray =zeros(2,length(timeTaken))
# TimePartnumArray[1,:] .= NumberParticles
# TimePartnumArray[2,:] .= timeTaken
# TimePartnumArray
# using JLD
# save("TimePartnumArray.jld", "TimePartnumArray", TimePartnumArray)

#savefig(aPlot,"TimeplotNumPart.png")
#plot!(LinRange(10,60000,1000), log.(LinRange(10,60000,1000)))
#0.0031*1000000*10/3600


# ----------------------------------------------------------------------------------------------------------------------------
# DELAUNAY #
# ----------------------------------------------------------------------------------------------------------------------------
points = zeros(size(transpose(SaveCoordArray[:,:])))

points += transpose(SaveCoordArray[:,:])

mesh = delaunay(points)

IndexList = 1:NPart
NumNeighbours = zeros(Int,NPart)
for iParticle in 1:NPart
    NumNeighbours[iParticle] = length(IndexList[mesh.vertex_neighbor_vertices[iParticle,:]])
end
NeighbourMatrix = zeros(Int,maximum(NumNeighbours),NPart)
for iParticle in 1:NPart
    NeighbourMatrix[1:NumNeighbours[iParticle],iParticle] = IndexList[mesh.vertex_neighbor_vertices[iParticle,:]]
end
# ----------------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------------
# HDF5 #
# ----------------------------------------------------------------------------------------------------------------------------
FullPath = PathName*FileName

fid = h5open(FullPath,"w")

h5File = create_group(fid, "InitGroup")

h5File["Coords"] = SaveCoordArray[:,:]

h5File["PolAng"] = rand(NPart)*2*pi

h5File["NeighbourMatrix"] = NeighbourMatrix

h5File["R"] = R

h5File["MaxNeighbour"] = maximum(NumNeighbours)

close(fid)
# ----------------------------------------------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------------------------------------------
# PARAMETERFILE #
# ----------------------------------------------------------------------------------------------------------------------------
open(NameListPath,"w") do ParameterFile
    println(ParameterFile, "&Parameters")
    println(ParameterFile, "    NumPart = "*string(NPart))
    println(ParameterFile, "    NumTimeSteps = "*string(NSimulationIterations))
    println(ParameterFile, "    SaveEvery = "*string(SaveEverySimulation))
    println(ParameterFile, "    MaxNeighbour = "*string(maximum(NumNeighbours)))
    println(ParameterFile, "    k = "*string(kSimulation))
    println(ParameterFile, "    kNonLin = "*string(kNonLin))
    println(ParameterFile, "    kBoundary = "*string(kBoundarySimulation))
    println(ParameterFile, "    xi = "*string(xiSimulation))
    println(ParameterFile, "    zeta = "*string(zetaSimulation))
    println(ParameterFile, "    Fa = "*string(FaSimulation))
    println(ParameterFile, "    b = "*string(b))
    println(ParameterFile, "    deltaT = "*string(deltaTSimulation))
    println(ParameterFile, "    IterMethod = "*"\""*IterMethodSimulation*"\"")
    println(ParameterFile, "    BaseDirName = "*"\""*BaseDirName*"\"")
    println(ParameterFile, "    InitFileName = "*"\""*InitFileName*"\"")
    println(ParameterFile, "    SaveFileName = "*"\""*SaveFileName*"\"")
    println(ParameterFile, "    BoundaryMethod = "*"\""*BoundaryMethodSimulation*"\"")
    println(ParameterFile, "    ElasticityMethod = "*"\""*ElasticityMethod*"\"")
    println(ParameterFile, "/")
end
# ----------------------------------------------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------------------------
# Plotting
#--------------------------------------------------------------------------------------------------------
# Plot gif

# Temporary time plot
# plot([100,200,300,500,1000],[40,84,140,226,534])
# scatter!([100,200,300,500,1000],[40,84,140,226,534])
# plot!(LinRange(100,1000,1000),LinRange(100,1000,1000).*0.5.+100)
if PlotEnd
    CircleArray = zeros(2,100)
    PhiCirc = LinRange(0,2*pi,100)
    CircleArray[1,:] = R*cos.(PhiCirc)
    CircleArray[2,:] = R*sin.(PhiCirc)


    function circle(x, y, r=1; n=30)
        θ = 0:360÷n:360
        Plots.Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
    end


    # TimeHop = 400
    # i = 1
    circles = circle.(SaveCoordArray[1,:,end],SaveCoordArray[2,:,end],SaveRadArray[:,end])


    plot_kwargs = (aspect_ratio=:equal, fontfamily="Helvetica", legend=false, line="red",
        color=:black, grid=false)


    aPlot =  Plots.plot(circles; plot_kwargs...)
    Plots.plot!(CircleArray[1,:],CircleArray[2,:],size=(800,800))
    display(aPlot)
end
if CirclePlot
    CircleArray = zeros(2,100)
    PhiCirc = LinRange(0,2*pi,100)
    CircleArray[1,:] = R*cos.(PhiCirc)
    CircleArray[2,:] = R*sin.(PhiCirc)


    function circle(x, y, r=1; n=30)
        θ = 0:360÷n:360
        Plots.Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
    end


    # TimeHop = 400
    # i = 1
    # circles = circle.(SaveCoordArray[1,:],SaveCoordArray[2,:],SaveRadArray[:])


    # plot_kwargs = (aspect_ratio=:equal, fontfamily="Helvetica", legend=false, line="red",
    #     color=:black, grid=false)


    # aPlot =  plot(circles; plot_kwargs...)
    # plot!(CircleArray[1,:],CircleArray[2,:],size=(800,800))
    # display(aPlot)
    TimeHop = 50
    i = 1
    anim = @animate for j ∈ 1:floor(Int,(NTimeSteps+Extrasteps)/500)
        i = 1+(j-1)*500
        circles = circle.(SaveCoordArray[1,:,i],SaveCoordArray[2,:,i],SaveRadArray[:,i])


        plot_kwargs = (aspect_ratio=:equal, fontfamily="Helvetica", legend=false, line="red",
            color=:black, grid=false, xlims=(-R-1,R+1), ylims=(-R-1,R+1))


        aPlot =  Plots.plot(circles; plot_kwargs...)
        Plots.plot!(CircleArray[1,:],CircleArray[2,:],size=(800,800))
        
    end
    gif(anim,"CurrentGif.gif",fps = 15)
end

#--------------------------------------------------------------------------------------------------------

if DelaunayPlot
    using GLMakie
    color = rand(size(mesh.points, 1))*0
    fig, ax, pl = Makie.poly(mesh.points, mesh.simplices, color=color, strokewidth=2, figure=(resolution=(800, 400),))
    display(fig)
end

if NeighbourPlot
    # quite unfinished
    i = 1

    function circle(x, y, r=1; n=30)
        θ = 0:360÷n:360
        Plots.Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
    end

    circles = circle.(SaveCoordArray[1,i,end],SaveCoordArray[2,i,end],SaveRadArray[i,end])

    plot_kwargs = (aspect_ratio=:equal, fontfamily="Helvetica", legend=false, line="red",
        color=:black, grid=false)


    aPlot =  plot(circles; plot_kwargs...)
    plot_kwargs = (aspect_ratio=:equal, fontfamily="Helvetica", legend=false, line=false,
        color=:black, grid=false)

    circles = circle.(SaveCoordArray[1,IndexList[mesh.vertex_neighbor_vertices[i,:]],end],SaveCoordArray[2,IndexList[mesh.vertex_neighbor_vertices[i,:]],end],SaveRadArray[IndexList[mesh.vertex_neighbor_vertices[i,:]],end])

    plot!(circles; plot_kwargs...)

    plot!(CircleArray[1,:],CircleArray[2,:],size=(800,800))
    display(aPlot)
    i += 1

    i = 1
    circles = circle.(SaveCoordArray[1,i,end],SaveCoordArray[2,i,end],SaveRadArray[i,end])


    plot_kwargs = (aspect_ratio=:equal, fontfamily="Helvetica", legend=false, line="red",
        color=:black, grid=false)


    aPlot =  plot(circles; plot_kwargs...)
    plot_kwargs = (aspect_ratio=:equal, fontfamily="Helvetica", legend=false, line=false,
        color=:black, grid=false)

    circles = circle.(SaveCoordArray[1,NeighbourMatrix[1:NumNeighbours[i],i],end],SaveCoordArray[2,NeighbourMatrix[1:NumNeighbours[i],i],end],SaveRadArray[NeighbourMatrix[1:NumNeighbours[i],i],end])

    plot!(circles; plot_kwargs...)

    plot!(CircleArray[1,:],CircleArray[2,:],size=(800,800))
    display(aPlot)
    i +=1
end
