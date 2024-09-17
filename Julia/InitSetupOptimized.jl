using Plots
using Delaunay
using HDF5


# ---------------------------------------------------------------
# SET PARAMETERS - SIMULATION -
NPart::Int64 = 100
NSimulationIterations::Int64 = 2000
SaveEverySimulation::Int64 = 2
kSimulation::Float64 = 5.0
kBoundarySimulation::Float64 = kSimulation
xiSimulation::Float64 = 1.0
zetaSimulation::Float64 = 1.0
FaSimulation::Float64 = 1
b::Float64 = 0.3
deltaTSimulation::Float64 = 0.01
IterMethodSimulation::String = "RK4"
BaseDirName::String = "/SaveParams"
FileName::String = "InitState.h5"
InitFileName::String = "../../HDF5Files/"*FileName
SaveFileName::String = "../../HDF5Files/SaveFiles.h5"
BoundaryMethod::String = "RepellingBoundary" #RepellingBoundary : AttractiveRepellingBoundary


# SET PARAMETERS - SAVE -
PathName::String = "./Fortran/HDF5Files/"
NameListPath::String = "./Fortran/Parameters/Parameters.nml"


# SET PARAMETERS - INITSTATE -
ρ::Float64 = 1.0
spread::Vector{Float64} = [0.85,1.15]
finalMeanRad::Float64 = 1.0
growSize::Float64 = 0.3
D::Float64 = 0.01
NTimeSteps::Int64 = 12000
Extrasteps::Int64 = 40000 # May need to be even longer
A::Float64 = 5.0
dt::Float64 = 0.02
saveEvery::Int64 = 1
BrownianMethod::String = "BrownianTemperaturePolar" # BrownianTemperaturePolar : BrownianTemperatureCartesian


# SET PARAMETERS - PLOTTING -
CirclePlot::Bool = false
DelaunayPlot::Bool = true
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
    RandDisp[1,:] = cos.(randAngle).*randRadius*sqrt(D*dt)
    RandDisp[2,:] = sin.(randAngle).*randRadius*sqrt(D*dt)
    return nothing
end

function BrownianTemperatureCartesian!(RandDisp::Matrix{Float64}, NPart::Int64, D::Float64, dt::Float64)
    # Set the coordinates using the cosine and sine of the angles, times the radius
    RandDisp[1,:] = (2*rand(NPart)-1)*sqrt(D*dt*3)
    RandDisp[2,:] = (2*rand(NPart)-1)*sqrt(D*dt*3)
    return nothing
end

function RepellingBoundary!(Coords::Array{Float64}, BoundaryForce::Matrix{Float64}, R::Float64, radCenter::Vector{Float64}, b::Float64, iParticle::Int64, iOld::Int64)
    # Initialize a force vector
    BoundaryForce[:,iParticle] .= 0.0
    # Check if the CENTER of the particles is outside the circle
    if radCenter[iParticle] > R
        # Add it to the force vector
        BoundaryForce[:,iParticle] += -(((radCenter[iParticle]-R))/radCenter[iParticle]).*[Coords[1,iParticle,iOld],Coords[2,iParticle,iOld]]
    end
    return nothing
end

function AttractiveRepellingBoundary!(Coords::Array{Float64}, BoundaryForce::Matrix{Float64}, R::Float64, radCenter::Vector{Float64}, b::Float64, iParticle::Int64, iOld::Int64)
    # Initialize a force vector
    BoundaryForce[:,iParticle] .= 0.0
    # Check if the CENTER of the particles is within 2b, the attractive/repulsive zone
    if (R-radCenter[iParticle]) < 2*b
        # Add it to the force vector
        BoundaryForce[:,iParticle] += -(((R-radCenter[iParticle])^2 -2*(R-radCenter[iParticle])*b)/radCenter[iParticle]).*[Coords[1,iParticle,iOld],Coords[2,iParticle,iOld]]
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



function InitializeSystem(NPart::Int64,ρ::Float64,D::Float64,spread::Vector{Float64},finalMeanRad::Float64,growSize::Float64,NTimeSteps::Int64,Extrasteps::Int64,A::Float64,dt::Float64,BrownianMethod!::Function,BoundaryMethod!::Function,b::Float64,save::Bool)

    # ------- Allocate -------

    # Initialize a temporary coordinate array
    CoordArray::Array{Float64} = zeros(2,NPart,2)

    # Initialize the radius of the particles
    radArray::Vector{Float64} = rand(NPart).*(spread[2]-spread[1]) .+ spread[1] .- growSize

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

    # Define an indexvector to be used later to find overlapping particles
    IndexVec::Vector{Int64} = Vector(1:NPart)

    # Set the incremental increase in radius equal to the growSize/ NTimeSteps
    IncrementRadIncrease::Float64 = growSize/NTimeSteps
    IncrementDReduce::Float64 = 2*D/(Extrasteps)
    
    # Initializing a vector to be filled with radii to each particle
    radCenter::Vector{Float64} = zeros(NPart)
    
    #Initializing the deltacoords matrix
    DeltaVec::Matrix{Float64} = zeros(2,NPart)

    # Initializing force vectors
    forceVec::Matrix{Float64} = zeros(2,NPart)
    BoundaryForce::Matrix{Float64} = zeros(2,NPart)

    CloseParticlesArray::Vector{Vector{Int64}} = []

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

    InitCloseParticles!(CloseParticlesArray, NPart)

    # Loop over all timesteps AND extrasteps
    for iTimeStep in 1:(NTimeSteps + Extrasteps)
        # Calculate the thermal displacements
        BrownianMethod!(RandDisp, randAngle, randRadius, NPart, D, dt)
        # Parallel loop over all particles
        Threads.@threads for iParticle in 1:NPart
            # Calc distance to all other particles and see if any are inside r_1 + r_2
            CloseParticles!(CloseParticlesArray, IndexVec, CoordArray, iOld, radArray, iParticle)
            # Empty vector to be used when adding the forces in x and y
            forceVec[:,iParticle] .= 0.0
            # Loop over the neighbours
            for iNeighbours in CloseParticlesArray[iParticle]
                # Calculate the distance difference
                DeltaVec[:,iParticle] = CoordArray[:,iParticle,iOld] - CoordArray[:,iNeighbours,iOld]
                # Add the force in the opposite direction to the force vector
                forceVec[:,iParticle] += abs(radArray[iParticle] + radArray[iNeighbours] - sqrt(DeltaVec[1,iParticle]^2 + DeltaVec[2,iParticle]^2))^(3/2) .* DeltaVec[:,iParticle]./sqrt(DeltaVec[1,iParticle]^2 + DeltaVec[2,iParticle]^2)
            end
            # Calculate the distance from the center of the circle
            radCenter[iParticle] = sqrt(CoordArray[1,iParticle,iOld]^2 + CoordArray[2,iParticle,iOld]^2)
        
            BoundaryMethod!(CoordArray, BoundaryForce, R, radCenter, b, iParticle, iOld)
            forceVec[:,iParticle] += BoundaryForce[:,iParticle]
            # Multiply this force with A and dt, and add this onto the cordinates of the particle
            CoordArray[:,iParticle,iNew] = CoordArray[:,iParticle,iOld] + forceVec[:,iParticle]*A*dt + RandDisp[:,iParticle]
        end
        # Increase radius if iTimeStep <= NTimeSteps
        if iTimeStep <= NTimeSteps
            radArray .+= IncrementRadIncrease
        end

        if iTimeStep > NTimeSteps && iTimeStep <= NTimeSteps+Extrasteps/2
            D -= IncrementDReduce
        end

        # Index Swap
        iOld = 3::Int64 - iOld
        iNew = 3::Int64 - iNew
        
        if iTimeStep%1000 == 0
            println("Percentage finished: ", iTimeStep*100/(NTimeSteps+Extrasteps),"%")
        end
        if save && iTimeStep%saveEvery == 0
            SaveRadArray[:,saveCounter] = radArray
            SaveCoordArray[:,:,saveCounter] = CoordArray[:,:,iOld]
            saveCounter += 1
        end
    end

    # Return
    if save
        return SaveCoordArray, SaveRadArray, R
    else
        return CoordArray[:,:,iNew], radArray, R
    end
#    return CoordArray[:,:,iNew], radArray, R
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

@time SaveCoordArray, SaveRadArray, R =  InitializeSystem(NPart,ρ,D,spread,finalMeanRad,growSize,NTimeSteps,Extrasteps,A,dt,BrownianMethodFunc!,BoundaryMethodFunc!,b, true)

# ----------------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------------
# DELAUNAY #
# ----------------------------------------------------------------------------------------------------------------------------
points = zeros(size(transpose(SaveCoordArray[:,:,end])))

points += transpose(SaveCoordArray[:,:,end])

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

h5File["Coords"] = SaveCoordArray[:,:,end]

h5File["PolAng"] = rand(NPart)*2*pi

h5File["NeighbourMatrix"] = NeighbourMatrix

h5File["R"] = R

close(fid)
# ----------------------------------------------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------------------------------------------
# PARAMETERFILE #
# ----------------------------------------------------------------------------------------------------------------------------
open(NameListPath,"w") do ParameterFile
    println(ParameterFile, "&Parameters")
    println(ParameterFile, "    NumPart = "*string(NPart))
    println(ParameterFile, "    NumTimeSteps = "*string(NSimulationIterations))
    println(ParameterFile, "    MaxNeighbour = "*string(maximum(NumNeighbours)))
    println(ParameterFile, "    SaveEvery = "*string(SaveEverySimulation))
    println(ParameterFile, "    k = "*string(kSimulation))
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
    println(ParameterFile, "    BoundaryMethod = "*"\""*BoundaryMethod*"\"")
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

# CircleArray = zeros(2,100)
# PhiCirc = LinRange(0,2*pi,100)
# CircleArray[1,:] = R*cos.(PhiCirc)
# CircleArray[2,:] = R*sin.(PhiCirc)


# function circle(x, y, r=1; n=30)
#     θ = 0:360÷n:360
#     Plots.Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
# end


# TimeHop = 400
# i = 1
# circles = circle.(SaveCoordArray[1,:,end],SaveCoordArray[2,:,end],SaveRadArray[:,end])


# plot_kwargs = (aspect_ratio=:equal, fontfamily="Helvetica", legend=false, line="red",
#     color=:black, grid=false)


# aPlot =  Plots.plot(circles; plot_kwargs...)
# Plots.plot!(CircleArray[1,:],CircleArray[2,:],size=(800,800))
# display(aPlot)

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
