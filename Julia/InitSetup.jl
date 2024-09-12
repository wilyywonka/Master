using Plots
using Delaunay
using HDF5


# ---------------------------------------------------------------
# SET PARAMETERS - SIMULATION -
NPart = 10000
NSimulationIterations = 100000
SaveEverySimulation = 100
kSimulation = 5
kBoundarySimulation = kSimulation
xiSimulation = 1.0
zetaSimulation = 1.0
FaSimulation = 0.2
b = 0.3
deltaTSimulation = 0.01
IterMethodSimulation = "RK4"
BaseDirName = "/SaveParams"
FileName = "InitState.h5"
InitFileName = "../../HDF5Files/"*FileName
SaveFileName = "../../HDF5Files/SaveFiles.h5"
BoundaryMethod = "RepellingBoundary" #RepellingBoundary : AttractiveRepellingBoundary


# SET PARAMETERS - SAVE -
PathName = "./Fortran/HDF5Files/"
NameListPath = "./Fortran/Parameters/Parameters.nml"


# SET PARAMETERS - INITSTATE -
ρ = 1
spread = [0.85,1.15]
finalMeanRad = 1
growSize = 0.3
D = 0.01
RandDist = "Uniform"
NTimeSteps = 12000
Extrasteps = 40000 # May need to be even longer
A = 5
dt = 0.02
saveEvery = 1
BrownianMethod = "BrownianTemperaturePolar" # BrownianTemperaturePolar : BrownianTemperatureCartesian


# SET PARAMETERS - PLOTTING -
CirclePlot = false
DelaunayPlot = true
NeighbourPlot = false

# ---------------------------------------------------------------

function InitCoords(NPart,R)
    # Define N random angles
    randAngle = rand(NPart)*2*pi
    # Define N random, uniformally distributed values. Multiply these with R², and take the square root.
    # This will ensure a uniformally distributed points, even when using polar coordinates.
    randRadius = sqrt.(rand(NPart)*(R^2))
    # Initialize the coordinate array    
    Coords = zeros(2,NPart)
    # Set the coordinates using the cosine and sine of the angles, times the radius
    Coords[1,:] = cos.(randAngle).*randRadius
    Coords[2,:] = sin.(randAngle).*randRadius
    # Return the coordinates
    return Coords
end

#MAY NEED CHANGES; NO DIVISION/MULTIPLICATION WITH E.G. √3
function BrownianTemperaturePolar(NPart, D, dt)
    # Define N random angles
    randAngle = rand(NPart)*2*pi
    # Define N random, uniformally distributed values and take the square root.
    # This will ensure a uniformally distributed points, even when using polar coordinates.
    randRadius = sqrt.(rand(NPart))
    # Initialize the coordinate array    
    RandDisp = zeros(2,NPart)
    # Set the coordinates using the cosine and sine of the angles, times the radius
    RandDisp[1,:] = cos.(randAngle).*randRadius*sqrt(D*dt)
    RandDisp[2,:] = sin.(randAngle).*randRadius*sqrt(D*dt)
    return RandDisp
end

function BrownianTemperatureCartesian(NPart, D, dt)
    # Initialize the coordinate array    
    RandDisp = zeros(2,NPart)
    # Set the coordinates using the cosine and sine of the angles, times the radius
    RandDisp[1,:] = (2*rand(NPart)-1)*sqrt(D*dt*3)
    RandDisp[2,:] = (2*rand(NPart)-1)*sqrt(D*dt*3)
    return RandDisp
end

function RepellingBoundary(Coords, R, radCenter, b)
    # Initialize a force vector
    BoundaryForce = zeros(2)
    # Check if the CENTER of the particles is outside the circle
    if radCenter > R
        # If so, calculate the force from the boundary
        V = (radCenter-R)
        # And add it to the force vector
        BoundaryForce += -(V/radCenter).*[Coords[1],Coords[2]]
    end
    return BoundaryForce
end

function AttractiveRepellingBoundary(Coords,R,radCenter,b)
    # Initialize a force vector
    BoundaryForce = zeros(2)
    # Check if the CENTER of the particles is within 2b, the attractive/repulsive zone
    if (R-radCenter) < 2*b
        # If so, calculate the force from the boundary
        V = (R-radCenter)^2 -2*(R-radCenter)*b
        # And add it to the force vector
        BoundaryForce += -(V/radCenter).*[Coords[1],Coords[2]]
    end
    return BoundaryForce
end

function InitializeSystem(NPart,ρ,D,spread,finalMeanRad,growSize,RandDist,NTimeSteps,Extrasteps,A,dt, saveEvery, BrownianMethod, BoundaryMethod, b,save=false)

    #Defining the radius of the boudary using the packing fraction, mean radius and the number of particles
    R = sqrt((NPart*(finalMeanRad)^2)/(ρ))

    # Initialize the radius of the particles using a set distribution
    # Can be expanded to include e.g. a truncated normal distribution if needed.
    if RandDist == "Uniform"
        radArray = rand(NPart).*(spread[2]-spread[1]) .+ spread[1] .- growSize
    end

    # Define an old and a new index, these will flip every iteration
    iNew = 2
    iOld = 1

    # Initialize a temporary coordinate array
    CoordArray= zeros(2,NPart,2)
    # Fill the coordinate array with initial coordinates given by the InitCoords function.
    CoordArray[:,:,iOld] = InitCoords(NPart,R)

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

    # Define an indexvector to be used later to find overlapping particles
    IndexVec = 1:NPart
    # Set the incremental increase in radius equal to the growSize/ NTimeSteps
    IncrementRadIncrease = growSize/NTimeSteps
    IncrementDReduce = 2*D/(Extrasteps)

    # Loop over all timesteps AND extrasteps
    for iTimeStep in 1:(NTimeSteps + Extrasteps)
        #Parallel loop over all particles
        RandDisp = BrownianMethod(NPart,D,dt)
        Threads.@threads for iParticle in 1:NPart
            # Calc distance to all other particles and see if any are inside r_1 + r_2
            CloseParticles = IndexVec[1:end .!= iParticle][findall((CoordArray[1,1:end .!= iParticle,iOld].-CoordArray[1,iParticle,iOld]).^2 .+ (CoordArray[2,1:end .!= iParticle,iOld].-CoordArray[2,iParticle,iOld]).^2 .< (radArray[1:end .!= iParticle].+radArray[iParticle]).^2)]
            # Empty vector to be used when adding the forces in x and y
            forceVec = zeros(2)
            # Loop over the neighbours
            for iNeighbours in CloseParticles
                # Calculate the distance difference
                DeltaVec = CoordArray[:,iParticle,iOld] - CoordArray[:,iNeighbours,iOld]
                # Add the force in the opposite direction to the force vector
                forceVec += abs(radArray[iParticle] + radArray[iNeighbours] - sqrt(DeltaVec[1]^2 + DeltaVec[2]^2))^(3/2) .* DeltaVec./sqrt(DeltaVec[1]^2 + DeltaVec[2]^2)
            end
            # Calculate the distance from the center of the circle
            radCenter = sqrt(CoordArray[1,iParticle,iOld]^2 + CoordArray[2,iParticle,iOld]^2)

            forceVec += BoundaryMethod(CoordArray[:,iParticle,iOld], R, radCenter, b)

            # Multiply this force with A and dt, and add this onto the cordinates of the particle
            CoordArray[:,iParticle,iNew] = CoordArray[:,iParticle,iOld] + forceVec*A*dt + RandDisp[:,iParticle]
        end

        # Increase radius if iTimeStep <= NTimeSteps
        if iTimeStep <= NTimeSteps
            radArray .+= IncrementRadIncrease
        end

        if iTimeStep > NTimeSteps && iTimeStep <= NTimeSteps+Extrasteps/2
            D -= IncrementDReduce
        end

        # Index Swap
        iOld = 3 - iOld
        iNew = 3 - iNew
        
        if iTimeStep%1000 == 0
            println("Percentage finished: ", iTimeStep*100/(NTimeSteps+Extrasteps),"%")
        end
        # Save config to array
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
end

# Setting the function variables equal to the functions
if BoundaryMethod == "AttractiveRepellingBoundary"
    BoundaryMethodFunc = AttractiveRepellingBoundary
elseif BoundaryMethod =="RepellingBoundary"
    BoundaryMethodFunc = RepellingBoundary
end
if BrownianMethod == "BrownianTemperaturePolar"
    BrownianMethodFunc = BrownianTemperaturePolar
elseif BrownianMethod =="BrownianTemperatureCartesian"
    BrownianMethodFunc = BrownianTemperatureCartesian
end

# ----------------------------------------------------------------------------------------------------------------------------
# CALC INITIAL VALS #
# ----------------------------------------------------------------------------------------------------------------------------
@time SaveCoordArray, SaveRadArray, R = InitializeSystem(NPart,ρ,D,spread,finalMeanRad,growSize,RandDist,NTimeSteps,Extrasteps,A,dt,saveEvery,BrownianMethodFunc,BoundaryMethodFunc,b,false)
# ----------------------------------------------------------------------------------------------------------------------------


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
    # i = 1#length(SaveRadArray[1,:])
    # for j in 1:floor(Int,(NTimeSteps+Extrasteps)/TimeHop)
    #     println(i)
    #     circles = circle.(SaveCoordArray[1,:,i],SaveCoordArray[2,:,i],SaveRadArray[:,i])


    #     plot_kwargs = (aspect_ratio=:equal, fontfamily="Helvetica", legend=false, line="red",
    #         color=:black, grid=false)


    #     aPlot =  plot(circles; plot_kwargs...)
    #     plot!(CircleArray[1,:],CircleArray[2,:],size=(800,800))
    #     display(aPlot)
    #     i += TimeHop
    # end
    TimeHop = 50
    i = 1
    anim = @animate for j ∈ 1:floor(Int,(NTimeSteps+Extrasteps)/500)
        i = 1+(j-1)*500
        circles = circle.(SaveCoordArray[1,:,i],SaveCoordArray[2,:,i],SaveRadArray[:,i])


        plot_kwargs = (aspect_ratio=:equal, fontfamily="Helvetica", legend=false, line="red",
            color=:black, grid=false, xlims=(-R-1,R+1), ylims=(-R-1,R+1))


        aPlot =  plot(circles; plot_kwargs...)
        plot!(CircleArray[1,:],CircleArray[2,:],size=(800,800))
        
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
