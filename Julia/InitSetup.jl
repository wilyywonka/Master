using Plots
using Delaunay
using HDF5

function InitCoords(NPart,R)
    randAngle = rand(NPart)*2*pi
    randRadius = sqrt.(rand(NPart)*(R^2))
    Coords = zeros(2,NPart)
    Coords[1,:] = cos.(randAngle).*randRadius
    Coords[2,:] = sin.(randAngle).*randRadius
    return Coords
end

function InitializeSystem(NPart,ρ,spread,finalMeanRad,growSize,RandDist,NTimeSteps,Extrasteps,A,dt, saveEvery,save=false)

    R = sqrt((NPart*(finalMeanRad)^2)/(ρ))

    if RandDist == "Uniform"
        radArray = rand(NPart).*(spread[2]-spread[1]) .+ spread[1] .- growSize
    end

    iNew = 2
    iOld = 1


    CoordArray= zeros(2,NPart,2)
    CoordArray[:,:,iOld] = InitCoords(NPart,R)

    if save
        SaveRadArray = zeros(NPart,floor(Int,(NTimeSteps+Extrasteps)/saveEvery)+1)
        SaveCoordArray = zeros(2,NPart,floor(Int,(NTimeSteps+Extrasteps)/saveEvery)+1)

        SaveRadArray[:,1] = radArray
        SaveCoordArray[:,:,1] = CoordArray[:,:,iOld]

        saveCounter = 2
    end

    IndexVec = 1:NPart
    IncrementRadIncrease = (growSize)/NTimeSteps


    for iTimeStep in 1:(NTimeSteps + Extrasteps)
        #Parallel
        Threads.@threads for iParticle in 1:NPart
            CloseParticles = IndexVec[1:end .!= iParticle][findall((CoordArray[1,1:end .!= iParticle,iOld].-CoordArray[1,iParticle,iOld]).^2 .+ (CoordArray[2,1:end .!= iParticle,iOld].-CoordArray[2,iParticle,iOld]).^2 .< (radArray[1:end .!= iParticle].+radArray[iParticle]).^2)]
            forceVec = zeros(2)
            for iNeighbours in CloseParticles
                DeltaVec = CoordArray[:,iParticle,iOld] - CoordArray[:,iNeighbours,iOld]
                forceVec += abs(radArray[iParticle] + radArray[iNeighbours] - sqrt(DeltaVec[1]^2 + DeltaVec[2]^2))^(3/2) .* DeltaVec./sqrt(DeltaVec[1]^2 + DeltaVec[2]^2)
            end
            radCenter = sqrt(CoordArray[1,iParticle,iOld]^2 + CoordArray[2,iParticle,iOld]^2)
            if radCenter > R - radArray[iParticle]
                V = A*(radCenter-(R-radArray[iParticle]))
                forceVec += -(V/radCenter).*[CoordArray[1,iParticle,iOld],CoordArray[2,iParticle,iOld]]
            end
            CoordArray[:,iParticle,iNew] = CoordArray[:,iParticle,iOld] + forceVec*A*dt
        end
        if iTimeStep <= NTimeSteps
            radArray .+= IncrementRadIncrease
        end
        iOld = 3 - iOld
        iNew = 3 - iNew
        if save && iTimeStep%saveEvery == 0
            SaveRadArray[:,saveCounter] = radArray
            SaveCoordArray[:,:,saveCounter] = CoordArray[:,:,iOld]
            saveCounter += 1
        end
    end
    if save
        return SaveCoordArray, SaveRadArray, R
    else
        return CoordArray[:,:,iNew], radArray[:,:,iNew], R[:,iNew]
    end
end


NPart = 500
ρ = 0.90
spread = [0.85,1.15]
finalMeanRad = 1
growSize = 0.3
RandDist = "Uniform"
NTimeSteps = 20000
Extrasteps = 30000 # May need to be even longer
A = 1
dt = 0.02
saveEvery = 1

SaveCoordArray, SaveRadArray, R = InitializeSystem(NPart,ρ,spread,finalMeanRad,growSize,RandDist,NTimeSteps,Extrasteps,A,dt,saveEvery,true)

# Plot
CircleArray = zeros(2,100)
PhiCirc = LinRange(0,2*pi,100)
CircleArray[1,:] = R*cos.(PhiCirc)
CircleArray[2,:] = R*sin.(PhiCirc)


function circle(x, y, r=1; n=30)
    θ = 0:360÷n:360
    Plots.Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
end



TimeHop = 400
i = 1#length(SaveRadArray[1,:])
for j in 1:floor(Int,(NTimeSteps+Extrasteps)/TimeHop)
    println(i)
    circles = circle.(SaveCoordArray[1,:,i],SaveCoordArray[2,:,i],SaveRadArray[:,i])


    plot_kwargs = (aspect_ratio=:equal, fontfamily="Helvetica", legend=false, line="red",
        color=:black, grid=false)


    aPlot =  plot(circles; plot_kwargs...)
    plot!(CircleArray[1,:],CircleArray[2,:],size=(800,800))
    display(aPlot)
    i += TimeHop
end


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




using GLMakie
color = rand(size(mesh.points, 1))*0
fig, ax, pl = Makie.poly(mesh.points, mesh.simplices, color=color, strokewidth=2, figure=(resolution=(800, 400),))
display(fig)

i = 1
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

#TODO
#Lag nabomatrisen og lagre både koordinater og nabomatrisen som HDF5