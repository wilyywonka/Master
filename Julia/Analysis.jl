using Plots
using HDF5

NumIter = 1000
NPart = 700
FaPlot = 0.5

R = h5read("./Fortran/HDF5Files/InitState.h5", "InitGroup/R")

FilePath = "./Fortran/HDF5Files/SaveFiles.h5"

SimulatedCoords = zeros(2, NPart, NumIter)
SimulatedPolAng = zeros(NPart, NumIter)

for i in 1:NumIter
    println("Reading: ", "SaveParams/Coords/"*string(i))
    SimulatedCoords[:,:,i] = h5read(FilePath, "SaveParams/Coords/"*string(i))
    println("Reading: ", "SaveParams/PolAng/"*string(i))
    SimulatedPolAng[:,i]   = h5read(FilePath, "SaveParams/PolAng/"*string(i))
end


CircleArray = zeros(2,100)
PhiCirc = LinRange(0,2*pi,100)
CircleArray[1,:] = R*cos.(PhiCirc)
CircleArray[2,:] = R*sin.(PhiCirc)

iPlot = 1
aPlot = Plots.scatter(SimulatedCoords[1,:,iPlot], SimulatedCoords[2,:,iPlot], size = (900,900), xlims = (-R-3,R+3), ylims = (-R-3,R+3), legend = false)
Plots.quiver!(SimulatedCoords[1,:,iPlot], SimulatedCoords[2,:,iPlot], quiver=(FaPlot.*cos.(SimulatedPolAng[:,iPlot]),FaPlot.*sin.(SimulatedPolAng[:,iPlot])))
Plots.plot!(CircleArray[1,:],CircleArray[2,:])
display(aPlot)
iPlot += 3


SimulationGif = true

if SimulationGif
    Hop = 1
    anim = @animate for j ∈ 1:floor(Int,(NumIter)/Hop)
        iPlot = Hop*j
        aPlot = Plots.scatter(SimulatedCoords[1,:,iPlot], SimulatedCoords[2,:,iPlot], size = (900,900), xlims = (-R-3,R+3), ylims = (-R-3,R+3), legend = false)
        Plots.quiver!(SimulatedCoords[1,:,iPlot], SimulatedCoords[2,:,iPlot], quiver=(FaPlot.*cos.(SimulatedPolAng[:,iPlot]),FaPlot.*sin.(SimulatedPolAng[:,iPlot])))
        Plots.plot!(CircleArray[1,:],CircleArray[2,:]) 
    end
    gif(anim,"SimulatedSystem.gif",fps = 15)
end
