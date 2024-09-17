using Plots
using HDF5

NumIter = 1000
NPart = 100
MarkerSizePlot = 2
ArrowLengthPlot = 0.8

R = h5read("./Fortran/HDF5Files/InitState.h5", "InitGroup/R")

FilePath = "./Fortran/HDF5Files/SaveFiles.h5"

SimulatedCoords = zeros(2, NPart, NumIter)
SimulatedPolAng = zeros(NPart, NumIter)
SimulatedDisplacement = zeros(2, NPart, NumIter)

for i in 1:NumIter
    println("Reading: ", "SaveParams/Coords/"*string(i))
    SimulatedCoords[:,:,i]       = h5read(FilePath, "SaveParams/Coords/"*string(i))
    println("Reading: ", "SaveParams/PolAng/"*string(i))
    SimulatedPolAng[:,i]         = h5read(FilePath, "SaveParams/PolAng/"*string(i))
    println("Reading: ", "SaveParams/Displacement/"*string(i))
    SimulatedDisplacement[:,:,i]   = h5read(FilePath, "SaveParams/Displacement/"*string(i))
end


CircleArray = zeros(2,100)
PhiCirc = LinRange(0,2*pi,100)
CircleArray[1,:] = R*cos.(PhiCirc)
CircleArray[2,:] = R*sin.(PhiCirc)

iPlot = 1
aPlot = Plots.scatter(SimulatedCoords[1,:,iPlot], SimulatedCoords[2,:,iPlot], size = (900,900), xlims = (-R-3,R+3), ylims = (-R-3,R+3), markersize=MarkerSizePlot, color = "black", label="Particle", legend=:topright)
Plots.plot!(SimulatedCoords[1,:,iPlot], SimulatedCoords[2,:,iPlot], quiver=(ArrowLengthPlot.*cos.(SimulatedPolAng[:,iPlot]),ArrowLengthPlot.*sin.(SimulatedPolAng[:,iPlot])), st=:quiver, color="red")
Plots.plot!([R*10, R*10+1],[R*10, R*10+1], color="red", label="Polarization")
Plots.plot!(SimulatedCoords[1,:,iPlot], SimulatedCoords[2,:,iPlot], quiver=(ArrowLengthPlot*SimulatedDisplacement[1,:,iPlot],ArrowLengthPlot*SimulatedDisplacement[2,:,iPlot]), st=:quiver, color="blue")
Plots.plot!([R*10, R*10+1],[R*10, R*10+1], color="blue", label="Displacement")
Plots.plot!(CircleArray[1,:],CircleArray[2,:], label="Boundary", color = "green")
display(aPlot)
iPlot += 2

SimulationGif = true

if SimulationGif
    Hop = 1
    anim = @animate for j ∈ 1:floor(Int,(NumIter)/Hop)
        iPlot = Hop*j
        aPlot = Plots.scatter(SimulatedCoords[1,:,iPlot], SimulatedCoords[2,:,iPlot], size = (900,900), xlims = (-R-3,R+3), ylims = (-R-3,R+3), markersize=MarkerSizePlot, color = "black", label="Particle", legend=:topright)
        Plots.plot!(SimulatedCoords[1,:,iPlot], SimulatedCoords[2,:,iPlot], quiver=(ArrowLengthPlot.*cos.(SimulatedPolAng[:,iPlot]),ArrowLengthPlot.*sin.(SimulatedPolAng[:,iPlot])), st=:quiver, color="red")
        Plots.plot!([R*10, R*10+1],[R*10, R*10+1], color="red", label="Polarization")
        Plots.plot!(SimulatedCoords[1,:,iPlot], SimulatedCoords[2,:,iPlot], quiver=(ArrowLengthPlot*SimulatedDisplacement[1,:,iPlot],ArrowLengthPlot*SimulatedDisplacement[2,:,iPlot]), st=:quiver, color="blue")
        Plots.plot!([R*10, R*10+1],[R*10, R*10+1], color="blue", label="Displacement")
        Plots.plot!(CircleArray[1,:],CircleArray[2,:], label="Boundary", color = "green")
    end
    gif(anim,"SimulatedSystem.gif",fps = 15)
end
