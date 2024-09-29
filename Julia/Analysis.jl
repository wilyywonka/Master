using GLMakie
using HDF5

NumIter = 2250
NPart = 50000
MarkerSizePlot = 2
ArrowLengthPlot = 0.9

R = h5read("./Fortran/HDF5Files/SaveFiles.h5", "SaveParams/R")

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

# seconds = 0:0.1:2
# measurements = [8.2, 8.4, 6.3, 9.5, 9.1, 10.5, 8.6, 8.2, 10.5, 8.5, 7.2,
#         8.8, 9.7, 10.8, 12.5, 11.6, 12.1, 12.1, 15.1, 14.7, 13.1]

#GLMakie.activate!(; screen_config...)


GLMakie.activate!(inline = false,focus_on_show=false)
fig = Figure(size = (950,800))
iPlot = Observable(1)

PartPointDat = @lift Point2f.(SimulatedCoords[1,:,$iPlot], SimulatedCoords[2,:,$iPlot])
ArrowDispDat = @lift Vec2f.(ArrowLengthPlot*SimulatedDisplacement[1,:,$iPlot],ArrowLengthPlot*SimulatedDisplacement[2,:,$iPlot])
ArrowPolDat = @lift Vec2f.(ArrowLengthPlot.*cos.(SimulatedPolAng[:,$iPlot]), ArrowLengthPlot.*sin.(SimulatedPolAng[:,$iPlot]))

ax = Axis(fig[1, 1])



GLMakie.scatter!(ax, PartPointDat, markersize = 2, color="black")
GLMakie.lines!(ax,CircleArray[1,:],CircleArray[2,:],color="green",label="Border")
GLMakie.arrows!(PartPointDat,ArrowPolDat,color="blue",arrowsize=5,label="Polarization")
GLMakie.arrows!(PartPointDat,ArrowDispDat,color="red",arrowsize=5,label="Displacement")
fig[1, 2] = GLMakie.Legend(fig,ax)
display(fig)


for i in 1:NumIter
    iPlot[] = i
    sleep(0.01)
end


# record(fig,"PolDispAnimationLong.gif",2:5:NumIter,framerate=30) do i
#     iPlot[]=i
# end



function angularChange(OldAng, NewAng)
    return abs.((abs.(OldAng.-NewAng.+pi)).%(2*pi) .-pi)
end

figDiffPol = Figure(size = (950,850))
iPlotPol = Observable(2)

PartPointDat = @lift Point2f.(SimulatedCoords[1,:,$iPlotPol], SimulatedCoords[2,:,$iPlotPol])
PartPolDiff = @lift angularChange(SimulatedPolAng[:,($iPlotPol-1)],SimulatedPolAng[:,$iPlotPol])


cmap = :matter

ax = Axis(figDiffPol[1, 1])

GLMakie.scatter!(ax, PartPointDat,color=PartPolDiff,colormap = cmap,markersize = 5)
GLMakie.lines!(ax,CircleArray[1,:],CircleArray[2,:],color="green",label="Border")
#figDiffPol[1, 2] = GLMakie.Legend(figDiffPol,ax)
Colorbar(figDiffPol[1, 2],colormap = cmap)

display(figDiffPol)


for i in 2:NumIter
    iPlotPol[] = i
    sleep(0.01)
end


# record(figDiffPol,"PolDiffAnimationLong.gif",2:5:NumIter,framerate=30) do i
#     iPlotPol[]=i
# end


figQuivCol = Figure(size = (950,800))
iPlotQuivCol = Observable(2)

PartPointDat = @lift Point2f.(SimulatedCoords[1,:,$iPlotQuivCol], SimulatedCoords[2,:,$iPlotQuivCol])
ArrowPolDat = @lift Vec2f.(ArrowLengthPlot.*cos.(SimulatedPolAng[:,$iPlotQuivCol]), ArrowLengthPlot.*sin.(SimulatedPolAng[:,$iPlotQuivCol]))
PolAngCol = @lift SimulatedPolAng[:,$iPlotQuivCol].%(2*pi)
PartPolDiff = @lift angularChange(SimulatedPolAng[:,($iPlotQuivCol-1)],SimulatedPolAng[:,$iPlotQuivCol])
ArrowDispDat = @lift Vec2f.(ArrowLengthPlot*SimulatedDisplacement[1,:,$iPlotQuivCol],ArrowLengthPlot*SimulatedDisplacement[2,:,$iPlotQuivCol])
axQuivCol = Axis(figQuivCol[1, 1])

GLMakie.scatter!(axQuivCol, PartPointDat,color=PartPolDiff,colormap = cmap,markersize = 7,label="Pol. Turn rate")
GLMakie.lines!(axQuivCol,CircleArray[1,:],CircleArray[2,:],color="green",label="Border")
GLMakie.arrows!(PartPointDat,ArrowPolDat,color ="blue",arrowsize=5,label="Polarization")
GLMakie.arrows!(PartPointDat,ArrowDispDat,color="red",arrowsize=5,label="Displacement")

figQuivCol[1, 2] = GLMakie.Legend(figQuivCol,axQuivCol)


display(figQuivCol)

for i in 2:NumIter
    iPlotQuivCol[] = i
    sleep(0.05)
end

# record(figQuivCol,"PolQuivAnimationLong.gif",2:5:NumIter,framerate=30) do i
#     iPlotQuivCol[]=i
# end


figCol = Figure(size = (1300,1100))
iPlotCol = Observable(1)

PartPointDat = @lift Point2f.(SimulatedCoords[1,:,$iPlotCol], SimulatedCoords[2,:,$iPlotCol])
ArrowPolDat = @lift Vec2f.(ArrowLengthPlot.*cos.(SimulatedPolAng[:,$iPlotCol]), ArrowLengthPlot.*sin.(SimulatedPolAng[:,$iPlotCol]))
PolAngCol = @lift SimulatedPolAng[:,$iPlotCol].%(2*pi)
ArrowDispCol = @lift Vec2f.(ArrowLengthPlot*SimulatedDisplacement[1,:,$iPlotCol],ArrowLengthPlot*SimulatedDisplacement[2,:,$iPlotCol])
DispAngCol = @lift atan.(SimulatedDisplacement[2,:,$iPlotCol],SimulatedDisplacement[1,:,$iPlotCol])

axCol = Axis(figCol[1, 1])



colorrangePol = [0, 2*pi]

GLMakie.lines!(axCol,CircleArray[1,:],CircleArray[2,:],color="green",label="Border")
#GLMakie.arrows!(PartPointDat,ArrowPolDat,color=PolAngCol,colorrange =colorrangePol,colormap = :cyclic_grey_15_85_c0_n256,arrowsize=5,label="Polarization")
GLMakie.arrows!(PartPointDat,ArrowPolDat,color=PolAngCol,colorrange =colorrangePol,colormap = :phase,arrowsize=5,label="Polarization")
#GLMakie.arrows!(PartPointDat,ArrowDispCol,color=DispAngCol,colorrange =colorrangePol,colormap = :phase,arrowsize=5,label="Displacement")

figCol[1, 2] = GLMakie.Legend(figCol,axCol)


display(figCol)

for i in 1:NumIter
    iPlotCol[] = i
    sleep(0.01)
end

# record(figCol,"TopDefAnnihilaion.gif",1:5:NumIter,framerate=10) do i
#     iPlotCol[]=i
# end


SimulationGif = true


# fig = Figure()
# # --- figure set-up: two panes side by side, right one with two sub-panes: 
# # --- left pane:
# ax = Axis(fig[1, 1]) 
# GLMakie.ylims!(ax, 0, 30)

# # --- right pane contains two sub-panes: [1,1]: top, [2,1]: bottom
# # --- SliderGrid: right pane, top sub-pane:
# sg = SliderGrid(  fig[1, 2][1, 1], 
#     (; label = "Voltage",     range = 0:0.1:10, format = "{:.1f}V", startvalue = 5.3),
#     (; label = "Current",     range = 0:0.1:20, format = "{:.1f}A", startvalue = 10.2),
#     (; label = "Resistance",  range = 0:0.1:50, format = "{:.1f}Ω", startvalue = 15.9),
#     width = 350,
#     tellheight = false,
#     )

# display(fig)

# function arrow0M!(x, y, u, v;lc=:black, la=2)
#     #for i in eachindex(x)
#     xMat = zeros(2,length(x))
#     yMat = zeros(2,length(x))
#     xMat[1,:] .= x
#     xMat[2,:] .= x.+u
#     yMat[1,:] .= y
#     yMat[2,:] .= y.+v

#     plot!(xMat, yMat, lc=lc,la=la,label=false)
#     #end
# end

# a = zeros(5)
# bT = collect(1:5)
# [a,a.+bT]'

# bTest = zeros(5,2)

# bTest[:,1] .= a
# bTest[:,2] .= a.+ bT


# Matrix([a,a.+bT])

# iPlot = 180
# aPlot = Plots.scatter(SimulatedCoords[1,:,iPlot], SimulatedCoords[2,:,iPlot], size = (1200,1200), xlims = (-R-3,R+3), ylims = (-R-3,R+3), markersize=MarkerSizePlot, color = "black", label="Particle", legend=:topright)
# # Plots.plot!(SimulatedCoords[1,:,iPlot], SimulatedCoords[2,:,iPlot], quiver=(ArrowLengthPlot.*cos.(SimulatedPolAng[:,iPlot]),ArrowLengthPlot.*sin.(SimulatedPolAng[:,iPlot])), st=:quiver, color="blue",arrow=arrow(0.00, 0.0))
# arrow0M!(SimulatedCoords[1,:,iPlot], SimulatedCoords[2,:,iPlot],ArrowLengthPlot.*cos.(SimulatedPolAng[:,iPlot]),ArrowLengthPlot.*sin.(SimulatedPolAng[:,iPlot]); lc=:blue, la=1)
# Plots.plot!([R*10, R*10+1],[R*10, R*10+1], color="blue", label="Polarization")
# # Plots.plot!(SimulatedCoords[1,:,iPlot], SimulatedCoords[2,:,iPlot], quiver=(ArrowLengthPlot*SimulatedDisplacement[1,:,iPlot],ArrowLengthPlot*SimulatedDisplacement[2,:,iPlot]), st=:quiver, color="red")
# # arrows!(SimulatedCoords[1,:,iPlot], SimulatedCoords[2,:,iPlot],ArrowLengthPlot*SimulatedDisplacement[1,:,iPlot],ArrowLengthPlot*SimulatedDisplacement[2,:,iPlot])

# arrow0M!(SimulatedCoords[1,:,iPlot], SimulatedCoords[2,:,iPlot],ArrowLengthPlot*SimulatedDisplacement[1,:,iPlot],ArrowLengthPlot*SimulatedDisplacement[2,:,iPlot]; lc=:red, la=1)
# Plots.plot!([R*10, R*10+1],[R*10, R*10+1], color="red", label="Displacement")
# Plots.plot!(CircleArray[1,:],CircleArray[2,:], label="Boundary", color = "green")
# display(aPlot)
# iPlot += 10

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
