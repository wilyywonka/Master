using GLMakie
using HDF5

NumIter::Int64 = 1000
NPart::Int64 = 1000
MarkerSizePlot::Int64 = 2
ArrowLengthPlot::Float64 = 0.9

FilePath::String = "./HDF5Files/SaveFiles.h5"

NormPlot::Bool = false
DiffPlot::Bool = false
QuivPolPlot::Bool = false
ColPlot::Bool = false

PolOrderAnalysis::Bool = true

if NormPlot|DiffPlot|QuivPolPlot|ColPlot

    R = h5read("./HDF5Files/SaveFiles.h5", "SaveParams/R")


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

end

if NormPlot

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
    # # record(fig,"PolDispAnimationLong.gif",2:5:NumIter,framerate=30) do i
    # #     iPlot[]=i
    # # end
end



function angularChange(OldAng, NewAng)
    return abs.((abs.(OldAng.-NewAng.+pi)).%(2*pi) .-pi)
end

if DiffPlot
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
    # # record(figDiffPol,"PolDiffAnimationLong.gif",2:5:NumIter,framerate=30) do i
    # #     iPlotPol[]=i
    # # end
end



if QuivPolPlot
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
        sleep(0.01)
    end
    # #end
    # # record(figQuivCol,"PolQuivAnimationLong.gif",2:5:NumIter,framerate=30) do i
    # #     iPlotQuivCol[]=i
    # # end
end


if ColPlot
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

    GLMakie.scatter!(axCol, PartPointDat,color=PolAngCol,colorrange =colorrangePol,colormap = :cyclic_grey_15_85_c0_n256,markersize = 15,label="Pol. Turn rate")
    GLMakie.arrows!(PartPointDat,ArrowPolDat,color=PolAngCol,colorrange =colorrangePol,colormap = :cyclic_grey_15_85_c0_n256,arrowsize=5,label="Polarization")

    #GLMakie.arrows!(PartPointDat,ArrowDispCol,color=DispAngCol,colorrange =colorrangePol,colormap = :phase,arrowsize=5,label="Displacement")

    figCol[1, 2] = GLMakie.Legend(figCol,axCol)


    display(figCol)

    for i in 1:NumIter
        iPlotCol[] = i
        sleep(0.0001)
    end
    # record(figCol,"TopDefAnnihilaion.gif",1:5:NumIter,framerate=10) do i
    #     iPlotCol[]=i
    # end
end



if PolOrderAnalysis
    NumParticlesList = [1500,2500,5000,10000,20000]
    NumCycles = 10
    OrderArray = zeros(NumIter,NumCycles,length(NumParticlesList))


    function PolarOrderCalculation(PhiCalc,NPart,NumIter)
        OrderVector = zeros(NumIter)
        for iTimestep in 1:NumIter
            # There will be no time-average, as this is 
            # M = |⟨cosϕ,sinϕ⟩|/NPart = √((Σᵢcosϕᵢ)² + (Σᵢsinϕᵢ)²)/N
            OrderVector[iTimestep] = sqrt((sum(cos.(PhiCalc[:,iTimestep]))^2 + sum(sin.(PhiCalc[:,iTimestep]))^2))/NPart
        end
        return OrderVector
    end

    for iNumPart in eachindex(NumParticlesList)
        
        # Set number of particles
        NPart = NumParticlesList[iNumPart]

        # Loop over number of cycles per number of particles
        for jNumCycles in 1:NumCycles
            
            # Remove the existing saved HDF5-file if it exists
            if isfile("./HDF5Files/SaveFiles.h5")
                run(`rm  ./HDF5Files/SaveFiles.h5`)
            end

            # Make initial system
            @time SaveCoordArray, SaveRadArray,R = InitializeSystem(NPart,ρ,D,spread,finalMeanRad,growSize,NTimeSteps,Extrasteps,ExtraSlots,MoveMultiplier,A,dt,BrownianMethodFunc!,BoundaryMethodFunc!,b, false)

            # Define neighbours
            NumNeighbours, NeighbourMatrix = DelaunayNeighbours(SaveCoordArray, NPart)

            # Write to HDF5
            WriteHDF5(PathName,FileName,SaveCoordArray,NPart,NeighbourMatrix,R,NumNeighbours)

            # Write to ParameterFile
            WriteParameterFile(NameListPath,NPart,NSimulationIterations,SaveEverySimulation,NumNeighbours,kSimulation,kNonLin,kBoundarySimulation,J,xiSimulation,zetaSimulation,FaSimulation,b,deltaTSimulation,
                IterMethodSimulation,BaseDirName,InitFileName,SaveFileName,BoundaryMethodSimulation,ElasticityMethod,PolarizationMethod)

            # Change directory to local run directory
            cd("FortranSerial/build/run/")
            # Run fortran simulation
            run(`time ./ActiveSolidsSimulation`)
            # Change directory back
            cd("../../../")

            # Initialize array to fill with Polarization Angles
            SimulatedPolAng = zeros(NPart, NumIter)

            # Read HDF5 file
            println("Reading: ", "SaveParams/PolAng/")
            for i in 1:NumIter
                SimulatedPolAng[:,i] = h5read(FilePath, "SaveParams/PolAng/"*string(i))
            end

            # Fill OrderArray with 
            OrderArray[:,jNumCycles,iNumPart] .= PolarOrderCalculation(SimulatedPolAng,NPart,NumIter)

            run(`rm  ./HDF5Files/SaveFiles.h5`)

        end

    end

    WriteOrderHDF5("OrderFileLonger.h5", OrderArray, NumParticlesList,NumIter)

    GLMakie.activate!(inline = false,focus_on_show=false)
    OrderFig = Figure(size = (950,800))
    axOrder = Axis(OrderFig[1,1])

    for i in 1:NumCycles
        if i == 1
            GLMakie.lines!(axOrder,1:NumIter,OrderArray[:,1,1],color="blue",label="N=1250")
            GLMakie.lines!(axOrder,1:NumIter,OrderArray[:,1,2],color="red",label="N=2500")
            GLMakie.lines!(axOrder,1:NumIter,OrderArray[:,1,3],color="green",label="N=5000")
            GLMakie.lines!(axOrder,1:NumIter,OrderArray[:,1,4],color="black",label="N=10000")
            OrderFig[1, 2] = GLMakie.Legend(OrderFig,axOrder)
        end
        GLMakie.lines!(axOrder,1:NumIter,OrderArray[:,i,1],color="blue")
        GLMakie.lines!(axOrder,1:NumIter,OrderArray[:,i,2],color="red")
        GLMakie.lines!(axOrder,1:NumIter,OrderArray[:,i,3],color="green")
        GLMakie.lines!(axOrder,1:NumIter,OrderArray[:,i,4],color="black")
    end
    display(OrderFig)
end


function WriteOrderHDF5(Filename, OrderArray, NumPartList,NumIter)
    fid = h5open(Filename,"w")

    h5File = create_group(fid, "OrderGroup")

    h5File["OrderArray"] = OrderArray

    h5File["NumPart"] = NumPartList

    h5File["NumIter"] = collect(1:NumIter)

    close(fid)
end

# WriteOrderHDF5("OrderFile.h5",OrderArray,[1500,2500,5000,10000],NumIter)

summedOrder = zeros(NumIter,length(NumParticlesList))
for i in 1:length(NumParticlesList)
    for j in 1:NumIter
        summedOrder[j,i] = mean(OrderArray[j,:,i])
    end
end


GLMakie.activate!(inline = false,focus_on_show=false)
OrderFigMean = Figure(size = (950,800))
axOrderMean = Axis(OrderFigMean[1,1])


GLMakie.lines!(axOrderMean,1:NumIter,summedOrder[:,1],color="blue",label="N=1250")
GLMakie.lines!(axOrderMean,1:NumIter,summedOrder[:,2],color="red",label="N=2500")
GLMakie.lines!(axOrderMean,1:NumIter,summedOrder[:,3],color="green",label="N=5000")
GLMakie.lines!(axOrderMean,1:NumIter,summedOrder[:,4],color="black",label="N=10000")
OrderFigMean[1, 2] = GLMakie.Legend(OrderFigMean,axOrderMean)

display(OrderFigMean)
# GLMakie.lines(1:NumIter,summedOrder[:,1],color="blue",label="N=1250")
# GLMakie.lines!(1:NumIter,summedOrder[:,2],color="red",label="N=2500")
# GLMakie.lines!(1:NumIter,summedOrder[:,3],color="green",label="N=5000")
# GLMakie.lines!(1:NumIter,summedOrder[:,4],color="black",label="N=10000")
