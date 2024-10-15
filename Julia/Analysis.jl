using GLMakie
using Plots
using Statistics
using HDF5

NumIter::Int64 = 4000
NPart::Int64 = 40000
MarkerSizePlot::Int64 = 2
ArrowLengthPlot::Float64 = 0.9

FilePath::String = "./HDF5Files/SaveFiles.h5"
H5FileName = "HDF5Files/CorrelationFileN40000.h5"

NormPlot::Bool = false
DiffPlot::Bool = false
QuivPolPlot::Bool = false
ColPlot::Bool = false

PolOrderAnalysis::Bool = false

NaiveCorrelationAnalysis::Bool = false
SimpleCorrelation::Bool = true

CorrLengthAnalysis::Bool = false
saveCorrelationFunction::Bool = false
CorrelationPlot::Bool = false

if NormPlot|DiffPlot|QuivPolPlot|ColPlot|SimpleCorrelation|NaiveCorrelationAnalysis

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

NumPartBases::Int64 = 1000
PartHop::Int64 = ceil(Int64, NPart/NumPartBases)

if length(collect(1:PartHop:NPart)) != NumPartBases
    error("Number of particle-bases is wrong. Expected number is: ", NumPartBases, "; got: ",length(collect(1:PartHop:NPart)), "\n Please choose a number of particle-bases divisible by number of particles!")
end
    

ΔBins::Float64 = 2.6 # 2*biggest radius, plus a little
NumBins::Int64 = ceil(Int64, 2*R/ΔBins + 10)
LBins::Float64 = NumBins*ΔBins


ExtraSlots::Int64 = 5
MoveMultiplier::Float64 = 1.5
PartSizeSpread::Vector{Float64} = [0.85,1.15]
SimulatedCoordsShift::Matrix{Float64} = zeros(2,NPart)


if SimpleCorrelation

    function AnalyzeCorrelationSimple(NumIter::Int64, NPart::Int64, ΔBins::Float64, NumBins::Int64, NumPartBases::Int64, 
        PartHop::Int64, SimulatedCoords::Array{Float64}, SimulatedPolAng::Matrix{Float64})

        CorrelationFunctionValues::Array{Float64} = zeros(NumBins,NumPartBases,NumIter)
        CorrelationFunctionMean::Matrix{Float64} = zeros(NumBins,NumIter)
        
        BinMidpointRadius = collect(1:NumBins)*ΔBins .- ΔBins/2
        NumParticleInBin::Array{Int64} = zeros(Int64,NumBins,NumPartBases,NumIter)
        PartBases::Vector{Int64} = collect(1:PartHop:NPart)
        RadPartMat::Matrix{Float64} = zeros(NPart,NumPartBases)
        CurrBinMat::Matrix{Int64} = zeros(NPart-1,NumPartBases)
        ParticleNumbers::Matrix{Int64} = zeros(Int64,NPart-1, NumPartBases)
        CounterVec::Vector{Int64} = ones(Int64,NumPartBases)
        PartBinSummed::Matrix{Int64} = zeros(Int64, NumBins, NumIter)

        for iIter in 1:NumIter
            #Parallel
            #Threads.@threads for iNumBase in 1:NumPartBases
            Threads.@threads for iNumBase in 1:NumPartBases

                iParticle::Int64 = PartBases[iNumBase]

                @view(CounterVec[iNumBase]) .= 1

                for jPart in 1:NPart
                    if jPart != iParticle
                        @view(ParticleNumbers[CounterVec[iNumBase],iNumBase]) .= jPart
                        @view(CounterVec[iNumBase]) .+= 1
                    end
                end

                for jNeighbour in 1:NPart-1
                    jParticle::Int64 = ParticleNumbers[jNeighbour,iNumBase]
                    @inbounds @view(RadPartMat[jNeighbour,iNumBase]) .= sqrt((SimulatedCoords[1,jParticle,iIter]-SimulatedCoords[1,iParticle,iIter]).^2 .+ (SimulatedCoords[2,iParticle,iIter]-SimulatedCoords[2,jParticle,iIter]).^2)
                    @inbounds @view(CurrBinMat[jNeighbour,iNumBase]) .= ceil(Int64,RadPartMat[jNeighbour,iNumBase]./ΔBins)
                    @inbounds @view(NumParticleInBin[CurrBinMat[jNeighbour,iNumBase],iNumBase,iIter]) .+= 1
                    @inbounds @view(CorrelationFunctionValues[CurrBinMat[jNeighbour,iNumBase],iNumBase,iIter]) .+= cos(SimulatedPolAng[iParticle,iIter]-SimulatedPolAng[jParticle,iIter])
                end
            end
            if iIter%100 == 0
                println("Iteration: ",iIter)
            end
        end

        for iIter in 1:NumIter
            for iBin in 1:NumBins
                SummedNumber = sum(NumParticleInBin[iBin, :, iIter])
                PartBinSummed[iBin,iIter] = round(Int64,SummedNumber/NumPartBases)
                if SummedNumber != 0
                    CorrelationFunctionMean[iBin,iIter] = sum(CorrelationFunctionValues[iBin, :, iIter])/SummedNumber
                else
                    CorrelationFunctionMean[iBin,iIter] = 0
                end
            end
        end
        

        return BinMidpointRadius, CorrelationFunctionMean, PartBinSummed
    end

    @time BinMidpointRadius, CorrelationValuesMean, NumPartBin = AnalyzeCorrelationSimple(NumIter, NPart, ΔBins, NumBins, NumPartBases, PartHop,SimulatedCoords, SimulatedPolAng)

end


if NaiveCorrelationAnalysis
    function NaiveCorrelation(NPart::Int64, NumIter::Int64, NumBins::Int64, MaxRadius::Float64, SimulatedCoords::Array{Float64}, SimulatedPolAng::Matrix{Float64}, NumPartBases::Int64, PartHop::Int64)
        CorrelationFunctionValues::Array{Float64} = zeros(NumBins,NumPartBases,NumIter)
        CorrelationFunctionMean::Matrix{Float64} = zeros(NumBins,NumIter)
        PartBases::Vector{Int64} = collect(1:PartHop:NPart)
        BinRadiusVector = LinRange(0, MaxRadius, NumBins+1)
        BinMidpointRadius = BinRadiusVector[2:end] .- (BinRadiusVector[2]-BinRadiusVector[1])/2
        PartWithinBinNum::Array{Int64} = zeros(Int64, NumBins, NumPartBases, NumIter)
        PartBinSummed::Matrix{Int64} = zeros(Int64, NumBins, NumIter)

        for iIter in 1:NumIter
            Threads.@threads for iParticle in 1:NumPartBases
                BasePartIdx::Int64 = PartBases[iParticle]
                Radii = sqrt.((SimulatedCoords[1,:,iIter].-SimulatedCoords[1,BasePartIdx,iIter]).^2 .+ (SimulatedCoords[2,:,iIter].-SimulatedCoords[2,BasePartIdx,iIter]).^2)
                for iBin in 1:NumBins
                    PartWithinBin::Int64 = 0
                    InnerR::Float64 = BinRadiusVector[iBin]
                    OuterR::Float64 = BinRadiusVector[iBin+1]
                    
                    IndicesInBin = findall(Radii .> InnerR .&& Radii .< OuterR)
                    for iNeighbour in IndicesInBin
                        if BasePartIdx != iNeighbour
                            PartWithinBin += 1
                            CorrelationFunctionValues[iBin,iParticle,iIter] += cos(SimulatedPolAng[BasePartIdx,iIter]-SimulatedPolAng[iNeighbour,iIter])
                        end
                    end
                    # Divide by number of particles
                    PartWithinBinNum[iBin,iParticle,iIter] = PartWithinBin
                end
            end
            if iIter%100 == 0
                println("Iteration: ",iIter)
            end
        end

        for iIter in 1:NumIter
            for iBin in 1:NumBins
                summedNumPart = sum(PartWithinBinNum[iBin,:,iIter])
                PartBinSummed[iBin,iIter] = round(Int64,summedNumPart/NumPartBases)
                if summedNumPart != 0
                    CorrelationFunctionMean[iBin,iIter] = sum(CorrelationFunctionValues[iBin,:,iIter])/summedNumPart
                end
            end
        end
        return BinMidpointRadius, CorrelationFunctionMean, PartBinSummed
    end
    @time BinMidpointRadiusNaive, CorrelationValuesMeanNaive, NumPartBinNaive = NaiveCorrelation(NPart, NumIter, NumBins, LBins, SimulatedCoords, SimulatedPolAng, NumPartBases, PartHop)
end

if saveCorrelationFunction
    function WriteCorrelationHDF5(Filename, CorrelationArray, NumBins, NumberParticles,  NumIter, BinMidpointRadius)
        fid = h5open(Filename,"w")
        
        h5File = create_group(fid, "CorrelationGroup")
        
        h5File["CorrelationArray"] = CorrelationArray
        
        h5File["NumPart"] = NumberParticles
        
        h5File["NumBins"] = NumBins
        
        h5File["BinMidpointRadius"] = BinMidpointRadius
        
        h5File["NumIter"] = NumIter
        
        close(fid)
    end

    WriteCorrelationHDF5(H5FileName, CorrelationValuesMean, NumBins, NPart,  NumIter, BinMidpointRadius)
end

if CorrelationPlot
    #WriteCorrelationHDF5("CorrelationFile.h5", CorrelationFunctionMean, NumBins, NPart,  NumIter, BinMidpointRadius)

    BinMidpointRadius = h5read("HDF5Files/CorrelationFile.h5", "CorrelationGroup/BinMidpointRadius")
    NumIter = h5read("HDF5Files/CorrelationFile.h5", "CorrelationGroup/NumIter")
    CorrelationFunctionMean = h5read("HDF5Files/CorrelationFile.h5", "CorrelationGroup/CorrelationArray")

    iPlot = Observable(1)
    Corrdata = @lift CorrelationFunctionMean[:,$iPlot]

    GLMakie.activate!(inline = false,focus_on_show=false)
    figCorr = Figure(size = (1200,900))
    axCorr = Axis(figCorr[1, 1],limits=(nothing,(-0.3, 1.05)))


    GLMakie.lines!(axCorr, BinMidpointRadius, Corrdata, label="Correlation")
    figCorr[1, 2] = GLMakie.Legend(figCorr,axCorr)
    display(figCorr)


    for i in 1:NumIter
        iPlot[] = i
        sleep(0.01)
    end
end

if CorrLengthAnalysis
    #
    # Corr <> corrcutoff
    # NumPart
    # 
    #
    #

    #
    #Outcomes
    # Corr is above, num is above -> success
    # Corr is above, num is below -> inf
    # Corr is below, num is above -> next
    # Corr is below, num is below -> Minrad
    #


    function CorrelationLengthAnalysis(BinMidpointRadius::Vector{Float64}, CorrelationValuesMean::Matrix{Float64}, NumPartBin::Matrix{Int64}, MinPartInBin::Int64,
         NumIter::Int64, MinCorrVic::Float64, CorrValCut::Float64, NumBins::Int64, InfStandin::Float64)

        CorrelationLength::Vector{Float64} = zeros(NumIter)
        
        for iIter in 1:NumIter
            if CorrelationValuesMean[1,iIter] > MinCorrVic
                for iBin in 2:NumBins
                    if iBin > floor(Int64,NumBins/2) && NumPartBin[iBin,iIter] < MinPartInBin
                        # Not enough particles are present in bin, and we are on the outer edge, corrlength is inf
                        # Inf
                        CorrelationLength[iIter] = InfStandin
                        break
                    end
                    # Either there are enough particles in bin or the bin is close to the centre.
                    if CorrelationValuesMean[iBin,iIter] < CorrValCut
                        # Success!
                        c::Float64 = ((BinMidpointRadius[iBin]-BinMidpointRadius[iBin-1])/(CorrelationValuesMean[iBin,iIter]-CorrelationValuesMean[iBin-1,iIter]))
                        CorrelationLength[iIter] = BinMidpointRadius[iBin] + (CorrValCut-CorrelationValuesMean[iBin,iIter])*c
                        break
                    end
                end
            end
        end
        return collect(1:NumIter), CorrelationLength
    end

    @time BinMidpointRadius, CorrelationValuesMean, NumPartBin = AnalyzeCorrelationSimple(NumIter, NPart, ΔBins, NumBins, NumPartBases, 
        PartHop,SimulatedCoords, SimulatedPolAng)

    InfStandin::Float64 = R*2.0::Float64
    MinPartInBin::Int64 = 15
    MinCorrVic::Float64 = 0.5
    CorrValCut::Float64 = 0.5

    @time TimeVector, CorrelationLength = CorrelationLengthAnalysis(BinMidpointRadius, CorrelationValuesMean, NumPartBin, MinPartInBin, NumIter, MinCorrVic, CorrValCut, NumBins, InfStandin)

    Plots.plot(TimeVector, CorrelationLength, size=(900,900))

end



if NormPlot

    #GLMakie.activate!(inline = false,focus_on_show=false)
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

    cmap = :matter

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
            @time SaveCoordArray, SaveRadArray,R = InitializeSystem(NPart,ρ,D,spread,finalMeanRad,growSize,NTimeSteps,Extrasteps,ExtraSlots,MoveMultiplier,A,dt,BrownianMethodFunc!,BoundaryMethodFunc!,b,SecondHeatingSteps,SecondHeatingDecrease)

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

    for i in 1:10
        if i == 1
            GLMakie.lines!(axOrder,1:NumIter,OrderArray[:,1,1],color="blue",label="N=1500")
            GLMakie.lines!(axOrder,1:NumIter,OrderArray[:,1,2],color="red",label="N=2500")
            GLMakie.lines!(axOrder,1:NumIter,OrderArray[:,1,3],color="green",label="N=5000")
            GLMakie.lines!(axOrder,1:NumIter,OrderArray[:,1,4],color="black",label="N=10000")
            GLMakie.lines!(axOrder,1:NumIter,OrderArray[:,1,5],color="black",label="N=20000")
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


# OrderArray = h5read("./OrderFileLonger.h5", "OrderGroup/OrderArray")
# NumParticlesList = h5read("./OrderFileLonger.h5", "OrderGroup/NumPart")
# NumIterList = h5read("./OrderFileLonger.h5", "OrderGroup/NumIter")

# summedOrder = zeros(NumIter,length(NumParticlesList))
# for i in 1:length(NumParticlesList)
#     for j in 1:NumIter
#         summedOrder[j,i] = mean(OrderArray[j,:,i])
#     end
# end


# GLMakie.activate!(inline = false,focus_on_show=false)
# OrderFigMean = Figure(size = (950,800))
# axOrderMean = Axis(OrderFigMean[1,1])


# GLMakie.lines!(axOrderMean,1:NumIter,summedOrder[:,1],color="blue",label="N=1500")
# GLMakie.lines!(axOrderMean,1:NumIter,summedOrder[:,2],color="red",label="N=2500")
# GLMakie.lines!(axOrderMean,1:NumIter,summedOrder[:,3],color="green",label="N=5000")
# GLMakie.lines!(axOrderMean,1:NumIter,summedOrder[:,4],color="black",label="N=10000")
# GLMakie.lines!(axOrderMean,1:NumIter,summedOrder[:,5],color="orange",label="N=20000")
# OrderFigMean[1, 2] = GLMakie.Legend(OrderFigMean,axOrderMean)

# display(OrderFigMean)

# α = Observable(1.0)

# GLMakie.activate!(inline = false,focus_on_show=false)
# OrderFigMean = Figure(size = (950,800))
# axOrderMean = Axis(OrderFigMean[1,1])



# NumIter1 = @lift (collect(1:NumIter).*0.001)./(NumParticlesList[1]^$α)
# NumIter2 = @lift (collect(1:NumIter).*0.001)./(NumParticlesList[2]^$α)
# NumIter3 = @lift (collect(1:NumIter).*0.001)./(NumParticlesList[3]^$α)
# NumIter4 = @lift (collect(1:NumIter).*0.001)./(NumParticlesList[4]^$α)
# NumIter5 = @lift (collect(1:NumIter).*0.001)./(NumParticlesList[5]^$α)

# #GLMakie.lines!(axOrderMean,NumIter1,summedOrder[:,1],color="blue",label="N=1500")
# #GLMakie.lines!(axOrderMean,NumIter2,summedOrder[:,2],color="red",label="N=2500")
# GLMakie.lines!(axOrderMean,NumIter3,summedOrder[:,3],color="green",label="N=5000")
# GLMakie.lines!(axOrderMean,NumIter4,summedOrder[:,4],color="black",label="N=10000")
# GLMakie.lines!(axOrderMean,NumIter5,summedOrder[:,5],color="orange",label="N=20000")
# OrderFigMean[1, 2] = GLMakie.Legend(OrderFigMean,axOrderMean)

# α[] = 1.0

# display(OrderFigMean)
