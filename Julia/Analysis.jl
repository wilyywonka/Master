using GLMakie
using Plots
using Statistics
using HDF5

NumIter::Int64 = 1428
NPart::Int64 = 10000
MarkerSizePlot::Int64 = 2
ArrowLengthPlot::Float64 = 0.9

FilePath::String = "./HDF5Files/SaveFiles2.h5"

NormPlot::Bool = false
DiffPlot::Bool = false
QuivPolPlot::Bool = true
ColPlot::Bool = false

PolOrderAnalysis::Bool = false

CorrelationAnalysis::Bool = false
NaiveCorrelationAnalysis::Bool = false
SimpleCorrelation::Bool = true

CorrelationPlot::Bool = true

if NormPlot|DiffPlot|QuivPolPlot|ColPlot|CorrelationAnalysis|SimpleCorrelation

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

NumPartBases::Int64 = 100
PartHop::Int64 = ceil(Int64, NPart/NumPartBases)

if length(collect(1:PartHop:NPart)) != NumPartBases
    error("Number of particle-bases is wrong. Expected number is: ", NumPartBases, "; got: ",length(collect(1:PartHop:NPart)), "\n Please choose a number of particle-bases divisible by number of particles!")
end
    

ΔBins = 2.0 # 2*biggest radius, plus a little
NumBins = ceil(Int64, 2*R/ΔBins + 10)
LBins = NumBins*ΔBins

NBins = NumBins
MaxBinRadiusFraq = LBins/R
ExtraSlots::Int64 = 5
MoveMultiplier::Float64 = 1.5
PartSizeSpread::Vector{Float64} = [0.85,1.15]
SimulatedCoordsShift::Matrix{Float64} = zeros(2,NPart)


CorrelationFunctionValues::Array{Float64} = zeros(NBins,NumPartBases,NumIter)
CorrelationFunctionMean::Matrix{Float64} = zeros(NBins,NumIter)



if CorrelationAnalysis
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
            @inbounds NeighbourIndices[IndexPlacement,iCell,jCell] = iParticle

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

    function CalculateCorrelation!(iIter::Int64, NBins::Int64, BinRadiusVector::Vector{Float64}, NeighbourIndices::Array{Int64}, NeighbourCount::Matrix{Int64}, CorrelationFunctionValues::Array{Float64}, l::Float64, lDiag::Float64,
        NCells::Int64, xValsMidCells::Vector{Float64}, yValsMidCells::Vector{Float64}, SimulatedCoordsShift::Matrix{Float64}, SimulatedPolAng::Matrix{Float64}, PartHop::Int64,NPart::Int64)

        Threads.@threads for iParticle in 1:PartHop:NPart
            CurrCellX = ceil(Int64,SimulatedCoordsShift[1,iParticle]/l)
            CurrCellY = ceil(Int64,SimulatedCoordsShift[2,iParticle]/l)
            #println(CurrCellX, ", ",CurrCellY)
            
            for iBin in 1:NBins
                PartWithinBin::Int64 = 0
                InnerR::Float64 = BinRadiusVector[iBin]
                OuterR::Float64 = BinRadiusVector[iBin+1]
                OuterComp::Float64 = 0.0
                InnerComp::Float64 = 0.0

                xSpan = max(1,ceil(Int64,(xValsMidCells[CurrCellX]-(OuterR+lDiag))/l)):min(NCells, ceil(Int64,(xValsMidCells[CurrCellX]+(OuterR+lDiag))/l))
                ySpan = max(1,ceil(Int64,(yValsMidCells[CurrCellY]-(OuterR+lDiag))/l)):min(NCells, ceil(Int64,(yValsMidCells[CurrCellY]+(OuterR+lDiag))/l))

                for xCell in xSpan
                    for yCell in ySpan
                        @inbounds CellRad::Float64 = sqrt((@view(xValsMidCells[xCell])[1]-@view(xValsMidCells[CurrCellX])[1])^2 + (@view(yValsMidCells[yCell])[1]-@view(yValsMidCells[CurrCellY])[1])^2)
                        OuterComp = CellRad-lDiag
                        InnerComp = CellRad+lDiag
                        if InnerComp > InnerR && OuterComp < OuterR
                            for iPartInCell in 1:NeighbourCount[xCell,yCell]
                                @inbounds iNeighbour::Int64 = NeighbourIndices[iPartInCell,xCell,yCell]
                                if iNeighbour != iParticle
                                    PartRad::Float64 = sqrt((@view(SimulatedCoordsShift[1,iParticle])[1]-@view(SimulatedCoordsShift[1,iNeighbour])[1])^2 + (@view(SimulatedCoordsShift[2,iParticle])[1]-@view(SimulatedCoordsShift[2,iNeighbour])[1])^2)
                                    if PartRad<OuterR && PartRad>InnerR
                                        PartWithinBin += 1
                                        @inbounds @view(CorrelationFunctionValues[iBin,iParticle,iIter]) .+= cos(SimulatedPolAng[iParticle,iIter]-SimulatedPolAng[iNeighbour,iIter])
                                    end
                                end
                            end
                        end
                    end
                end
                if PartWithinBin != 0
                    # Divide by number of particles
                    CorrelationFunctionValues[iBin,iParticle,iIter] *= 1/PartWithinBin
                end
            end
        end
        return nothing
    end

    function AnalyzeCorrelation(NPart,PartHop,NBins,NumIter,MaxBinRadiusFraq,ExtraSlots::Int64,MoveMultiplier::Float64,PartSizeSpread::Vector{Float64},
        SimulatedCoordsShift::Matrix{Float64},CorrelationFunctionValues::Array{Float64},CorrelationFunctionMean::Matrix{Float64})
        
        BinRadiusVector::Vector{Float64} = collect(LinRange(0, R*MaxBinRadiusFraq,NBins+1))
        BinMidpointRadius::Vector{Float64} = BinRadiusVector[2:end] .- (BinRadiusVector[2]-BinRadiusVector[1])/2

        GridSpanDef::Vector{Float64} = [0,2*MoveMultiplier*R]
        NCells::Int64 = Int(floor((GridSpanDef[2]-GridSpanDef[1])/(2*PartSizeSpread[2])))
    
        xVec::Vector{Float64} = collect(LinRange(GridSpanDef[1],GridSpanDef[2],NCells + 1))
        yVec::Vector{Float64} = collect(LinRange(GridSpanDef[1],GridSpanDef[2],NCells + 1))

    
        l::Float64 = xVec[2]-xVec[1]
        xValsMidCells::Vector{Float64} = xVec[2:end] .- l/2
        yValsMidCells::Vector{Float64} = yVec[2:end] .- l/2

        lDiag::Float64 = sqrt(2)*l

        NeighbourCount::Matrix{Int64} = zeros(Int64,NCells,NCells)

        DivVec::Vector{Float64} = zeros(2)


        RefPart::Int64 = argmin(SimulatedCoords[1,:,1].^2 .+ SimulatedCoords[2,:,1].^2)
        SimulatedCoordsShift[1,:] .= @view(SimulatedCoords[1,:,1]) .+ (MoveMultiplier*R - SimulatedCoords[1,RefPart,1])
        SimulatedCoordsShift[2,:] .= @view(SimulatedCoords[2,:,1]) .+ (MoveMultiplier*R - SimulatedCoords[2,RefPart,1])

        NeighbourIndices, MaxPart = InitNeighbours(NeighbourCount,NCells,SimulatedCoordsShift,l,DivVec,ExtraSlots)

        # Update every 1:NumIter
        for iIter in 1:NumIter
            println("Iteration:", iIter)
            # GridSpanX::Vector{Float64} = [SimulatedCoords[1,RefPart,iIter]-MoveMultiplier*R,SimulatedCoords[1,RefPart,iIter]+MoveMultiplier*R]
            # GridSpanY::Vector{Float64} = [SimulatedCoords[2,RefPart,iIter]-MoveMultiplier*R,SimulatedCoords[2,RefPart,iIter]+MoveMultiplier*R]

            SimulatedCoordsShift[1,:] .= @view(SimulatedCoords[1,:,iIter]) .+ (MoveMultiplier*R - SimulatedCoords[1,RefPart,iIter])
            SimulatedCoordsShift[2,:] .= @view(SimulatedCoords[2,:,iIter]) .+ (MoveMultiplier*R - SimulatedCoords[2,RefPart,iIter])

            MaxPart, NeighbourIndices = UpdateNeighbours(NeighbourIndices,NeighbourCount,NCells,SimulatedCoordsShift,l,DivVec,MaxPart,ExtraSlots)
            # println(NeighbourIndices)
            # Split into bases and parallel
            CalculateCorrelation!(iIter, NBins, BinRadiusVector, NeighbourIndices, NeighbourCount, CorrelationFunctionValues, l, lDiag,
                NCells, xValsMidCells, yValsMidCells, SimulatedCoordsShift, SimulatedPolAng, PartHop, NPart)
            for iBin in 1:NBins
                CorrelationFunctionMean[iBin,iIter] = mean(CorrelationFunctionValues[iBin,:,iIter])
            end
        end
        return BinMidpointRadius, CorrelationFunctionMean
    end
    
    @profview BinMidpointRadius,CorrelationFunctionMean = AnalyzeCorrelation(NPart,PartHop,NBins,NumIter,MaxBinRadiusFraq,ExtraSlots,MoveMultiplier,PartSizeSpread,
        SimulatedCoordsShift,CorrelationFunctionValues,CorrelationFunctionMean)    

end

if SimpleCorrelation

    function AnalyzeCorrelationSimple(NumIter::Int64, NPart::Int64, ΔBins::Float64, NumBins::Int64, NumPartBases::Int64, 
        PartHop::Int64, SimulatedCoords::Array{Float64}, SimulatedPolAng::Matrix{Float64})

        CorrelationFunctionValues::Array{Float64} = zeros(NBins,NumPartBases,NumIter)
        CorrelationFunctionMean::Matrix{Float64} = zeros(NBins,NumIter)
        
        BinMidpointRadius = collect(1:NumBins)*ΔBins/2
        NumParticleInBin::Array{Int64} = zeros(Int64,NumBins,NumPartBases,NumIter)
        PartBases::Vector{Int64} = collect(1:PartHop:NPart)
        RadPartMat::Matrix{Float64} = zeros(NPart,NumPartBases)
        CurrBinMat::Matrix{Int64} = zeros(NPart-1,NumPartBases)
        ParticleNumbers::Matrix{Int64} = zeros(Int64,NPart-1, NumPartBases)
        CounterVec::Vector{Int64} = ones(Int64,NumPartBases)


        for iIter in 1:NumIter
            #Parallel
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
                    if iIter==NumIter && CurrBinMat[jNeighbour,iNumBase] == 1 && cos(SimulatedPolAng[iParticle,iIter]-SimulatedPolAng[jParticle,iIter])<0.5
                        println("ΔBins: ", ΔBins, " rad: ", RadPartMat[jNeighbour,iNumBase], " | cos: ", cos(SimulatedPolAng[iParticle,iIter]-SimulatedPolAng[jParticle,iIter]))
                        println("iNumBase ",iNumBase, " jParticle ",jParticle)
                    end
                    @inbounds @view(CorrelationFunctionValues[CurrBinMat[jNeighbour,iNumBase],iNumBase,iIter]) .+= cos(SimulatedPolAng[iParticle,iIter]-SimulatedPolAng[jParticle,iIter])
                end
            end
            println("Iteration: ",iIter)
        end

        for iIter in 1:NumIter
            for iBin in 1:NumBins
                SummedNumber = sum(NumParticleInBin[iBin, :, iIter])

                if SummedNumber != 0
                    CorrelationFunctionMean[iBin,iIter] = sum(CorrelationFunctionValues[iBin, :, iIter])/SummedNumber
                else
                    CorrelationFunctionMean[iBin,iIter] = 0
                end
            end
        end
        

        return BinMidpointRadius, CorrelationFunctionMean
    end
    
#    BinMidpointRadius, CorrelationFunctionMean1 = AnalyzeCorrelationSimple(10, NPart, ΔBins, NumBins, NumPartBases, 
 #       PartHop, CorrelationFunctionValues, CorrelationFunctionMean,SimulatedCoords, SimulatedPolAng)

    #@profview BinMidpointRadius, CorrelationFunctionMean2 = AnalyzeCorrelationSimple2(10, NPart, ΔBins, NumBins, NumPartBases, 
     #   PartHop, CorrelationFunctionValues, CorrelationFunctionMean,SimulatedCoords, SimulatedPolAng)


    @time BinMidpointRadius, CorrelationFunctionMean = AnalyzeCorrelationSimple(NumIter, NPart, ΔBins, NumBins, NumPartBases, 
        PartHop,SimulatedCoords, SimulatedPolAng)

end
# iIter = 1428
# iNumBase = 29
# jParticle = 
# aPlot = Plots.scatter(SimulatedCoords[1,:,iIter],SimulatedCoords[2,:,iIter],color="blue")
# Plots.scatter!([SimulatedCoords[1,jParticle,iIter]],[SimulatedCoords[2,jParticle,iIter]],color="red")
# Plots.scatter!([SimulatedCoords[1,iNumBase,iIter]],[SimulatedCoords[2,iNumBase,iIter]],color="green")
# Plots.quiver!(SimulatedCoords[1,:,iIter],SimulatedCoords[2,:,iIter],quiver=(cos.(SimulatedPolAng[:,iIter]), sin.(SimulatedPolAng[:,iIter])))
# Plots.xlims!([SimulatedCoords[1,jParticle,iIter]-3,SimulatedCoords[1,jParticle,iIter]+3])
# Plots.ylims!([SimulatedCoords[2,jParticle,iIter]-3,SimulatedCoords[2,jParticle,iIter]+3])
# display(aPlot)

a = GLMakie.lines(BinMidpointRadius,CorrelationFunctionMean[:,1428])

# if NaiveCorrelationAnalysis
#     # NBins = 20
#     # MaxBinRadiusFraq = 0.5
#     # NumPartBases::Int64 = 1000

#     # CorrelationFunctionValues::Array{Float64} = zeros(NBins,NumPartBases,NumIter)
#     # CorrelationFunctionMean::Matrix{Float64} = zeros(NBins,NumIter)

#     BinRadiusVector = LinRange(0, R*MaxBinRadiusFraq,NBins+1)
#     BinMidpointRadius = BinRadiusVector[2:end] .- (BinRadiusVector[2]-BinRadiusVector[1])/2

#     for iIter in 1:NumIter
#         for iParticle in 1:NumPartBases
#             Radii = sqrt.((SimulatedCoords[1,:,iIter].-SimulatedCoords[1,iParticle,iIter]).^2 .+ (SimulatedCoords[2,:,iIter].-SimulatedCoords[2,iParticle,iIter]).^2)
#             for iBin in 1:NBins
#                 PartWithinBin::Int64 = 0
#                 InnerR::Float64 = BinRadiusVector[iBin]
#                 OuterR::Float64 = BinRadiusVector[iBin+1]
                
#                 IndicesInBin = findall(Radii .> InnerR .&& Radii .< OuterR)
#                 for iNeighbour in IndicesInBin
#                     if iParticle != iNeighbour
#                         PartWithinBin += 1
#                         CorrelationFunctionValues[iBin,iParticle,iIter] += cos(SimulatedPolAng[iParticle,iIter]-SimulatedPolAng[iNeighbour,iIter])
#                     end
#                 end
#                 if PartWithinBin != 0
#                     # Divide by number of particles
#                     CorrelationFunctionValues[iBin,iParticle,iIter] *= 1/PartWithinBin
#                 end
#             end
#         end
#         for iBin in 1:NBins
#             CorrelationFunctionMean[iBin,iIter] = mean(CorrelationFunctionValues[iBin,:,iIter])
#         end
#         println("Iteration: ",iIter)
#     end
#     Plots.plot(BinMidpointRadius,CorrelationFunctionMean[:,1])
# end


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

WriteCorrelationHDF5("HDF5Files/CorrelationFile.h5", CorrelationFunctionMean, NumBins, NPart,  NumIter, BinMidpointRadius)

if CorrelationPlot
    #WriteCorrelationHDF5("CorrelationFile.h5", CorrelationFunctionMean, NBins, NPart,  NumIter, BinMidpointRadius)

    BinMidpointRadius = h5read("HDF5Files/CorrelationFile.h5", "CorrelationGroup/BinMidpointRadius")
    NumIter = h5read("HDF5Files/CorrelationFile.h5", "CorrelationGroup/NumIter")
    CorrelationFunctionMean = h5read("HDF5Files/CorrelationFile.h5", "CorrelationGroup/CorrelationArray")

    iPlot = Observable(1)
    Corrdata = @lift CorrelationFunctionMean[:,$iPlot]

    GLMakie.activate!(inline = false,focus_on_show=false)
    figCorr = Figure(size = (950,800))
    axCorr = Axis(figCorr[1, 1],limits=(nothing,(-0.05, 1.05)))


    GLMakie.lines!(axCorr, BinMidpointRadius, Corrdata, label="Correlation")
    figCorr[1, 2] = GLMakie.Legend(figCorr,axCorr)
    display(figCorr)


    for i in 1:NumIter
        iPlot[] = i
        sleep(0.01)
    end
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


if true
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

