using Plots
using HDF5

function DefGrid(Configuration, R, NPart, SmallRad,length,NumBase,nullpoint,Center)
    CoordArray = zeros(2,NPart)
    if Configuration=="RandomRS"
        NumPlaced = 0
        while NumPlaced < NPart
            coord = (rand(2).-0.5).*2*R.+Center
    
            if (coord[1]-Center[1])^2 + (coord[2]-Center[2])^2 < R^2
                CoordArray[:,NumPlaced+1] = coord
                NumPlaced += 1
            end
        end
    end

    if Configuration=="RandomRS_Vicinity"
        ErrorSet = 0
        success = false
        while true
            CoordArray = zeros(2,NPart)
            NumPlaced = 0
            ErrorTries = 0
            while NumPlaced < NPart
                coord = (rand(2).-0.5).*2*R.+Center
                if (coord[1]-Center[1])^2 + (coord[2]-Center[2])^2 < R^2
                    if NumPlaced > 0
                        DiffArray = CoordArray[:,1:NumPlaced].-coord
                        DiffLengths = DiffArray[1,:].^2 .+ DiffArray[2,:].^2
                        if all(DiffLengths.>SmallRad^2)
                            CoordArray[:,NumPlaced+1] = coord
                            NumPlaced += 1
                            ErrorTries = 0
                        else
                            ErrorTries += 1
                        end   
                    else
                        CoordArray[:,NumPlaced+1] = coord
                        NumPlaced += 1
                        ErrorTries = 0     
                    end
                end
                if ErrorTries > 1000
                    ErrorSet += 1
                    break
                end
                if NumPlaced == NPart
                    success = true
                    break
                end 
            end
            if success
                break
            end
            if ErrorSet > 100
                error("Too many tries, try again")
            end
        end
    end
    
    if Configuration=="HexagonalGrid"
        BaseIndex = 1
        Layer = 1
        angle = pi*60/180
        CoordArray[:,1] = [nullpoint,nullpoint]
        for i in 2:NPart
            if BaseIndex >= NumBase
                if Layer%2 == 1
                    CoordArray[:,i] = CoordArray[:,i-1] + length.*[cos(angle),sin(angle)]
                else
                    CoordArray[:,i] = CoordArray[:,i-1] + length.*[cos(2*angle),sin(2*angle)]
                end
                Layer = 3 - Layer
                BaseIndex = 1
            else
                CoordArray[:,i] = CoordArray[:,i-1] - (Layer-1.5)*2*[length,0]
                BaseIndex += 1
            end
        end
    end
    return CoordArray    
end

function FindNeighbours(NeighbouringType,CoordArray,NeighbourRad,NPart)
    if NeighbouringType == "Radius"
        IndexVec = 1:NPart
        NumNeighbours = zeros(Int64,NPart)
        for i in 1:NPart
            NumNeighbours[i] = length(findall((CoordArray[1,:][1:end .!= i].-CoordArray[1,i]).^2 + (CoordArray[2,:][1:end .!= i].-CoordArray[2,i]).^2 .< NeighbourRad^2))
        end
        if minimum(NumNeighbours) <= 1
            error("Paricle(s) present with no neighbours")
        end
        NeighbourMatrix = zeros(Int64,maximum(NumNeighbours),NPart)
        for i in 1:NPart
            if NumNeighbours[i] > 0
                NeighbourMatrix[1:NumNeighbours[i],i] = IndexVec[1:end .!= i][findall((CoordArray[1,1:end .!= i].-CoordArray[1,i]).^2 .+ (CoordArray[2,1:end .!= i].-CoordArray[2,i]).^2 .< NeighbourRad^2)]
            end
        end
    end
    return NeighbourMatrix
end

R = 6
CircR = 20
Center = [2,3]

NPart = 200
SmallRad = 0.637
Length = 0.3
NumBase = 50
NullPoint = -2*R/3
# Possibilities RandomRS, RandomRS_Vicinity, HexagonalGrid
Configuration = "RandomRS_Vicinity"

CoordArray = DefGrid(Configuration, R, NPart, SmallRad, Length, NumBase, NullPoint,Center)

NeighbourMatrix = FindNeighbours("Radius",CoordArray,SmallRad+0.5,NPart)

println(size(NeighbourMatrix))

CircleArray = zeros(2,100)
PhiCirc = LinRange(0,2*pi,100)
CircleArray[1,:] = CircR*cos.(PhiCirc)
CircleArray[2,:] = CircR*sin.(PhiCirc)


j = 24

plot(CircleArray[1,:],CircleArray[2,:],size=(800,800))
scatter!(CoordArray[1,:],CoordArray[2,:])
scatter!([CoordArray[1,j]],[CoordArray[2,j]],c="red")
scatter!(CoordArray[1,NeighbourMatrix[:,j][findall(NeighbourMatrix[:,j] .!= 0)]],CoordArray[2,NeighbourMatrix[:,j][findall(NeighbourMatrix[:,j] .!= 0)]],c="green")


