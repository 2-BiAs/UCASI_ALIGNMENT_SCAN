function meshOutput = New_Mesh(varargin)
    
    //Generate New Mesh structure
    
    meshOutput.mVertices = [];
    meshOutput.iIndices = [];
    
    select argn(2)
    case 1;

        iX_Count = varargin(1).iX_Count;
        iY_Count = varargin(1).iY_Count;

        for i = 0:iY_Count - 1
            for j = 1:iX_Count
                iTriBuffer($+1, 1) = i * iX_Count + j;
                iTriBuffer($, 2) = i * iX_Count + j + 1;
                iTriBuffer($, 3) = (i + 1) * iX_Count + j;

                iTriBuffer($+1, 1) = (i + 1) * iX_Count + j + 1;
                iTriBuffer($, 2) = (i + 1) * iX_Count + j;
                iTriBuffer($, 3) = i * iX_Count + j + 1;
            end
        end
        
        mVertices = GridToList(varargin(1).mX, varargin(1).mY)

        meshOutput.mVertices = mVertices;
        meshOutput.iIndices = iTriBuffer;

    case 2;
        meshOutput.mVertices = varargin(1);
        meshOutput.iIndices = varargin(2);
    case 3;
    else
    end
endfunction
