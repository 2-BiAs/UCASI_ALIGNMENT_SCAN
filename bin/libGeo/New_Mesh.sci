function meshOutput = New_Mesh(varargin)
    
    //Generate New Mesh structure
    
    meshOutput.mVertices = [];
    meshOutput.iIndices = [];
    
    select argn(2)
    case 2;
        meshOutput.mVertices = varargin(1);
        meshOutput.iIndices = varargin(2);
    else
    end
endfunction
