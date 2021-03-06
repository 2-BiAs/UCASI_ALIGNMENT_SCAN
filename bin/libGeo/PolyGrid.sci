function meshOutput = PolyGrid(pmListPolys, fSpacingNom)
    
    meshOutput = New_Mesh();
    plCombinedMasks = [];
    
    for i = 1:length(pmListPolys)
        plCombinedMasks($ + 1: $ + size(pmListPolys(i).plPoints, 1), 1:2) = pmListPolys(i).plPoints;
    end
    
    rectBoundMain = GetBoundingBox(plCombinedMasks);
    clear plCombinedMasks;
    
    gridTemp = New_Grid(rectBoundMain, fSpacingNom, "STAGGERED");
    
    meshTemp = New_Mesh(gridTemp);
    
    iTriCount = size(meshTemp.iIndices, 1);
    
    iIndexBuffer = [1:size(meshTemp.mVertices,1)];
    
    for i = 1:length(pmListPolys)
        
        iMaskIndices = [];
        
        for j = 1:size(meshTemp.mVertices, 1)
            if PointInPolygon(pmListPolys(i).plPoints, meshTemp.mVertices(j,:)) then
                iMaskIndices(1, $+1) = j; 
            end
        end
        
        select pmListPolys(i).sType
        case 'UNION'
            iIndexBuffer = union(iMaskIndices, iIndexBuffer);
        case 'DIFF'
            iIndexBuffer = setdiff(iIndexBuffer, iMaskIndices);
        case 'INT'
            iIndexBuffer = intersect(iMaskIndices, iIndexBuffer);
        else
        end
        
    end
    
    for i = 1:iTriCount
        iTriVertsInBuffer = intersect(meshTemp.iIndices(i,:), iIndexBuffer);
        if length(iTriVertsInBuffer) == 3 then
            iMaskedTriBuffer($+1, 1:3) = meshTemp.iIndices(i, 1:3); 
        end
    end
    
    //Fix Index Offsets Cause by vertex removal
    for i = 1:size(meshTemp.mVertices,1)
        for j = 1:length(iIndexBuffer)
            if iIndexBuffer(j) == i then
                iIndexLUT(i) = j;
                break;
            else
                iIndexLUT(i) = 0;
            end
        end
    end
    
    for i = 1:size(iMaskedTriBuffer,1)
        iMaskedTriBuffer(i,1) = iIndexLUT(iMaskedTriBuffer(i,1));
        iMaskedTriBuffer(i,2) = iIndexLUT(iMaskedTriBuffer(i,2));
        iMaskedTriBuffer(i,3) = iIndexLUT(iMaskedTriBuffer(i,3));
    end
    
    meshOutput.iIndices = iMaskedTriBuffer;
    
    for i = 1:length(iIndexBuffer)
        meshOutput.mVertices($+1, 1:2) = meshTemp.mVertices(iIndexBuffer(i), 1:2);
    end
    
endfunction
