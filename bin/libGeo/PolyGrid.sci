function meshOutput = PolyGrid(pmListPolys, fSpacingNom)
    
    meshOutput = New_Mesh();
    plCombinedMasks = [];
    
    for i = 1:length(pmListPolys)
        plCombinedMasks($ + 1: $ + size(pmListPolys(i).plPoints, 1), 1:2) = pmListPolys(i).plPoints;
    end
    
    rectBoundMain = GetBoundingBox(plCombinedMasks);
    clear plCombinedMasks;
    
    iColumnCount = int(rectBoundMain(2,1) / fSpacingNom);
    iRowCount = int(rectBoundMain(2,2) / fSpacingNom);
    
    fGridLeft = rectBoundMain(1,1);
    fGridRight = rectBoundMain(1,1) + rectBoundMain(2,1);
    
    fGridTop = rectBoundMain(1,2);
    fGridBottom = rectBoundMain(1,2) - rectBoundMain(2,2);
    
    fX = linspace(fGridLeft, fGridRight, iColumnCount);
    fY = linspace(fGridTop, fGridBottom, iRowCount);
    
    //mprintf("\nGenerate starting grid, and call GridToList.\n");
    //pause
    [X, Y] = ndgrid(fX, fY);
    mTempGrid = GridToList(X, Y);
    
    iTriBuffer = [];
    
    for i = 0:iRowCount - 1
        for j = 1:iColumnCount
            iTriBuffer($+1, 1) = i * iColumnCount + j;
            iTriBuffer($, 2) = i * iColumnCount + j + 1;
            iTriBuffer($, 3) = (i + 1) * iColumnCount + j;
            
            iTriBuffer($+1, 1) = (i + 1) * iColumnCount + j + 1;
            iTriBuffer($, 2) = (i + 1) * iColumnCount + j;
            iTriBuffer($, 3) = i * iColumnCount + j + 1;
        end
    end
    
    iTriCount = size(iTriBuffer, 1);
    
    iIndexBuffer = [1:size(mTempGrid,1)];
    
    for i = 1:length(pmListPolys)
        
        iMaskIndices = [];
        
        for j = 1:size(mTempGrid, 1)
            if PointInPolygon(pmListPolys(i).plPoints, mTempGrid(j,:)) then
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
        iTriVertsInBuffer = intersect(iTriBuffer(i,:), iIndexBuffer);
        if length(iTriVertsInBuffer) == 3 then
            iMaskedTriBuffer($+1, 1:3) = iTriBuffer(i, 1:3); 
        end
    end
    
    //Fix Index Offsets Cause by vertex removal
    for i = 1:size(mTempGrid,1)
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
        meshOutput.mVertices($+1, 1:2) = mTempGrid(iIndexBuffer(i), 1:2);
    end
    
endfunction
