function gridOutput = New_Grid(rectBoundingBox, fNomSpacing, sType)
    
    iColumnCount = int(rectBoundMain(2,1) / fSpacingNom);
    iRowCount = int(rectBoundMain(2,2) / fSpacingNom);
    
    fGridLeft = rectBoundMain(1,1);
    fGridRight = rectBoundMain(1,1) + rectBoundMain(2,1);
    
    fGridTop = rectBoundMain(1,2);
    fGridBottom = rectBoundMain(1,2) - rectBoundMain(2,2);
    
    fX = linspace(fGridLeft, fGridRight, iColumnCount);
    fY = linspace(fGridTop, fGridBottom, iRowCount);
    
    [X, Y] = ndgrid(fX, fY);
    
    select sType
    case "SQUARE"
    //do nothing
    case "STAGGERED"
    
        fWidth = abs(fGridRight - fGridLeft);
        fActualSpacingX = fWidth / iColumnCount;
        pause
        for i = 1:size(X,1)
            for j = 1:size(X,2)
                if modulo(j, 2) == 0 then
                    X(i,j) = X(i,j) - fActualSpacingX / 4;
                else
                    X(i,j) = X(i,j) + fActualSpacingX / 4;
                end
            end
        end
    else
        pause
    end
    
    gridOutput.mX = X;
    gridOutput.mY = Y;
    gridOutput.iX_Count = iColumnCount;
    gridOutput.iY_Count = iRowCount;
    gridOutput.sType = sType;
    
endfunction
