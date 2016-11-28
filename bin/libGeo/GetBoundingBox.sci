function rectOutput = GetBoundingBox(plPoints)
    
    xMin = min(plPoints(:,1));
    xMax = max(plPoints(:,1));
    yMin = min(plPoints(:,2));
    yMax = max(plPoints(:,2));
    
    rectOutput = [xMin, yMax; xMax - xMin, yMax - yMin]; //[Left Top; Width Hieght]
endfunction
