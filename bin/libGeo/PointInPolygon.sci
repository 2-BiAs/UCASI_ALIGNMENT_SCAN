function [bIsInside] = PointInPolygon(pgInput, vPoint)
//*****************************************************************************
// function: bIsInside_poly - this function returns a Yes (1), No(0)
// depending on if a point is bIsInside a  polygon or not
//
// Inputs:
//    xpol, ypol are the coordinates of a simple closed polygon 
//    vPoint(1), vPoint(1) are the coordinates of the point
//
// Outputs: 1 (bIsInside polygon), 0 (not bIsInside polygon)
//
//*****************************************************************************

//********************************************
// show the polygon
do_the_plot = 0

if(do_the_plot)
    scf(9999);
    // it is assumed that this is a closed polygon
    plot2d(xtemp,ytemp)
end
//********************************************

iN = size(pgInput, 2)


bIsInside = %F
j = iN; // j is the previous vertice 
i = 1
while  i <= iN

    if ((((pgInput(i, 2) <= vPoint(1)) & (vPoint(1) < pgInput(j, 2)))|((pgInput(j, 2) <= vPoint(1)) & (vPoint(1) < pgInput(i, 2)))) & ..
       (vPoint(1) < ((pgInput(j, 1) - pgInput(i, 1))/(pgInput(j, 2) - pgInput(i, 2))) * (vPoint(1) - pgInput(i, 2)) + pgInput(i, 1)))

          bIsInside = ~bIsInside;
    end
    i = i + 1;
    j = i - 1;
end

endfunction
