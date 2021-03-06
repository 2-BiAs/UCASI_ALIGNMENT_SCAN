function mOutput = GridToList(mX, mY)
    //Takes a Grid Generated by ndgrid, or such, and converts to a list
    
    [iN, iM] = size(mX);

    for i = 1:1:iM
        for j = 1:1:iN
            mOutput((i - 1) * iN + j, 1:2) = [mX(j, i) mY(j, i)];
        end
    end
      
endfunction
