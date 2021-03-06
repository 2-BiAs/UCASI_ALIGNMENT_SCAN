function mOutput = GetSags(mPoints, spParams)
    //Generates an array of sag values given an array of [u, v] points
    //And surface parameter structure defining general polynomial freeform surface
    
    format('e', 16);

    //Build Polynomial Terms
    sPolynomialTerms = '';
    iNM = size(spParams.mPolyCoef);
    iN = iNM(1);
    iM = iNM(2); 
    for i = 1:iN
        for j = 1:iM
            if spParams.mPolyCoef(i, j) ~= 0 then
                sPolynomialTerms = sPolynomialTerms + '+ ' + string(spParams.mPolyCoef(i,j)) +...
                ' * (x - fXs * ones(x)).^(' + string(i) + ' - 1) .* (y - fYs * ones(y)).^(' + string(j) + ' - 1)';
            end
        end
    end
    
    //Set Surface Parameters
    fXs = 0;
    fYs = 0;
    fConicConstant = spParams.fCC;
    
    //Build Conic Equation 
    sConic = 'z = (x.^2 - 2 * x * fXs + fXs^2 * ones(x) + y.^2 - 2 * y * fYs + fYs^2 * ones(y)) .*...
    (R * (ones(x) + sqrt(ones(x) - (1 + fConicConstant) * (x.^2 - 2 * x * fXs + fXs^2 * ones(x) + y.^2 - 2 * y * fYs + fYs^2 * ones(y)) / R^2))) .^ (-1) ';
    
    //Define Surface Funtion
    deff('z = Z(x, y, R)', sConic + sPolynomialTerms);
    
    x = mPoints(:,1);
    y = mPoints(:,2);
    
    format('v', 10);
    
    mOutput = Z(x, y, spParams.fROC);
    
endfunction
