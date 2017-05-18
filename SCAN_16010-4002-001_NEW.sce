clear
clc
clf(1); clf(2); clf(3); clf(4)
format('v', 10)
csvDefault("eol", "windows")
csvDefault("blank", "on")

[ok, sPartNumber] = getvalue("Input Part Number", "PN", list("str", 1), ["16010-4002-001"]);
if ~ok then
    abort
end

[ok, sSerialNumber] = getvalue("Input Serial Number", "SN", list("str", 1), ["default"]);
if ~ok then
    abort
end

directory = uigetdir(pwd(), "Select Working Directoy");
status = chdir(directory)
realpath = cd(directory)
if ~status then
  disp('failed to change directory')
  abort
end

if ~isdir(directory) then
    [status, err] = mkdir(sSerialNumber, 'imgs')
    if status ~= 1 then
       disp(err)
       abort
    end
end

//exec('WriteMatrixCSV.sci');
//exec('PointMatrixToList.sci');

//Define Surface Parameters
fConicConstant = -1.374;
fRadiusOfCurvature = 38.448;    //(mm)
iConcavity = +1;    //-1 = concave //+1 = convex //no 0!
fOffset = 0.0; // mm

fPolynomialCoefficients = [0];
//...
//    [           0,             0,  2.8265179e-6, -4.8964451e-7;...
//                0,  6.5191331e-5, -1.8966279e-7,             0;...
//    -5.9920375e-5, -7.0557011e-7,             0,             0;...
//    -2.6125317e-7,             0,             0,             0];
//

fRadiusOfCurvature = abs(fRadiusOfCurvature) * -iConcavity;
fPolynomialCoefficients = fPolynomialCoefficients * -iConcavity;

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//Define Scan Parameters
fOD = 70;   //(mm)
fID = 28;

iN_R = 24; //20;  //Number of radial steps
iN_T = 16; //20;  //Number of angular steps

fXs = 0;
fYs = 0;

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////



//Build Polynomial Terms
sPolynomialTerms = '';
iNM = size(fPolynomialCoefficients);
iN = iNM(1);
iM = iNM(2); 
for i = 1:iN
    for j = 1:iM
        if fPolynomialCoefficients(i, j) ~= 0 then
            sPolynomialTerms = sPolynomialTerms + '+ ' + string(fPolynomialCoefficients(i,j)) +...
            ' * (x - fXs * ones(x)).^(' + string(i) + ' - 1) .* (y - fYs * ones(y)).^(' + string(j) + ' - 1)';
        end
    end
end

//Build Polynomial Normal Terms
sPnX = '';
sPnY = '';
iNM = size(fPolynomialCoefficients);
iN = iNM(1);
iM = iNM(2); 
for i = 1:iN
    for j = 1:iM
        if fPolynomialCoefficients(i, j) ~= 0 then
            sPnX = sPnX + '- ' + string(fPolynomialCoefficients(i,j)) +...
            ' * ((x - fXs * ones(x)).^' + string(i-1) +' * ' + string(i-1) + ') .* (y - fYs * ones(y)).^' + string(j-1) + ' .* (x - fXs * ones(x)).^(-1)';            
            
            sPnY = sPnY + '- ' + string(fPolynomialCoefficients(i,j)) +...
            ' * ((x - fXs * ones(x)).^' + string(i-1) + ' * ' + string(j-1) + ') .* (y - fYs * ones(y)).^' + string(j-1) + ' .* (y - fYs * ones(y)).^(-1)';
        end
    end
end

//Build Conic Equation 
sConic = 'z = (x.^2 - 2 * x * fXs + fXs^2 * ones(x) + y.^2 - 2 * y * fYs + fYs^2 * ones(y)) .*...
    (R * (ones(x) + sqrt(ones(x) - (1 + fConicConstant) * (x.^2 - 2 * x * fXs + fXs^2 * ones(x) + y.^2 - 2 * y * fYs + fYs^2 * ones(y)) / R^2))) .^ (-1) ';

//Build Conic Normals
sCnX = 'nx = -(2 * x - 2 * fXs * ones(x)) .* (R * (ones(x) + sqrt(ones(1)-(x.^2 - 2 * x * fXs + fXs^2 * ones(x) + y.^2 - 2 * y * fYs + fYs^2 * ones(y)) * (1 + fConicConstant) / R^2))).^(-1)-...
(1 / 2) * (x.^2 - 2 * x * fXs + fXs^2 * ones(x) + y.^2 - 2 * y * fYs + fYs^2 * ones(y)) .* (1 + fConicConstant) .* (2 * x - 2 * fXs * ones(x)) .*...
((R^3 * (ones(x) + sqrt(ones(x) - (x.^2 - 2 * x * fXs + fXs^2 * ones(x) + y.^2 - 2 * y * fYs + fYs^2 * ones(y)) * (1 + fConicConstant) / R^2)).^2)...
.* sqrt(ones(x) - (x.^2 - 2 * x * fXs + fXs^2 * ones(x) + y.^2 - 2 * y * fYs + fYs^2 * ones(y)) * (1 + fConicConstant) / R^2)).^(-1) ';

sCnY = 'ny = -(2 * y - 2 * fYs * ones(y)) .* (R * (ones(x) + sqrt(ones(1)-(x.^2 - 2 * x *fXs+fXs^2 * ones(x) + y.^2 - 2 * y * fYs + fYs^2 * ones(y)) * (1 + fConicConstant) / R^2))).^(-1)-...
(1 / 2) * (x.^2 - 2 * x *fXs + fXs^2 * ones(x) + y.^2 - 2 * y *fYs+fYs^2 * ones(y)) * (1 + fConicConstant) .* (2 * y - 2 *fYs* ones(y)) .*...
((R^3 * (ones(x) + sqrt(ones(x) - (x.^2 - 2 * x *fXs+fXs^2 * ones(x) + y.^2 - 2 * y *fYs+fYs^2 * ones(y)) * (1 + fConicConstant) / R^2)).^2)...
.* sqrt(ones(x) - (x.^2 - 2 * x *fXs+fXs^2 * ones(x) + y.^2 - 2 * y *fYs+fYs^2 * ones(y)) * (1 + fConicConstant) / R^2)).^(-1) ';

//Define Surface Funtions
deff('z = Z(x, y, R)', sConic + sPolynomialTerms);
deff('x = X(r, t)', 'x = r'' * cos(t)');
deff('y = Y(r, t)', 'y = r'' * sin(t)');
//deff('x = X(u, v)', 'x = (u * cos(fTheta))'' * ones(v) - ones(u'') * (v * sin(fTheta))');
//deff('y = Y(u, v)', 'y = (u * sin(fTheta))'' * ones(v) + ones(u'') * (v * cos(fTheta))');

//Define Normal Fuctions
deff('nx = nX(x, y, R)', sCnX + sPnX);
deff('ny = nY(x, y, R)', sCnY + sPnY);


//Define Domain in the Polar Parameters (r, t)
r = linspace(fID / 2, fOD / 2, iN_R);
t = linspace(0, 2 * %pi, iN_T);

x = X(r, t);
y = Y(r, t);
z = Z(x, y, fRadiusOfCurvature);

NX = nX(x, y, fRadiusOfCurvature);
NY = nY(x, y, fRadiusOfCurvature);
NZ = ones(NX);

NX = NX;   //Flip normals? or not
NY = NY;
NZ = NZ;


fig_the = scf(1);
[xf, yf, zf] = nf3d(x', y', z'); 
xset('colormap', jetcolormap(128));
plot3d1(xf, yf, zf, theta = 300, alpha = 60, leg = "@@", flag = [-1 6 4]);
axes_the = gca();
axes_the.axes_reverse = ["off", "off", "off"];
//facets_the = axes_the.children(1);
//facets_the.hiddencolor = -1;
//colorbar(min(z), max(z));


INV_NMAG = (NX.^2 + NY.^2 + NZ.^2).^(-1/2);

NX = NX .* INV_NMAG;
NY = NY .* INV_NMAG;
NZ = NZ .* INV_NMAG;

[a, b] = size(x);
for i = 1:a
    x_list((i - 1) * b + 1 : i * b) = x(i, :)  
end

[a, b] = size(y);
for i = 1:a
    y_list((i - 1) * b + 1 : i * b) = y(i, :)  
end

[a, b] = size(z);
for i = 1:a
    z_list((i - 1) * b + 1 : i * b) = z(i, :)  
end

[a, b] = size(NX);
for i = 1:a
    NX_list((i - 1) * b + 1 : i * b) = NX(i, :)  
end

[a, b] = size(NY);
for i = 1:a
    NY_list((i - 1) * b + 1 : i * b) = NY(i, :)  
end

[a, b] = size(NZ);
for i = 1:a
    NZ_list((i - 1) * b + 1 : i * b) = NZ(i, :)  
end

xarrows([x_list; x_list + NX_list*10], [y_list; y_list + NY_list*10], [z_list; z_list + NZ_list*10], 25)

//xs2pdf(gcf(), pwd() + '\' + string(sSerialNumber) + '\imgs\' + string(sPartNumber) + '_'...
//    + string(sSerialNumber) + '_SURFACE.pdf');

/////////////////////////////////Apply normal offset to target points

x_list = x_list + fOffset * NX_list;
y_list = y_list + fOffset * NY_list;
z_list = z_list + fOffset * NZ_list;

///////////////////////////////////////////////////////////////// SAVE TARGET POINTS
sFileToSave = 'TARGET_POINTS.csv';
sTempFile = TMPDIR + "\" + sFileToSave;
fFile = mopen(sTempFile, 'wt');

for i=1:length(x_list) - 1
    mfprintf(fFile, "%.6f, %.6f, %.6f, %.6f, %.6f, %.6f\n", x_list(i), y_list(i), z_list(i), NX_list(i), NY_list(i), NZ_list(i));
end

mfprintf(fFile, "%.6f, %.6f, %.6f, %.6f, %.6f, %.6f", x_list($), y_list($), z_list($), NX_list($), NY_list($), NZ_list($));

mclose(fFile);

dos('move ' + sTempFile + ' ' + pwd() + '\' + sFileToSave);
/////////////////////////////////////////////////////////////////



M_Actual = csvRead('ACTUAL_POINTS.csv');


[a b] = size(M_Actual)
for i=1:iN_R
    x_actual(i, 1:iN_T) = M_Actual((i - 1) * iN_T + 1:i * iN_T, 1)';
    y_actual(i, 1:iN_T) = M_Actual((i - 1) * iN_T + 1:i * iN_T, 2)';
    z_actual(i, 1:iN_T) = M_Actual((i - 1) * iN_T + 1:i * iN_T, 3)';
end

[a, b] = size(x_actual);
for i = 1:a
    x_actual_list((i - 1) * b + 1 : i * b) = x_actual(i, :)  
end

[a, b] = size(y_actual);
for i = 1:a
    y_actual_list((i - 1) * b + 1 : i * b) = y_actual(i, :)  
end

[a, b] = size(z_actual);
for i = 1:a
    z_actual_list((i - 1) * b + 1 : i * b) = z_actual(i, :)  
end

z_error = z_actual - Z(x_actual, y_actual, fRadiusOfCurvature);

fig_error2 = scf(3);
xset('colormap', jetcolormap(128));
[xf, yf, zf] = nf3d(x_actual', y_actual', z_error' * 1000); 
plot3d1(xf, yf, zf, theta = 300, alpha = 60, leg = "@@", flag = [-1 6 4]);
colorbar(min(z_error * 1000), max(z_error * 1000),,fmt="%.3f");

a = gca();
a.view = "2d";
a.axes_reverse = ["off", "off", "off"];

format('v', 10)
a.x_label.text = "$\Large X \tt(mm)$";
a.y_label.text = "$\Large Y \tt(mm)$";
a.title.text = "$\text{\begin{gather}\Huge{Surface \ Error}$" + ...
    "$\\ \LARGE{" + sPartNumber + " \ SN" + sSerialNumber + "}$" + ...
    "$\\ \small{R_{nom} = " + string(fRadiusOfCurvature) + "mm}$" +...
     "\end{gather}}$";

mkdir('/' + sSerialNumber);
mkdir('/' + sSerialNumber + '/imgs');
xs2pdf(gcf(), pwd() + '\' + string(sSerialNumber) + '\imgs\' + string(sPartNumber) + '_'...
    + string(sSerialNumber) + '_FIGURE_ERROR_REAL.pdf');


function zt = ZT(x, y, parameters)
//    zt = Z(x - parameters(1) + parameters(4) * sign(x),...
//           y - parameters(2) + parameters(5) * sign(y),...
//           fRadiusOfCurvature) + parameters(3);
    zt = Z(x - parameters(1) + parameters(4) * cos(atan(y, x)),...
           y - parameters(2) + parameters(4) * sin(atan(y, x)),...
           fRadiusOfCurvature) + parameters(3);
endfunction

function ze = ZE(parameters, x_exp, y_exp, z_exp)
    ze = ZT(x_exp, y_exp, parameters) - z_exp; 
endfunction

param_noms = [0, 0, 0, 0, 0];
//param_binf = [-%inf, -%inf, -%inf, -.020, -%inf];
//param_bsup = [%inf, %inf, %inf, -.019, %inf];

[f, param_opt] =  leastsq(list(ZE, x_actual_list', y_actual_list', z_actual_list'), param_noms);
//[f, param_opt] =  leastsq(list(ZE, x_actual_list', y_actual_list', z_actual_list'),'b', param_binf, param_bsup, param_noms);


//z_error_ls_fit_all = z_actual - (Z(x_actual - param_opt(1), y_actual - param_opt(2), param_opt(4)) + param_opt(3));
//z_error_ls_fit_all = - z_error_ls_fit_all; //This makes the sign correct so that a positive error is in the outward surface normal dirrection

z_error_ls = z_actual - (Z(x_actual - param_opt(1), y_actual - param_opt(2), fRadiusOfCurvature)); 
z_error_ls_xyz = z_actual - (Z(x_actual - param_opt(1), y_actual - param_opt(2), fRadiusOfCurvature) + param_opt(3));
z_error_ls_xyz_tool = z_actual - (Z(x_actual - param_opt(1) + param_opt(4) * cos(atan(y_actual, x_actual)), y_actual - param_opt(2) + param_opt(4) * sin(atan(y_actual, x_actual)), fRadiusOfCurvature) + param_opt(3));
//z_error_ls = -z_error_ls; //This makes the sign correct so that a positive error is in the outward surface normal dirrection

fig_error2 = scf(4);
xset('colormap', jetcolormap(128));
[xf, yf, zf] = nf3d(x_actual', y_actual', z_error_ls' * 1000); 
plot3d1(xf, yf, zf, theta = 300, alpha = 60, leg = "@@", flag = [-1 6 4]);
colorbar(min(z_error_ls * 1000), max(z_error_ls * 1000),,fmt="%.3f");

e = gce();
e.parent.title.text = "$\Large Z_{err} \tt(\mu m)$";

a = gca();
a.view = "2d";
a.axes_reverse = ["off", "off", "off"];

format('v', 10)
a.x_label.text = "$\Large X \tt(mm)$";
a.y_label.text = "$\Large Y \tt(mm)$";
a.title.text = "$\text{\begin{gather}\Huge{Surface \ Error}$" + ...
    "$\\ \LARGE{" + sPartNumber + " \ SN" + sSerialNumber + "}$" + ...
    "$\\ \normalsize{R_{nom} = " + string(fRadiusOfCurvature) + "mm}$" +...
    "$\\ \normalsize{X_{offset} = " + string(round(param_opt(1)*1000)/1000) + "mm}$" +...
    "$\\ \normalsize{Y_{offset} = " + string(round(param_opt(2)*1000)/1000) + "mm}$" +...
     "\end{gather}}$";

//facets_a = a.children(1);
//facets_a.hiddencolor = -1;

mkdir('/' + sSerialNumber);
mkdir('/' + sSerialNumber + '/imgs');
xs2pdf(gcf(), pwd() + '\' + string(sSerialNumber) + '\imgs\' + string(sPartNumber) + '_'...
    + string(sSerialNumber) + '_FIGURE_ERROR_FIT_XY.pdf');

//fig_error2 = scf(5);
//xset('colormap', jetcolormap(128));
//[xf, yf, zf] = nf3d(x_actual', y_actual', z_error_ls_fit_all' * 1000); 
//plot3d1(xf, yf, zf, theta = 300, alpha = 60, leg = "@@", flag = [-1 6 4]);
//colorbar(min(z_error_ls_fit_all * 1000), max(z_error_ls_fit_all * 1000),,fmt="%.3f");
//
//e = gce();
//e.parent.title.text = "$\Large Z_{err} \tt(\mu m)$";
//
//a = gca();
//a.view = "2d";
//a.axes_reverse = ["off", "off", "off"];
//
//format('v', 10)
//a.x_label.text = "$\Large X \tt(mm)$";
//a.y_label.text = "$\Large Y \tt(mm)$";
//a.title.text = "$\text{\begin{gather}\Huge{Surface \ Error}$" + ...
//    "$\\ \LARGE{" + sPartNumber + " \ SN" + sSerialNumber + "}$" + ...
//    "$\\ \small{R_{bf} = " + string(round(param_opt(4)*1000)/1000) + "mm}$" +...
//    "$\\ \small{X_{offset} = " + string(round(param_opt(1)*1000)/1000) + "mm}$" +...
//    "$\\ \small{Y_{offset} = " + string(round(param_opt(2)*1000)/1000) + "mm}$" +...
//    "$\\ \small{Z_{offset} = " + string(round(param_opt(3)*1000)/1000) + "mm}$" +...
//     "\end{gather}}$";
//
////facets_a = a.children(1);
////facets_a.hiddencolor = -1;
//
//mkdir('/' + sSerialNumber);
//mkdir('/' + sSerialNumber + '/imgs');
//xs2pdf(gcf(), pwd() + '\' + string(sSerialNumber) + '\imgs\' + string(sPartNumber) + '_'...
//    + string(sSerialNumber) + '_FIGURE_ERROR_FIT_XYZR.pdf');
    
    ////////////////////////////////////////////////////////
    
    fig_error2 = scf(7);
xset('colormap', jetcolormap(128));
[xf, yf, zf] = nf3d(x_actual', y_actual', z_error_ls_xyz_tool' * 1000); 
plot3d1(xf, yf, zf, theta = 300, alpha = 60, leg = "@@", flag = [-1 6 4]);
colorbar(min(z_error_ls_xyz_tool * 1000), max(z_error_ls_xyz_tool * 1000),,fmt="%.3f");

e = gce();
e.parent.title.text = "$\Large Z_{err} \tt(\mu m)$";

a = gca();
a.view = "2d";
a.axes_reverse = ["off", "off", "off"];

format('v', 10)
a.x_label.text = "$\Large X \tt(mm)$";
a.y_label.text = "$\Large Y \tt(mm)$";
a.title.text = "$\text{\begin{gather}\Huge{Surface \ Error}$" + ...
    "$\\ \LARGE{" + sPartNumber + " \ SN" + sSerialNumber + "}$" + ...
    "$\\ \normalsize{X_{shift} = " + string(round(param_opt(4)*1000)/1000) + "mm}$" +...
    "$\\ \normalsize{Y_{shift} = " + string(round(param_opt(5)*1000)/1000) + "mm}$" +...
    "$\\ \normalsize{X_{offset} = " + string(round(param_opt(1)*1000)/1000) + "mm}$" +...
    "$\\ \normalsize{Y_{offset} = " + string(round(param_opt(2)*1000)/1000) + "mm}$" +...
    "$\\ \normalsize{Z_{offset} = " + string(round(param_opt(3)*1000)/1000) + "mm}$" +...
     "\end{gather}}$";

//facets_a = a.children(1);
//facets_a.hiddencolor = -1;

mkdir('/' + sSerialNumber);
mkdir('/' + sSerialNumber + '/imgs');
xs2pdf(gcf(), pwd() + '\' + string(sSerialNumber) + '\imgs\' + string(sPartNumber) + '_'...
    + string(sSerialNumber) + '_FIGURE_ERROR_FIT_XYZT.pdf');
    
    //////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//
//
//[M N] = size(x)
//
//for i = 1 : M - 1
//    for j = 1 : N - 1
//        grad_zu(i, j) = ((z_error_ls_fit_all(i, j) - z_error_ls_fit_all(i + 1, j)) / sqrt((x_actual(i, j) - x_actual(i + 1, j)) ^ 2 + (y_actual(i, j) - y_actual(i + 1, j)) ^ 2) + (z_error_ls_fit_all(i, j + 1) - z_error_ls_fit_all(i + 1, j + 1)) / sqrt((x_actual(i, j + 1) - x_actual(i + 1, j + 1)) ^ 2 + (y_actual(i, j + 1) - y_actual(i + 1, j + 1)) ^ 2)) / 2;
//        
//        grad_zv(i, j) = ((z_error_ls_fit_all(i, j) - z_error_ls_fit_all(i, j + 1)) / sqrt((x_actual(i, j) - x_actual(i, j + 1)) ^ 2 + (y_actual(i, j) - y_actual(i, j + 1)) ^ 2) + (z_error_ls_fit_all(i + 1, j) - z_error_ls_fit_all(i + 1, j + 1)) / sqrt((x_actual(i + 1, j) - x_actual(i + 1, j + 1)) ^ 2 + (y_actual(i + 1, j) - y_actual(i + 1, j + 1)) ^ 2)) / 2;
//         
//        
//        x_grad(i,j) = (x_actual(i, j) + x_actual(i, j + 1) + x_actual(i + 1, j) + x_actual(i + 1, j + 1)) / 4;
//        y_grad(i,j) = (y_actual(i, j) + y_actual(i, j + 1) + y_actual(i + 1, j) + y_actual(i + 1, j + 1)) / 4;    
//    end
//end
//
//grad_mag = (grad_zu.^2 + grad_zv.^2).^(1/2) * 50 * 10^3;
////grad_mag = grad_zv * (50 * 10^3);
//slope_error = scf(6);
//xset('colormap', jetcolormap(128));
//[xf, yf, zf] = nf3d(x_grad', y_grad', grad_mag'); 
//plot3d1(xf, yf, zf, theta = 300, alpha = 60, leg = "@@", flag = [-1 6 4]);
//colorbar(min(grad_mag), max(grad_mag),,fmt="%.3f");
//
//e = gce();
//e.parent.title.text = "$\Large ||{\nabla Z_{err}}|| \tt\Big(\frac{\mu m}{50 mm}\Big)$";
//
//a = gca();
//a.view = "2d";
//a.axes_reverse = ["off", "on", "on"];
//
//format('v', 10)
//a.x_label.text = "$\Large X \tt(mm)$";
//a.y_label.text = "$\Large Y \tt(mm)$";
//a.title.text = "$\text{\begin{gather}\Huge{Slope \ Error}$" + "$\\ \LARGE{' + sPartNumber  +' \ SN" + sSerialNumber + "}$" +...
//    "$\\ \LARGE{R_{bf} = " + string(param_opt(4)) + " mm}\end{gather}}$";
//
//facets_a = a.children(1);
//facets_a.hiddencolor = -1;
//
//mkdir('/' + sSerialNumber);
//mkdir('/' + sSerialNumber + '/imgs');
//xs2pdf(gcf(), pwd() + '\' + string(sSerialNumber) + '\imgs\' + string(sPartNumber) + '_'...
//    + string(sSerialNumber) + '_SLOPE_ERROR_FIT_XYZR.pdf');

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
    
vTS = datevec(now());
sTimeStampTemp = string(vTS(1:5));

for i = 1 : size(sTimeStampTemp, "*")
    sTimeStamp(i) = strcat(tokens(sTimeStampTemp(i), " "), "");
end

/////////////////////////////////////////////Write RAW points to output folder

sFileToSave = 'Surface_Points_' + [sTimeStamp(1)+sTimeStamp(2)+sTimeStamp(3)+sTimeStamp(4)+sTimeStamp(5)] + '.csv';
sTempFile = TMPDIR + "\" + sFileToSave;
fFile = mopen(sTempFile, 'wt');

for i=1:size(M_Actual, 1) - 1
    mfprintf(fFile, "%f.6, %f.6, %f.6, %f.6, %f.6, %f.6\n", M_Actual(i, 1), M_Actual(i, 2), M_Actual(i, 3), M_Actual(i, 4), M_Actual(i, 5), M_Actual(i, 6));
end
mfprintf(fFile, "%f.6, %f.6, %f.6, %f.6, %f.6, %f.6", M_Actual($, 1), M_Actual($, 2), M_Actual($, 3), M_Actual($, 4), M_Actual($, 5), M_Actual($, 6));
mclose(fFile);

dos('move ' + sTempFile + ' ' + pwd() + '\' + string(sSerialNumber) + '\' + sFileToSave);
