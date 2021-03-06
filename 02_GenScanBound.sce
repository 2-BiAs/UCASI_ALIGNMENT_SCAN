//Generate uCasi Mirror Scan Boundary

alpha = 54.7*%pi/180;
theta = [linspace(%pi/2 - alpha, %pi/2 + alpha, 8)' ; linspace(3*%pi/2 - alpha, 3*%pi/2 + alpha, 8)'];
R = 26;

mPoints = [R * cos(theta), R * sin(theta)];
mPoints2 = mPoints ./ 2;

bmMask = PolygonMask(mPoints, 'INT'); //Generate intersection mask with the filled backround
bmMask2 = PolygonMask(mPoints2, 'DIFF');

meshGrid = PolyGrid(list(bmMask, bmMask2), 1);

spSurfParams = SurfaceParameters(100, 0, [0]);
fSags = GetSags(meshGrid.mVertices, spSurfParams);

meshGrid.mVertices(:, 3) = fSags(:);

PlotTriMesh(meshGrid);

//plot(mGridPoints(:,1), mGridPoints(:,2), '*');
//axesMain = gca();
//axesMain.isoview = "on"
