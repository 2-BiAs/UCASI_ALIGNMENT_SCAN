//Generate Mesh For 16010-4002-001

theta = linspace(0, 2 * %pi, 25)';
theta($) = [];

OD = 70;
ID = 28;

mPoints = [OD / 2 * cos(theta), OD / 2 * sin(theta)];
mPoints2 = [ID / 2 * cos(theta), ID / 2 * sin(theta)];

bmMask = PolygonMask(mPoints, 'INT'); //Generate intersection mask with the filled backround
bmMask2 = PolygonMask(mPoints2, 'DIFF');


meshGrid = PolyGrid(list(bmMask, bmMask2), 1);

spSurfParams = SurfaceParameters(-38.448, -1.374, [0]);
fSags = GetSags(meshGrid.mVertices, spSurfParams);

meshGrid.mVertices(:, 3) = fSags(:);

PlotTriMesh(meshGrid);

//plot(mGridPoints(:,1), mGridPoints(:,2), '*');
axesMain = gca();
axesMain.isoview = "on"
