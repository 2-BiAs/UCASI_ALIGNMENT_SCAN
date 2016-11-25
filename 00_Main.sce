
exec('01_Initialize.sce');

//sp = SurfaceParameters(100, 0, [0])
//z = GetSags([-5:.1:5, -5:.1:5], sp)

[X Y] = ndgrid([-5:1:5], [-5:1:5]);

P = GridToList(X, Y)
