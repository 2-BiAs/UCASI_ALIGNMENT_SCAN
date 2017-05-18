function axOutput = PlotTriMesh(TriMesh)
    
    iTriCount = size(TriMesh.iIndices, 1);
    
    
    for i = 1:iTriCount
        
        //Indices for single triangle
        iT = TriMesh.iIndices(i, 1:3);
        
        //mprintf("i = %d\n", i)
        try
            X(1:3, $+1) = [TriMesh.mVertices(iT(1),1); TriMesh.mVertices(iT(2),1); TriMesh.mVertices(iT(3),1)];
            Y(1:3, $+1) = [TriMesh.mVertices(iT(1),2); TriMesh.mVertices(iT(2),2); TriMesh.mVertices(iT(3),2)];
            Z(1:3, $+1) = [TriMesh.mVertices(iT(1),3); TriMesh.mVertices(iT(2),3); TriMesh.mVertices(iT(3),3)];
        catch
            printf('I be pausin here (in PlotTriMesh)');
            pause;
        end
    end
    
    figCurrent = scf();
    
    figCurrent.color_map = jetcolormap(128);
    
    axCurrent = gca();
    axCurrent.isoview = "on"
    
    plot3d(X, Y, Z);
    
    disp('fart')
    eSurf = gce();
    //pause
    tlistSurfData=eSurf.data;
    tlistNewData = tlist(["3d" "x" "y" "z" "color"], tlistSurfData.x, tlistSurfData.y, tlistSurfData.z, tlistSurfData.z); 
    eSurf.data = tlistNewData;
    eSurf.color_mode=-1;  //remove the facets boundary
    eSurf.color_flag=3;   // interpolated shading mode
    eSurf.cdata_mapping = 'scaled';

    
    axOutput = gca();
        
endfunction
