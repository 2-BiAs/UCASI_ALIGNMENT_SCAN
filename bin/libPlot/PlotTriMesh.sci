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
            pause;
        end
    end
    
    figCurrent = scf();
    
    figCurrent.color_map = jetcolormap(128);
    
    axCurrent = gca();
    axCurrent.isoview = "on"
    
    plot3d(X, Y, Z);
    
    eSurf = gce();
     
    
    eSurf.color_flag=1; //color according to z
    //eSurf.color_mode=-2;  //remove the facets boundary by setting color_mode to white color
//    eSurf.color_flag=2; //color according to given colors
//    eSurf.color_mode = -1; // put the facets boundary back by setting
    
//    eSurf.color_flag=3; // interpolated shading mode
    
    
    
    axOutput = gca();
        
endfunction
