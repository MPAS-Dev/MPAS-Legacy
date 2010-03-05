

% This script open a grid.nc file and writes out the grid description
% in ascii format to be read in by OpenDX

clear all

% begin periodic parameters
doPeriodic = 1
dc = 1000.0
nx = 200
ny = 200
% end periodic parameters

doWrite = 1
doVor = 0
doTri = 0
doEdge = 1

ncid = netcdf.open('../grid.nc','nc_nowrite');

if (doVor == 1)
    
    xV_id = netcdf.inqVarID(ncid,'xVertex');
    yV_id = netcdf.inqVarID(ncid,'yVertex');
    zV_id = netcdf.inqVarID(ncid,'zVertex');
    nEdgesOnCell_id = netcdf.inqVarID(ncid,'nEdgesOnCell');
    verticesOnCell_id = netcdf.inqVarID(ncid,'verticesOnCell');
    areaCell_id = netcdf.inqVarID(ncid,'areaCell');

    xV=netcdf.getVar(ncid, xV_id);
    yV=netcdf.getVar(ncid, yV_id);
    zV=netcdf.getVar(ncid, zV_id);
    nEdgesOnCell=netcdf.getVar(ncid, nEdgesOnCell_id);
    verticesOnCell=netcdf.getVar(ncid, verticesOnCell_id);
    areaCell = netcdf.getVar(ncid, areaCell_id);

    xC_id = netcdf.inqVarID(ncid,'xCell');
    yC_id = netcdf.inqVarID(ncid,'yCell');
    zC_id = netcdf.inqVarID(ncid,'zCell');

    xC=netcdf.getVar(ncid, xC_id);
    yC=netcdf.getVar(ncid, yC_id);
    zC=netcdf.getVar(ncid, zC_id);
    
    work=size(nEdgesOnCell(:,1));
    nCells=work(1)

    if (doWrite == 1)
    system('rm -f ../dx/vor.position.data');
    system('rm -f ../dx/vor.edge.data');
    system('rm -f ../dx/vor.loop.data');
    system('rm -f ../dx/vor.face.data');
    system('rm -f ../dx/vor.area.data');

    iloop=0;
    iedge=0;
    for i=1:nCells
     dlmwrite('../dx/vor.face.data', i-1, '-append');
     dlmwrite('../dx/vor.area.data', areaCell(i), ...
        'precision', '%18.10e', '-append');
       dlmwrite('../dx/vor.loop.data', iloop, ...
        'precision', '%10i', '-append');
       edge(1:nEdgesOnCell(i)) = iedge;
     
     for j=1:nEdgesOnCell(i)
         x(1,j) = xV(verticesOnCell(j,i));
         x(2,j) = yV(verticesOnCell(j,i));
         x(3,j) = zV(verticesOnCell(j,i));
     end;
     
     if (doPeriodic == 1);
         for j=1:nEdgesOnCell(i);
             dx = x(1,j)-xC(i);
             dy = x(2,j)-yC(i);
             if(abs(dx) > 0.1*nx*dc);
                 if(dx > 0);, x(1,j) = x(1,j) - nx*dc;, end;
                 if(dx < 0);, x(1,j) = x(1,j) + nx*dc;, end;
             end;
             if(abs(dy) > 0.1*ny*dc*sqrt(3)/2);
                 if(dy > 0);, x(2,j) = x(2,j) - sqrt(3)*nx*dc/2;, end;
                 if(dy < 0);, x(2,j) = x(2,j) + sqrt(3)*nx*dc/2;, end;
             end;
         end;
     end;
     
     for j=1:nEdgesOnCell(i)
         dlmwrite('../dx/vor.position.data', x(:,j), 'delimiter', '\t', ...
             'precision', '%18.10e', '-append');
         edge(j) = iedge + j - 1;
       end;
       dlmwrite('../dx/vor.edge.data', edge(1:nEdgesOnCell(i)), ...
        'delimiter', '\t', 'precision', '%10i', '-append')
       iloop = iloop + nEdgesOnCell(i);
      iedge = iedge + nEdgesOnCell(i);
    end;
    
    end;

end;

if (doTri == 1)

    xC_id = netcdf.inqVarID(ncid,'xCell');
    yC_id = netcdf.inqVarID(ncid,'yCell');
    zC_id = netcdf.inqVarID(ncid,'zCell');
    nCellsOnVertex = 3;
    cellsOnVertex_id = netcdf.inqVarID(ncid, 'cellsOnVertex');
    areaTriangle_id = netcdf.inqVarID(ncid,'areaTriangle');

    xC=netcdf.getVar(ncid, xC_id);
    yC=netcdf.getVar(ncid, yC_id);
    zC=netcdf.getVar(ncid, zC_id);
    cellsOnVertex=netcdf.getVar(ncid, cellsOnVertex_id);
    areaTriangle = netcdf.getVar(ncid, areaTriangle_id);
    
    xV_id = netcdf.inqVarID(ncid,'xVertex');
    yV_id = netcdf.inqVarID(ncid,'yVertex');
    zV_id = netcdf.inqVarID(ncid,'zVertex');

    xV=netcdf.getVar(ncid, xV_id);
    yV=netcdf.getVar(ncid, yV_id);
    zV=netcdf.getVar(ncid, zV_id);

    work=size(cellsOnVertex);
    nVertices = work(:,2)

    if (doWrite == 1)
    system('rm -f ../dx/tri.position.data');
    system('rm -f ../dx/tri.edge.data');
    system('rm -f ../dx/tri.loop.data');
    system('rm -f ../dx/tri.face.data');
    system('rm -f ../dx/tri.area.data');
    
    iloop=0;
    iedge=0;
    for i=1:nVertices
     dlmwrite('../dx/tri.face.data', i-1, '-append');
     dlmwrite('../dx/tri.area.data', areaTriangle(i), ...
        'precision', '%18.10e', '-append');
     dlmwrite('../dx/tri.loop.data', iloop, ...
        'precision', '%10i', '-append');
     edge(1:3) = iedge;
     for j=1:nCellsOnVertex
         x(1,j) = xC(cellsOnVertex(j,i));
         x(2,j) = yC(cellsOnVertex(j,i));
         x(3,j) = zC(cellsOnVertex(j,i));
     end;
     
     if (doPeriodic == 1);
         for j=1:nCellsOnVertex;
             dx = x(1,j)-xV(i);
             dy = x(2,j)-yV(i);
             if(abs(dx) > 0.1*nx*dc);
                 if(dx > 0);, x(1,j) = x(1,j) - nx*dc;, end;
                 if(dx < 0);, x(1,j) = x(1,j) + nx*dc;, end;
             end;
             if(abs(dy) > 0.1*ny*dc*sqrt(3)/2);
                 if(dy > 0);, x(2,j) = x(2,j) - sqrt(3)*nx*dc/2;, end;
                 if(dy < 0);, x(2,j) = x(2,j) + sqrt(3)*nx*dc/2;, end;
             end;
         end;
     end;
     
     for j=1:nCellsOnVertex;
         dlmwrite('../dx/tri.position.data', x(:,j), 'delimiter', '\t', ...
             'precision', '%18.10e', '-append')
         edge(j) = iedge + j - 1;
     end;
     dlmwrite('../dx/tri.edge.data', edge(1:3), ...
         'delimiter', '\t', 'precision', '%10i', '-append')
     iloop = iloop + 3;
     iedge = iedge + 3;
    end;
    
    end;

end;

if (doEdge == 1)
    
    if (doWrite == 1)
        system('rm -f ../dx/edge.position.data');
        system('rm -f ../dx/normal.data');
        system('rm -f ../dx/tangent.data');
    end;
    
    [nEdgesName, nEdgesLength] = netcdf.inqDim(ncid,1);
    xE_id = netcdf.inqVarID(ncid,'xEdge');
    yE_id = netcdf.inqVarID(ncid,'yEdge');
    zE_id = netcdf.inqVarID(ncid,'zEdge');
    nCellsOnEdge = 2;
    nEdges = nEdgesLength;
    cellsOnEdge_id = netcdf.inqVarID(ncid, 'cellsOnEdge');

    xE=netcdf.getVar(ncid, xE_id);
    yE=netcdf.getVar(ncid, yE_id);
    zE=netcdf.getVar(ncid, zE_id);
    cellsOnEdge=netcdf.getVar(ncid, cellsOnEdge_id);
    
    xC_id = netcdf.inqVarID(ncid,'xCell');
    yC_id = netcdf.inqVarID(ncid,'yCell');
    zC_id = netcdf.inqVarID(ncid,'zCell');

    xC=netcdf.getVar(ncid, xC_id);
    yC=netcdf.getVar(ncid, yC_id);
    zC=netcdf.getVar(ncid, zC_id);
    
    for i=1:nEdges
      
        j1 = cellsOnEdge(1,i);
        j2 = cellsOnEdge(2,i);
        iCell1 = min(j1,j2);
        iCell2 = max(j1,j2);
        
        x(1) = xC(iCell2)-xC(iCell1);
        x(2) = yC(iCell2)-yC(iCell1);
        x(3) = zC(iCell2)-zC(iCell1);
        
        normal = x ./ sqrt(x(1)^2 + x(2)^2 + x(3)^2);
        
        x(1) = xE(i); x(2) = yE(i); x(3) = zE(i);
        tangent(1) = x(2).*normal(3) - x(3).*normal(2);
        tangent(2) = x(3).*normal(1) - x(1).*normal(3);
        tangent(3) = x(1).*normal(2) - x(2).*normal(1);
    
    
        if (doWrite == 1)
        
            dlmwrite('../dx/edge.position.data', xE(i), ...
             'precision', '%18.10e', '-append')
       
            dlmwrite('../dx/edge.position.data', yE(i), ...
             'precision', '%18.10e', '-append')
       
            dlmwrite('../dx/edge.position.data', zE(i), ...
             'precision', '%18.10e', '-append')
       
            dlmwrite('../dx/normal.data', normal(1), ...
             'precision', '%18.10e', '-append')
       
            dlmwrite('../dx/normal.data', normal(2), ...
             'precision', '%18.10e', '-append')
       
            dlmwrite('../dx/normal.data', normal(3), ...
             'precision', '%18.10e', '-append')
       
            dlmwrite('../dx/tangent.data', tangent(1), ...
             'precision', '%18.10e', '-append')
       
            dlmwrite('../dx/tangent.data', tangent(2), ...
             'precision', '%18.10e', '-append')
       
            dlmwrite('../dx/tangent.data', tangent(3), ...
             'precision', '%18.10e', '-append')
    
    end;
    
end;    
    
end;

netcdf.close(ncid)