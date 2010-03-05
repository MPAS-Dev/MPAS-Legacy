
% This script open a output.nc file and writes out the output
% in ascii format to be read in by OpenDX

clear all

eps = 1.0e-12

ncid = netcdf.open('../output.nc','nc_nowrite')

doThickness = 1;
doKE = 1;
doVorticity = 1;
doVelocity = 1;

[nCellsName, nCellsLength] = netcdf.inqDim(ncid,1);
[nEdgesName, nEdgesLength] = netcdf.inqDim(ncid,2);
[nVerticesName, nVerticesLength] = netcdf.inqDim(ncid,5);
[nVertLevelsName, nVertLevelsLength] = netcdf.inqDim(ncid,8);
[nTracersName, nTracersLength] = netcdf.inqDim(ncid,9);
[TimeName, TimeLength] = netcdf.inqDim(ncid,0);

thicknessID = netcdf.inqVarID(ncid,'h');
work =  netcdf.getVar(ncid,thicknessID);
[thicknessName,xtype,dimids,natts] = netcdf.inqVar(ncid,thicknessID);
thickness=work;

keID = netcdf.inqVarID(ncid,'ke');
work =  netcdf.getVar(ncid,keID);
[keName,xtype,dimids,natts] = netcdf.inqVar(ncid,keID);
ke=work;

vorticityID = netcdf.inqVarID(ncid,'vorticity');
work =  netcdf.getVar(ncid,vorticityID);
[vorticityName,xtype,dimids,natts] = netcdf.inqVar(ncid,vorticityID);
vorticity=work;

uID = netcdf.inqVarID(ncid,'u');
work =  netcdf.getVar(ncid,uID);
[uName,xtype,dimids,natts] = netcdf.inqVar(ncid,uID);
u=work;

vID = netcdf.inqVarID(ncid,'v');
work =  netcdf.getVar(ncid,vID);
[vName,xtype,dimids,natts] = netcdf.inqVar(ncid,vID);
v=work;


if (doThickness == 1)
system('rm -f ../dx/h.*.*.data')
for iLevel=1:nVertLevelsLength
for iTime=0:TimeLength-1
    stringTime = int2str(iTime)
    stringVert = int2str(iLevel)
    FileName = strcat('../dx/', thicknessName, '.', ...
        stringVert, '.', stringTime, '.', 'data')
    for iCell=1:nCellsLength
      data = thickness(iLevel,iCell,iTime+1);
      if abs(data) < eps, data=0, end
      dlmwrite(FileName, data, ...
         'precision', '%18.10e', '-append')
    end
end
end
end

if (doKE == 1)
system('rm -f ../dx/ke.*.*.data')
for iLevel=1:nVertLevelsLength
for iTime=0:TimeLength-1
    stringTime = int2str(iTime)
    stringVert = int2str(iLevel)
    FileName = strcat('../dx/', keName, '.', ...
        stringVert, '.', stringTime, '.', 'data')
    for iCell=1:nCellsLength
      data = ke(iLevel,iCell,iTime+1);
      if abs(data) < eps, data=0;, end
      dlmwrite(FileName, data, ...
         'precision', '%18.10e', '-append')
    end
end
end
end

if (doVorticity == 1)
system('rm -f ../dx/vorticity.*.*.data')
for iLevel=1:nVertLevelsLength
for iTime=0:TimeLength-1
    stringTime = int2str(iTime)
    stringVert = int2str(iLevel)
    FileName = strcat('../dx/', vorticityName, '.', ...
        stringVert, '.', stringTime, '.', 'data')
    for iVertex=1:nVerticesLength
      data = vorticity(iLevel,iVertex,iTime+1);
      if abs(data) < eps, data=0;, end
      dlmwrite(FileName, data, ...
         'precision', '%18.10e', '-append')
    end
end
end
end

if (doVelocity == 1)
system('rm -f ../dx/u.*.*.data')
for iLevel=1:nVertLevelsLength
for iTime=0:TimeLength-1
    stringTime = int2str(iTime)
    stringVert = int2str(iLevel)
    FileName = strcat('../dx/', uName, '.', ...
        stringVert, '.', stringTime, '.', 'data')
    for iEdge=1:nEdgesLength
      data = u(iLevel,iEdge,iTime+1);
      if abs(data) < eps, data=0;, end;
      dlmwrite(FileName, data, ...
         'precision', '%18.10e', '-append')
    end
end
end
end

if (doVelocity == 1)
system('rm -f ../dx/v.*.*.data')
for iLevel=1:nVertLevelsLength
for iTime=0:TimeLength-1
    stringTime = int2str(iTime)
    stringVert = int2str(iLevel)
    FileName = strcat('../dx/', vName, '.', ...
        stringVert, '.', stringTime, '.', 'data')
    for iEdge=1:nEdgesLength
      data = v(iLevel,iEdge,iTime+1);
      if abs(data) < eps, data=0;, end;
      dlmwrite(FileName, data, ...
         'precision', '%18.10e', '-append')
    end
end
end
end


