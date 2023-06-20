function [ ScanLSgridPhi ] = Initialize2(nelx,nely,nelz)

ScanLSgridPhi = ones( nely + 2,nelx+2,nelz);
ScanLSgridPhi( :, 2:2:nelx+2,:) = 0;

m = size( 2:2:nelx+2,2);

for i = 1:m
    ScanLSgridPhi( :, (i-1)*2+1 ,: ) = (-1)^i*rand(nely+2,1,nelz);
end

for i = 1:nelz
    ScanLSgridPhi = Reinitialize( ScanLSgridPhi,1,1,80);
    
end
