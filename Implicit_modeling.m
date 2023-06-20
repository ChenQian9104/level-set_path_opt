clc;clear;
nelx = 60; nely = 30; dx = 1; dy = 1;
nelz_base = 5; nelz_bracket = 20; nelz_cylinder = 11;
nelz = nelz_base + nelz_bracket + nelz_cylinder;
nele = nelx*nely*nelz;
[x y z] = meshgrid( -0.5:nelx + 0.5, -0.5:nely + 0.5,1:nelz);
r = 4;
%% Define component level set
GeoLSgridPhi_base = min( min( min( x, nelx - x), min( y, nely-y) ), min( z-1, nelz_base - z) );
GeoLSgridPhi_hole1 = min( ( sqrt( (x-6).^2 + (y-6).^2 ) - r ), min( z-1, nelz_base - z ) );
GeoLSgridPhi_hole2 = min( ( sqrt( (x-54).^2 + (y-6).^2)- r  ), min( z-1, nelz_base - z ) );
GeoLSgridPhi_hole3 = min( ( sqrt( (x-54).^2 + (y-24).^2) - r ), min(z-1, nelz_base - z ) );
GeoLSgridPhi_hole4 = min( ( sqrt( (x-6).^2 + (y-24).^2) - r), min( z-1, nelz_base - z ) );

GeoLSgridPhi_boss1 = min( min( min( x - 20, 24 - x), min( y - 4, 26-y) ), min( z-nelz_base, nelz_base + nelz_bracket - z) );
GeoLSgridPhi_boss2 = min( min( min( x - 36, 40 - x), min( y - 4, 26-y) ), min( z-nelz_base, nelz_base + nelz_bracket - z) );


GeoLSgridPhi_cylinder1 = min ( ( 11 - sqrt( (y-15).^2 + (z-25).^2) ), min( x - 20, 24 - x) );
GeoLSgridPhi_cylinder2 = min ( ( 11 - sqrt( (y-15).^2 + (z-25).^2) ), min( x - 36, 40 - x) );

GeoLSgridPhi_hole5 = min ( ( 5-sqrt( (y-15).^2 + (z-25).^2) ), min( x - 20, 24 - x) );
GeoLSgridPhi_hole6 = min ( ( 5-sqrt( (y-15).^2 + (z-25).^2) ), min( x - 36, 40 - x) );


%% Assemble the base plate
GeoLSgridPhi_base = min( GeoLSgridPhi_base, GeoLSgridPhi_hole1 );
GeoLSgridPhi_base = min( GeoLSgridPhi_base, GeoLSgridPhi_hole2 );
GeoLSgridPhi_base = min( GeoLSgridPhi_base, GeoLSgridPhi_hole3 );
GeoLSgridPhi_base = min( GeoLSgridPhi_base, GeoLSgridPhi_hole4 );


GeoLSgridPhi_part = max( GeoLSgridPhi_base, GeoLSgridPhi_boss1);
GeoLSgridPhi_part = max( GeoLSgridPhi_part, GeoLSgridPhi_boss2);

GeoLSgridPhi_part = max( GeoLSgridPhi_part, GeoLSgridPhi_cylinder1);
GeoLSgridPhi_part = max( GeoLSgridPhi_part, GeoLSgridPhi_cylinder2);

GeoLSgridPhi_part = min( GeoLSgridPhi_part, -GeoLSgridPhi_hole5);
GeoLSgridPhi_part = min( GeoLSgridPhi_part, -GeoLSgridPhi_hole6);

GeoLSgridPhi_Assemble = GeoLSgridPhi_part(2:end-1,2:end-1,:);

[Q1,Q2] = meshgrid(1:nelx,1:nely);
S = GeoLSgridPhi_Assemble(:,:,2);
figure(2)
contour(Q1,Q2,S,[0:1:10],'ShowText','on')

struc = zeros(nely,nelx, nelz);
struc( GeoLSgridPhi_Assemble >= 0 ) = 1;



%%  ==================== Finite element method ========================= %%
nodeIncrement = (nelx+1)*(nely+1);
eleIncrement = nelx*nely;

nnode = (nelx+1)*(nely+1)*nelz;

m=0;
for j = 1:nely

    for i = 1:nelx
        m = m+1;       
        node1 = i + (j-1)*(nelx+1);
        node2 = node1+1;
        node3 = node2 + (nelx+1);
        node4 = node1 + (nelx+1);
        node5 = node1 + nodeIncrement;
        node6 = node2 + nodeIncrement;
        node7 = node3 + nodeIncrement;
        node8 = node4 + nodeIncrement;
        edofMat(m,:) = [ node1, node2, node3, node4, node5, node6,node7,node8];
    end
end

for k = 2:nelz
    for i = 1:eleIncrement
        edofMat( i + (k-1)*eleIncrement,:) = edofMat(i + (k-2)*eleIncrement,:) + nodeIncrement;
    end
end

