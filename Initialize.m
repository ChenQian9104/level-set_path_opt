%% initialize sacnning path
function [ ScanLSgridPhi ] = Initialize(nelx,nely,nelz)
[x1, y1, z1] = meshgrid( -0.5:nelx+0.5, -0.5:nely+0.5, 0.5:nelz-0.5 );   % meshgrid for fixed geometry
ScanLSgridX = x1(:);
ScanLSgridY = y1(:);
ScanLSgridZ = z1(:);

%cx = [ -0.5  10  20   30  40   50   60   70   80   90   100.5    -0.5  10  20   30  40   50   60   70   80   90   100.5  -0.5  10  20   30  40   50   60   70   80   90   100.5];
%cy = [   10  10  10   0   0     0   0   10   10   10    10       25   25  25   25  25   25   25   25   25   25    25     40   40  40   40   50   50   50   40   40   40    40];

cx = [   5  15  30  45  55       10  30  50       5  15  30  45  55 ];
cy = [   10 10  10  10  10       15  15  15       25 25  25  25  25   ];

for i = 1:length(cx)
    temPhi(:,i) = sqrt( ( ScanLSgridX(:) - cx(i) ).^2 +( ScanLSgridY(:) - cy(i) ).^2 ) - 2;   
end

LSgridPhi_Scan =  min( temPhi.').';
LSgridPhi_Scan((ScanLSgridX - min(ScanLSgridX))  .* (ScanLSgridX - max(ScanLSgridX))  .* (ScanLSgridY - max(ScanLSgridY)) .* (ScanLSgridY - min(ScanLSgridY)) <= 100*eps) = 5;

ScanLSgridPhi = reshape(LSgridPhi_Scan,[nely+2,nelx+2,nelz]);   % level set function for fixed geometry
