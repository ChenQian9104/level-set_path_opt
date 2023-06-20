function [Theta1,Theta2] = ThetaScanningPath(ScanLSgridPhi,dx,dy)
nely = size( ScanLSgridPhi,1)-2;
nelx = size( ScanLSgridPhi,2)-2;
nelz = size( ScanLSgridPhi,3);
Theta1 = zeros(nely,nelx,nelz);
Theta2 = zeros(nely+2,nelx+2,nelz);


for k = 1:nelz
    Matrix = Matrix4diff( ScanLSgridPhi(:,:,k ) );
    dphidy = ( Matrix.j_plus_1(2:end-1,2:end-1) - Matrix.j_minus_1(2:end-1,2:end-1) )/dy;
    dphidx = ( Matrix.i_plus_1(2:end-1,2:end-1) - Matrix.i_minus_1(2:end-1,2:end-1) )/dx;
    for j = 1:nely
        for i = 1:nelx
            if dphidx(j,i) == 0
                Theta1(j,i,k) = 0;
            else
                Theta1(j,i,k) = pi/2 + atan( dphidy(j,i)/dphidx(j,i) );
            end
        end
    end   
end
Theta2(2:end-1,2:end-1,:) = Theta1;

