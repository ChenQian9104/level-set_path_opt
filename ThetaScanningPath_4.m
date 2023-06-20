function [Theta1,Theta2,Dtheta_dphi_i_minus_1,Dtheta_dphi_i_plus_1,Dtheta_dphi_j_minus_1,Dtheta_dphi_j_plus_1] = ThetaScanningPath_4(ScanLSgridPhi,dx,dy)
nely = size( ScanLSgridPhi,1);
nelx = size( ScanLSgridPhi,2);
nelz = size( ScanLSgridPhi,3);
Theta1 = zeros(nely-2,nelx-2,nelz);
Theta2 = zeros(nely,nelx,nelz);

dtheta_dphi_i_minus_1 = zeros(nely,nelx,nelz);    % theta(i-1,j)/phi(i,j)
dtheta_dphi_i_plus_1 = zeros(nely,nelx,nelz);     % theta(i+1,j)/phi(i,j)
dtheta_dphi_j_minus_1 = zeros(nely,nelx,nelz);    % theta(i,j-1)/phi(i,j)
dtheta_dphi_j_plus_1 = zeros(nely,nelx,nelz);     % theta(i,j+1)/phi(i,j)


Dtheta_dphi_i_minus_1 = zeros(nely-2,nelx-2,nelz);    % theta(i-1,j)/phi(i,j)
Dtheta_dphi_i_plus_1 = zeros(nely-2,nelx-2,nelz);     % theta(i+1,j)/phi(i,j)
Dtheta_dphi_j_minus_1 = zeros(nely-2,nelx-2,nelz);    % theta(i,j-1)/phi(i,j)
Dtheta_dphi_j_plus_1 = zeros(nely-2,nelx-2,nelz);     % theta(i,j+1)/phi(i,j)


for k = 1:nelz
    Matrix = Matrix4diff( ScanLSgridPhi(:,:,k ) );
    dphidy = ( Matrix.j_plus_1(:,:) - Matrix.j_minus_1(:,:) )/dy;
    dphidx = ( Matrix.i_plus_1(:,:) - Matrix.i_minus_1(:,:) )/dx;
    for j = 2:nely-1
        for i = 2:nelx-1
           
            Theta2(j,i,k) = pi/2 + atan( dphidy(j,i)/dphidx(j,i) );           
            
            if i > 2
                if dphidx(j,i-1)^2 + dphidy(j,i-1)^2 <= 1e-3
                    dtheta_dphi_i_minus_1(j,i,k) = 0;
                else
                    dtheta_dphi_i_minus_1(j,i,k) = -dphidy(j,i-1)/(dphidx(j,i-1)^2 + dphidy(j,i-1)^2); 
                end
            end
               
            if i < nelx-1
                
                if dphidx(j,i+1)^2 + dphidy(j,i+1)^2 <= 1e-3
                    dtheta_dphi_i_plus_1(j,i,k) = 0;
                else                                        
                    dtheta_dphi_i_plus_1(j,i,k) = dphidy(j,i+1)/(dphidx(j,i+1)^2 + dphidy(j,i+1)^2);                  
                end                
                             
            end
            
            
            if j > 2 
                if dphidx(j-1,i)^2 + dphidy(j-1,i)^2 <= 1e-3
                    dtheta_dphi_j_minus_1(j,i,k) = 0;
                else
                    dtheta_dphi_j_minus_1(j,i,k) = dphidx(j-1,i)/(dphidx(j-1,i)^2 + dphidy(j-1,i)^2);
                end
                
            end
            
            if j < nely-1
                if dphidx(j+1,i)^2 + dphidy(j+1,i)^2 <= 1e-3
                    dtheta_dphi_j_plus_1(j,i,k) = 0;
                else                    
                    dtheta_dphi_j_plus_1(j,i,k) =  -dphidx(j+1,i)/(dphidx(j+1,i)^2 + dphidy(j+1,i)^2);           
                end

            end            
            
        end
    end   
end
Theta1(:,:,:) = Theta2(2:end-1,2:end-1,:);



Dtheta_dphi_i_minus_1(:,:,:) = dtheta_dphi_i_minus_1(2:end-1,2:end-1,:);    
Dtheta_dphi_i_plus_1(:,:,:)  = dtheta_dphi_i_plus_1(2:end-1,2:end-1,:);      
Dtheta_dphi_j_minus_1(:,:,:) = dtheta_dphi_j_minus_1(2:end-1,2:end-1,:);    
Dtheta_dphi_j_plus_1(:,:,:)  = dtheta_dphi_j_plus_1(2:end-1,2:end-1,:);  



