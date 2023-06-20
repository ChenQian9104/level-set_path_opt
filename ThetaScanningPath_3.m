function [Theta1,Theta2,Dtheta_dphi_i_minus_1,Dtheta_dphi_i_plus_1,Dtheta_dphi_j_minus_1,Dtheta_dphi_j_plus_1] = ThetaScanningPath_3(ScanLSgridPhi,dx,dy)
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
            if abs( dphidx(j,i) ) <= 1e-3
                if dphidy(j,i) >=0
                    Theta2(j,i,k) = pi;
                else
                    Theta2(j,i,k) = 0;
                end
            else
                Theta2(j,i,k) = pi/2 + atan( dphidy(j,i)/dphidx(j,i) );
            end
            
            if i > 2 
                
                if abs( dphidx(j,i-1) ) <= 1e-1
                    dtheta_dphi_i_minus_1(j,i,k) = 0;
                else
                    A  = dphidy(j,i-1)/dphidx(j,i-1);
                    %dtheta_dphi_i_minus_1(j,i,k) = -1/( 1 + A^2 )*A^2/dphidy(j,i-1);
                    dtheta_dphi_i_minus_1(j,i,k) = 1/( 1 + A^2 )*-A/dphidx(j,i-1);
                    if abs( dtheta_dphi_i_minus_1(j,i,k) ) >=100
                       dtheta_dphi_i_minus_1(j,i,k) = sign( dtheta_dphi_i_minus_1(j,i,k) )*100;
                    end
                end
                
            end
            
            if i < nelx-1
                
                if abs( dphidx(j,i+1) ) <= 1e-1
                    dtheta_dphi_i_plus_1(j,i,k) = 0;
                else
                    A  = dphidy(j,i+1)/dphidx(j,i+1);
                    %dtheta_dphi_i_plus_1(j,i,k) = 1/( 1 + A^2 )*A^2/dphidy(j,i+1);
                    dtheta_dphi_i_plus_1(j,i,k) = 1/( 1 + A^2 )*A/dphidx(j,i+1);
                    if abs( dtheta_dphi_i_plus_1(j,i,k) ) >=100
                       dtheta_dphi_i_plus_1(j,i,k) = sign(  dtheta_dphi_i_plus_1(j,i,k) )*100;
                    end
                end                
                             
            end
            
            
            if j > 2 
                if abs( dphidx(j-1,i) ) <= 1e-1
                    dtheta_dphi_j_minus_1(j,i,k) = 0;
                else
                    A = dphidy(j-1,i)/dphidx(j-1,i);
                    dtheta_dphi_j_minus_1(j,i,k) = 1/( 1 + A^2 )/dphidx(j-1,i);
                    if abs( dtheta_dphi_j_minus_1(j,i,k) ) >=100
                       dtheta_dphi_j_minus_1(j,i,k) = sign(   dtheta_dphi_j_minus_1(j,i,k) )*100;
                    end                    
                end
                
            end
            
            if j < nely-1
                if abs( dphidx(j+1,i) ) <= 1e-1
                    dtheta_dphi_j_plus_1(j,i,k) = 0;
                else
                    A = dphidy(j+1,i)/dphidx(j+1,i);
                    dtheta_dphi_j_plus_1(j,i,k) = -1/( 1 + A^2 )/dphidx(j+1,i);
                    if abs( dtheta_dphi_j_plus_1(j,i,k) ) >=100
                       dtheta_dphi_j_plus_1(j,i,k) = sign(   dtheta_dphi_j_plus_1(j,i,k) )*100;
                    end                      
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


%=======  Filter out unreasonable value ======%
for k = 1:nelz 
    median1 = median( median( abs(Dtheta_dphi_i_minus_1(:,:,k) ) ) );
    median2 = median( median( abs(Dtheta_dphi_i_plus_1(:,:,k) ) ) );
    median3 = median( median( abs(Dtheta_dphi_j_minus_1(:,:,k) ) ) );
    median4 = median( median( abs(Dtheta_dphi_j_plus_1(:,:,k) ) ) );
    for j = 1:nely-2
        
        for i = 1:nelx-2
            if abs( Dtheta_dphi_i_minus_1(j,i,k)/median1 ) >= 5
                Dtheta_dphi_i_minus_1(j,i,k) = 0;
            end
            
            if abs( Dtheta_dphi_i_plus_1(j,i,k)/median2 ) >= 5
                Dtheta_dphi_i_plus_1(j,i,k) = 0;
            end
            
            if abs( Dtheta_dphi_j_minus_1(j,i,k)/median3 ) >= 5
                Dtheta_dphi_j_minus_1(j,i,k) = 0;
            end
            
            if abs( Dtheta_dphi_j_plus_1(j,i,k)/median4 ) >= 5
                Dtheta_dphi_j_plus_1(j,i,k) = 0;
            end            
            
        end
    end 
end




