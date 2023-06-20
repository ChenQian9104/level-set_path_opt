function [shapeSens] = Sensitivity( Be,C,struc,U,edofMat,depsilon,Dtheta_dphi_i_minus_1,Dtheta_dphi_i_plus_1,Dtheta_dphi_j_minus_1,Dtheta_dphi_j_plus_1)
nely = size( Dtheta_dphi_i_minus_1,1);
nelx = size( Dtheta_dphi_i_minus_1,2);
nelz = size( Dtheta_dphi_i_minus_1,3);

shapeSens = zeros( nely, nelx, nelz);

ele = 0;
%=== Sensitivity for the first layer scanning path level-set function ====%
for elz = 1:nelz
    for ely = 1:nely
        for elx = 1:nelx
            ele = ele + 1;
            
            if elx >= 2
                w1 = 2*U( edofMat( ele-1,:)');
                term1 = struc(ely,elx-1,elz)*w1'*Be'*C*depsilon(:,ele-1)*Dtheta_dphi_i_minus_1(ely,elx,elz);
            else
                term1 = 0;
            end
            
            if elx <= nelx-1
                w1 = 2*U( edofMat(ele+1,:)' );
                term2 = struc( ely, elx+1, elz)*w1'*Be'*C*depsilon(:,ele+1)*Dtheta_dphi_i_plus_1(ely,elx,elz);
            else
                term2 = 0;
            end
            
            if ely >= 2
                w1 = 2*U( edofMat(ele-nelx,:)' );
                term3 = struc(ely-1,elx,elz)*w1'*Be'*C*depsilon(:,ele-nelx)*Dtheta_dphi_j_minus_1(ely,elx,elz);
            else
                term3 = 0;
            end
            
            if ely <= nely - 1
                w1 = 2*U( edofMat( ele+nelx,: )' );
                term4 = struc(ely+1,elx,elz)*w1'*Be'*C*depsilon(:,ele+nelx)*Dtheta_dphi_j_plus_1(ely,elx,elz);
            else
                term4 = 0;
            end
            shapeSens(ely,elx,elz) = ( term1 + term2 + term3 + term4 );
        end
    end
end



%{
for ely = 1:nely
    for elx = 1:nelx 
        ele = ele + 1;
        
        if elx >= 2
            w1 = W( edofMat(ele-1,:)' );  % adjoint variables
            depsilon = [ depx_dtheta1(ely, elx-1,1); depy_dtheta1(ely,elx-1,1); depz_dtheta1(ely,elx-1,1);0;0;0];
            term1 = struc(ely,elx-1,1)*w1'*Be'*C*depsilon*Dtheta_dphi_i_minus_1(ely,elx,1);
        else
            term1 = 0;
        end
        
        if elx <= nelx-1
            w1 = W( edofMat(ele+1,:)' );  % adjoint variables
            depsilon = [ depx_dtheta1(ely, elx+1,1); depy_dtheta1(ely,elx+1,1); depz_dtheta1(ely,elx+1,1);0;0;0];
            term2 = struc(ely,elx+1,1)*w1'*Be'*C*depsilon*Dtheta_dphi_i_plus_1(ely,elx,1);            
        else
            term2 = 0;
        end
        
        if ely >= 2
            w1 = W( edofMat(ele-nelx,:)' );  % adjoint variables
            depsilon = [ depx_dtheta1(ely-1, elx,1); depy_dtheta1(ely-1,elx,1); depz_dtheta1(ely-1,elx,1);0;0;0];
            term3 = struc(ely-1,elx,1)*w1'*Be'*C*depsilon*Dtheta_dphi_j_minus_1(ely,elx,1);            
        else
            term3 = 0;
        end

        if ely <= nely-1
            w1 = W( edofMat(ele+nelx,:)' );  % adjoint variables
            depsilon = [ depx_dtheta1(ely+1, elx,1); depy_dtheta1(ely+1,elx,1); depz_dtheta1(ely+1,elx,1);0;0;0];
            term4 = struc(ely+1,elx,1)*w1'*Be'*C*depsilon*Dtheta_dphi_j_plus_1(ely,elx,1);            
        else
            term4 = 0;
        end        
        
        shapeSens(ely,elx,1) = -( term1 + term2 + term3 + term4 );
        
    end
end

for elz = 2:nelz
    for ely = 1:nely
        for elx = 1:nelx
            ele = ele + 1;
            
            if elx >= 2
                w1 = W( edofMat(ele-1-nelx*nely,:)' );
                w2 = W( edofMat(ele-1,:)' );
                depsilon1 = [ depx_dtheta2(ely, elx-1,elz-1); depy_dtheta2(ely,elx-1,elz-1); depz_dtheta2(ely,elx-1,elz-1);0;0;0];
                depsilon2 = [ depx_dtheta1(ely, elx-1,elz); depy_dtheta1(ely,elx-1,elz); depz_dtheta1(ely,elx-1,elz);0;0;0];
                term1 = ( struc(ely,elx-1,elz-1)*w1'*Be'*C*depsilon1 + struc(ely,elx-1,elz)*w2'*Be'*C*depsilon2)*Dtheta_dphi_i_minus_1(ely,elx,elz);
            else
                term1 = 0;
            end
            
            if elx <= nelx-1
                w1 = W( edofMat(ele+1-nelx*nely,:)' );
                w2 = W( edofMat(ele+1,:)' );
                depsilon1 = [ depx_dtheta2(ely, elx+1,elz-1); depy_dtheta2(ely,elx+1,elz-1); depz_dtheta2(ely,elx+1,elz-1);0;0;0];
                depsilon2 = [ depx_dtheta1(ely, elx+1,elz); depy_dtheta1(ely,elx+1,elz); depz_dtheta1(ely,elx+1,elz);0;0;0];
                term2 = ( struc(ely,elx+1,elz-1)*w1'*Be'*C*depsilon1 + struc(ely,elx+1,elz)*w2'*Be'*C*depsilon2)*Dtheta_dphi_i_plus_1(ely,elx,elz);
            else
                term2 = 0;
            end
            
            if ely >=2
                w1 = W( edofMat(ele-nelx-nelx*nely,:)' );
                w2 = W( edofMat(ele-nelx,:)' );
                depsilon1 = [ depx_dtheta2(ely-1, elx,elz-1); depy_dtheta2(ely-1,elx,elz-1); depz_dtheta2(ely-1,elx,elz-1);0;0;0];
                depsilon2 = [ depx_dtheta1(ely-1, elx,elz); depy_dtheta1(ely-1,elx,elz); depz_dtheta1(ely-1,elx,elz);0;0;0];
                term3 = ( struc(ely-1,elx,elz-1)*w1'*Be'*C*depsilon1 + struc(ely-1,elx,elz)*w2'*Be'*C*depsilon2)*Dtheta_dphi_j_minus_1(ely,elx,elz);                
            else
                term3 = 0;
            end
            
            if ely <= nely-1
                w1 = W( edofMat(ele+nelx-nelx*nely,:)' );
                w2 = W( edofMat(ele+nelx,:)' );
                depsilon1 = [ depx_dtheta2(ely+1, elx,elz-1); depy_dtheta2(ely+1,elx,elz-1); depz_dtheta2(ely+1,elx,elz-1);0;0;0];
                depsilon2 = [ depx_dtheta1(ely+1, elx,elz); depy_dtheta1(ely+1,elx,elz); depz_dtheta1(ely+1,elx,elz);0;0;0];
                term4 = ( struc(ely+1,elx,elz-1)*w1'*Be'*C*depsilon1 + struc(ely+1,elx,elz)*w2'*Be'*C*depsilon2)*Dtheta_dphi_j_plus_1(ely,elx,elz);                
            else
                term4 = 0;
            end
            
            shapeSens(ely,elx,elz) = -( term1 + term2 + term3 + term4 );
            
        end
    end
end
%}