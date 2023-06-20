%% ====== Scanning path optimization for bearing bracket =============== %%
clc;clear;

%% ======================== DEFINE THE MATERIAL PROPERTY =============== %%
E = 2100; % Young's modulus
nu = 0.3; % Poisson's ratio
C = zeros(6,6);
E0 = E/(1+nu)/(1-2*nu);
G = E/(1+nu)/2;

C(1,1) = E0*(1-nu);  C(1,2) = E0*nu;    C(1,3) = E0*nu;
C(2,1) = E0*nu;      C(2,2) = E0*(1-nu);C(2,3) = E0*nu;
C(3,1) = E0*nu;      C(3,2) = E0*nu;    C(3,3) = E0*(1-nu);
C(4,4) = G; C(5,5) = G; C(6,6) = G;

%% ====== Implicit modelling of the bearing bracket with base  ========= %%
nelx = 60; nely = 30; dx = 1; dy = 1; L = 1;
nelz_base = 5; nelz_bracket = 20; nelz_cylinder = 11;
nelz = nelz_base + nelz_bracket + nelz_cylinder;
nele = nelx*nely*nelz;
[x y z] = meshgrid( -0.5:nelx + 0.5, -0.5:nely + 0.5,0.5:nelz-0.5);
r = 4;
GeoLSgridPhi_base = min( min( min( x, nelx - x), min( y, nely-y) ), min( z, nelz_base - z) );
GeoLSgridPhi_hole1 = min( ( sqrt( (x-6).^2 + (y-6).^2 ) - r ), min( z, nelz_base - z ) );
GeoLSgridPhi_hole2 = min( ( sqrt( (x-54).^2 + (y-6).^2)- r  ), min( z, nelz_base - z ) );
GeoLSgridPhi_hole3 = min( ( sqrt( (x-54).^2 + (y-24).^2) - r ), min(z, nelz_base - z ) );
GeoLSgridPhi_hole4 = min( ( sqrt( (x-6).^2 + (y-24).^2) - r), min( z, nelz_base - z ) );

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

struc = zeros(nely,nelx, nelz);
struc( GeoLSgridPhi_Assemble >= 0 ) = 1;



%%  ============= Finite element and node information ================== %%
nodeIncrement = (nelx+1)*(nely+1);
eleIncrement = nelx*nely;
nnode = ( nelx + 1 )*( nely + 1 )*( nelz +  1);
ndofs = 3*nnode;

Be = Be_3D(L);
edofMat = zeros(nele,8);  % connectivity matrix: [node1 node2 node3 node4 .... node8] - element 1

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

edofMat1 = zeros(nele,24);
for i = 1:8
    edofMat1(:,3*i-2) = 3*edofMat(:,i)-2;
    edofMat1(:,3*i-1) = 3*edofMat(:,i)-1;
    edofMat1(:,3*i) = 3*edofMat(:,i);
end

%% ========================== BOUNDARY CONDITION ======================= %%
fixednid = [1:nodeIncrement]'; % Node IDs
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
F = sparse(ndofs,1);
U = zeros(ndofs,1);
freedofs = setdiff(1:ndofs,fixeddof);

%% ===== Redefine layer-wise level set function of Geo and Path ======== %%
GeoLSgridPhi = GeoLSgridPhi_Assemble;
ScanLSgridPhi = Initialize(nelx,nely,nelz); % [nelx+2,nely+2,nelz]

for i = 1:nelz
    GeoLSgridPhi(:,:,i) = Reinitialize( GeoLSgridPhi(:,:,i),1,1,20);
    ScanLSgridPhi(:,:,i) = Reinitialize( ScanLSgridPhi(:,:,i), 1, 1,80);
end

S1 = ScanLSgridPhi(:,:,1);
q1 = 1:nelx+2;
q2 = 1:nely+2;
[Q1,Q2] = meshgrid(q1,q2);
figure(1)
contour(Q1,Q2,S1,[-10:1:10],'ShowText','on')

S2 = GeoLSgridPhi(:,:,2);
S3 = GeoLSgridPhi(:,:,10);
S4 = GeoLSgridPhi(:,:,20);
q1 = 1:nelx;
q2 = 1:nely;
[Q3,Q4] = meshgrid(q1,q2);
figure(2)
contour(Q3,Q4,S2,[0:1:10],'ShowText','on')
figure(3)
contour(Q3,Q4,S3,[0:1:10],'ShowText','on')
figure(4)
contour(Q3,Q4,S4,[0:1:10],'ShowText','on')


%% ======================== OPTIMIZATION =============================== %%
Theta = zeros(nely,nelx,nelz);   % angle of scanning path of each element
epsilon = zeros(6,nele);         % Inherent strains of each element

[Theta,~,~,~,~,~] = ThetaScanningPath_4(ScanLSgridPhi,dx,dy); 
KE = E*lk_H8(nu);
xPhys = ones(nele,1);
m=0;
for k = 1:nelz
    for j = 1:nely
        for i = 1: nelx
            m = m + 1;
            if struc( j,i,k ) == 0
                xPhys(m,1) = 0.001;
            end
        end
    end
end

iK = reshape(kron(edofMat1,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat1,ones(1,24))',24*24*nele,1);
sK = reshape(KE(:)*xPhys(:)',24*24*nele,1);

K = sparse(iK,jK,sK); K = (K+K')/2;
F_th = sparse(ndofs,1); % DEFINE THERMAL STRESS

object = [];
LLL = []; UUU=[];

for iterNum = 1:100
% update the F vector caused by scanning pattern related inherent strain 
    tic
    [epsilon,depsilon] = Inherent_strain( Theta );    % Inherent strain vector  
    F_th = sparse(ndofs,1);
    ele = 0;
    for elz = 1:nelz
        for ely = 1:nely
            for elx = 1:nelx
                ele = ele + 1;
                f_th = struc(ely,elx,elz)*Be'*C*epsilon(:,ele)*L^3;
                F_th(edofMat1(ele,:)',1) = F_th(edofMat1(ele,:)',1) + f_th;
            end
        end
    end
    
    U(freedofs,:) = K(freedofs,freedofs)\( F(freedofs,:) + F_th(freedofs,:) );
    compliance = ( F + F_th)'*U;
    object = [object compliance];
    toc
    
    [~,~,Dtheta_dphi_i_minus_1,Dtheta_dphi_i_plus_1,Dtheta_dphi_j_minus_1,Dtheta_dphi_j_plus_1] = ThetaScanningPath_4(ScanLSgridPhi,dx,dy);
    %[shapeSens] = Sensitivity2( Be,C,struc,U,edofMat1,depsilon,ScanLSgridPhi,Dtheta_dphi_i_minus_1,Dtheta_dphi_i_plus_1,Dtheta_dphi_j_minus_1,Dtheta_dphi_j_plus_1);
    [shapeSens] = Sensitivity( Be,C,struc,U,edofMat1,depsilon,Dtheta_dphi_i_minus_1,Dtheta_dphi_i_plus_1,Dtheta_dphi_j_minus_1,Dtheta_dphi_j_plus_1);
    
    if iterNum <=  50
        shapeSens = sign(shapeSens).*( abs( shapeSens).^0.3);

    elseif iterNum <= 100
        shapeSens = sign(shapeSens).*( abs( shapeSens).^0.5 );
    elseif iterNum <= 150
        shapeSens = sign(shapeSens).*( abs( shapeSens).^1.0 );
    else 
        shapeSens = sign(shapeSens).*( abs( shapeSens) );
    end
    
  
    for k = 1:nelz
        Sens  = shapeSens(:,:,k);
        Sens_median = median( median( abs(Sens) ) );
        if abs( Sens_median ) >= 0.01
            
            for j = 1:nely
                for i = 1:nelx
                    if abs( Sens(j,i)/Sens_median ) >= 5
                        shapeSens(j,i,k) = 0;
                    end
                end
            end
        end     
    end
    
    %=====================================================================%
   
    if iterNum <= 50
        [ScanLSgridPhi] = evolve(shapeSens,ScanLSgridPhi,4);
    elseif iterNum <= 150
        [ScanLSgridPhi] = evolve(shapeSens,ScanLSgridPhi,2);
    else
        [ScanLSgridPhi] = evolve(shapeSens,ScanLSgridPhi,1);
    end
  
    
    %{
    if mod( iterNum,1) == 0
        

       
        for i = 1:nelz
            S = ScanLSgridPhi(:,:,i);
            S(S>6) = -1; S(S<-6) = 1;
            for m = 1:10
                S = Reinitialize(S,1,1,80);
            end
            
            ScanLSgridPhi(:,:,i) = S;
           % ScanLSgridPhi(:,:,i) = Reinitialize( ScanLSgridPhi(:,:,i), 1, 1,80);
        end
     
        
    end
    
    %}
    

    
    
    if mod( iterNum,1) == 0
        for i = 1:5   
            S = ScanLSgridPhi(:,:,i);
            %S = Reinitialize( ScanLSgridPhi(:,:,i), 1, 1,80);
            if iterNum <=20
                
                lower = ceil( min( S(:) ) ); 
                upper = floor( max( S(:) ) );
                LLL = [LLL lower]; UUU = [UUU upper];
                const = [ lower:1:upper ];
                [perimeter] = Boundary_integration( S,  const );
                [ ~,pos] = max( perimeter);
                if mod( iterNum,600) == 0
                    S = S - const(1,pos);
                else
                    if abs( upper) >= abs( lower)
                        S = S - upper/10;
                    else
                        S = S - lower/10;
                    end
                    
                end
                S = Reinitialize( S, 1, 1,80);
                %S = S - const(1,pos);
            
            else
                S = Reinitialize( S, 1, 1,80);
            end
            
            
            
            %{
            for m = 1:10
                S( S>(upper + const(1,pos) )/2  ) = -1; S( S <( lower + const(1,pos) )/2  ) = 1;
                S = Reinitialize(S,1,1,80);
            end
            ScanLSgridPhi(:,:,i) = S;
            %}
            
            ScanLSgridPhi(:,:,i) = S;
        end
 
    end  
   
    
    

    [Theta,~,~,~,~,~] = ThetaScanningPath_4(ScanLSgridPhi,dx,dy);
    
    S1 = ScanLSgridPhi(:,:,5);
    S2 = ScanLSgridPhi(:,:,10);
    S3 = ScanLSgridPhi(:,:,20);
    S4 = ScanLSgridPhi(:,:,30);
    q1 = 1:nelx+2;
    q2 = 1:nely+2;
    [Q1,Q2] = meshgrid(q1,q2);
    figure(1)
    contour(Q1,Q2,S1,[-10:1:10],'ShowText','on')
    figure(2)
    contour(Q1,Q2,S2,[-10:1:10],'ShowText','on')
    figure(3)
    contour(Q1,Q2,S3,[-10:1:10],'ShowText','on')
    figure(4)
    contour(Q1,Q2,S4,[-10:1:10],'ShowText','on')
    
    Sens1 = shapeSens(:,:,15);
    for j = 1:nely 
        Sens( nely - j + 1,:) =  Sens1(j,:);
    end
    figure(5)
    figure(5)
    imagesc(Sens);colorbar('EastOutside');
    axis equal; axis tight; axis off;
    figure(6)
    plot( 1:iterNum, object);

    

end


