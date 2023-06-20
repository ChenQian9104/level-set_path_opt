function [depx_dtheta1, depy_dtheta1, depz_dtheta1, depx_dtheta2,depy_dtheta2, depz_dtheta2] = D_strain(theta)
nely = size(theta,1); nelx = size(theta,2); nelz = size(theta,3); %epsilon = zeros(6,nelx*nely*nelz);
fx = zeros(nely,nelx,nelz); fy =zeros(nely,nelx,nelz); fz = zeros(nely,nelx,nelz); 
theta0 = zeros(nely,nelx,nelz); fxy = zeros(nely,nelx,nelz);
depx_dtheta1 = zeros(nely,nelx,nelz); depy_dtheta1 = zeros(nely,nelx,nelz);  depz_dtheta1 = zeros(nely,nelx,nelz);
depx_dtheta2 = zeros(nely,nelx,nelz); depy_dtheta2 = zeros(nely,nelx,nelz);  depz_dtheta2 = zeros(nely,nelx,nelz);

dfx_dtheta1 = zeros(nely,nelx,nelz); dfx_dtheta2 = zeros(nely,nelx,nelz);
dfy_dtheta1 = zeros(nely,nelx,nelz); dfy_dtheta2 = zeros(nely,nelx,nelz);
dfz_dtheta1 = zeros(nely,nelx,nelz); dfz_dtheta2 = zeros(nely,nelx,nelz);

dfxy_dtheta1 = zeros(nely,nelx,nelz);  dfxy_dtheta2 = zeros(nely,nelx,nelz);


dtheta0_dtheta1 = zeros(nely,nelx,nelz);
dtheta0_dtheta2 = zeros(nely,nelx,nelz);

%===== Curve fitting for inherent strain components with delta_theta =====%
delta_theta1 = [ 0 45 67 90 113 135 180]*pi/180;
strain_x1 = [ -0.0145 -0.012 -0.013 -0.0136   -0.013 -0.012 -0.0145];
strain_y1 = [ -0.0065 -0.008 -0.007 -0.007  -0.007 -0.008 -0.0065];
strain_z1 = [ 0.012 0.0126 0.0125 0.0125 0.0125 0.0126 0.012];
p1 = polyfit( delta_theta1, strain_x1, 6);
p2 = polyfit( delta_theta1, strain_y1, 6);
p3 = polyfit( delta_theta1, strain_z1, 6);

for elz = 1:nelz-1
    for ely = 1:nely
        for elx = 1:nelx
            delta_theta = abs( theta( ely, elx, elz ) - theta( ely, elx, elz+1) );
            fx( ely, elx, elz ) = p1(1)*delta_theta^6 + p1(2)*delta_theta^5 + p1(3)*delta_theta^4 + p1(4)*delta_theta^3+ p1(5)*delta_theta^2 + p1(6)*delta_theta + p1(7);
            fy( ely, elx, elz ) = p2(1)*delta_theta^6 + p2(2)*delta_theta^5 + p2(3)*delta_theta^4 + p2(4)*delta_theta^3+ p2(5)*delta_theta^2 + p2(6)*delta_theta + p2(7);
            fz( ely, elx, elz ) = p3(1)*delta_theta^6 + p3(2)*delta_theta^5 + p3(3)*delta_theta^4 + p3(4)*delta_theta^3+ p3(5)*delta_theta^2 + p3(6)*delta_theta + p3(7);  
            
            fxy(ely, elx, elz ) = sqrt( fx(ely,elx,elz)^2 + fy(ely,elx,elz)^2 );
            theta0(ely, elx, elz ) = pi + atan( fy( ely, elx, elz )/fx( ely, elx, elz ) );
            
             
            if theta( ely, elx, elz) - theta( ely, elx, elz+1) >= 0
                dfx_dtheta1( ely, elx, elz) = 6*p1(1)*delta_theta^5 + 5*p1(2)*delta_theta^4 + 4*p1(3)*delta_theta^3 + ...
                    3*p1(4)*delta_theta^2 + 2*p1(5)*delta_theta + p1(6);
                dfx_dtheta2( ely, elx, elz) = -dfx_dtheta1(ely,elx,elz);
                %dfx_dtheta2( ely, elx, elz) = -6*p1(1)*delta_theta^5 - 5*p1(2)*delta_theta^4 - 4*p1(3)*delta_theta^3 - ...
                %    3*p1(4)*delta_theta^2 - 2*p1(5)*delta_theta - p1(6);
                
                dfy_dtheta1( ely, elx, elz) = 6*p2(1)*delta_theta^5 + 5*p2(2)*delta_theta^4 + 4*p2(3)*delta_theta^3 + ...
                    3*p2(4)*delta_theta^2 + 2*p2(5)*delta_theta + p2(6);
                dfy_dtheta2( ely, elx, elz) = -dfy_dtheta1(ely,elx,elz);
                
                dfz_dtheta1( ely, elx, elz) = 6*p3(1)*delta_theta^5 + 5*p3(2)*delta_theta^4 + 4*p3(3)*delta_theta^3 + ...
                    3*p3(4)*delta_theta^2 + 2*p3(5)*delta_theta + p3(6);
                dfz_dtheta2( ely, elx, elz) = -dfz_dtheta1(ely,elx,elz);                
                        
                depz_dtheta1(ely,elx,elz) = dfz_dtheta1(ely,elx,elz);
                depz_dtheta2(ely,elx,elz) = dfz_dtheta2(ely,elx,elz);
                
            else
                
                dfx_dtheta1( ely, elx, elz) = -6*p1(1)*delta_theta^5 - 5*p1(2)*delta_theta^4 - 4*p1(3)*delta_theta^3 - ...
                    3*p1(4)*delta_theta^2 - 2*p1(5)*delta_theta - p1(6);    
                dfx_dtheta2( ely, elx, elz) = -dfx_dtheta1( ely, elx, elz);
                
                dfy_dtheta1( ely, elx, elz) = -6*p2(1)*delta_theta^5 - 5*p2(2)*delta_theta^4 - 4*p2(3)*delta_theta^3 - ...
                    3*p2(4)*delta_theta^2 - 2*p2(5)*delta_theta - p2(6);    
                dfy_dtheta2( ely, elx, elz) = -dfy_dtheta1( ely, elx, elz);  
                
                dfz_dtheta1( ely, elx, elz) = -6*p3(1)*delta_theta^5 - 5*p3(2)*delta_theta^4 - 4*p3(3)*delta_theta^3 - ...
                    3*p3(4)*delta_theta^2 - 2*p3(5)*delta_theta - p3(6);    
                dfz_dtheta2( ely, elx, elz) = -dfz_dtheta1( ely, elx, elz);                
                
            end
            
            dfxy_dtheta1( ely, elx, elz ) = ( fx(ely,elx,elz)*dfx_dtheta1(ely,elx,elz) +...
                fy(ely,elx,elz)*dfy_dtheta1(ely,elx,elz) )/fxy(ely,elx,elz);
            
            dfxy_dtheta2( ely, elx, elz ) = ( fx(ely,elx,elz)*dfx_dtheta2(ely,elx,elz) +...
                fy(ely,elx,elz)*dfy_dtheta2(ely,elx,elz) )/fxy(ely,elx,elz);
            
            dtheta0_dtheta1( ely, elx, elz ) = ( dfy_dtheta1(ely, elx, elz)*fx( ely, elx, elz ) - dfx_dtheta1(ely,elx,elz)*fy(ely,elx,elz) )/...
                ( fx( ely, elx, elz )^2 + fy( ely, elx, elz )^2 );
            
            dtheta0_dtheta2( ely, elx, elz ) = ( dfy_dtheta2(ely, elx, elz)*fx( ely, elx, elz ) - dfx_dtheta2(ely,elx,elz)*fy(ely,elx,elz) )/...
                ( fx( ely, elx, elz )^2 + fy( ely, elx, elz )^2 );            
            
            
            if cos( theta0(ely, elx, elz ) + theta(ely,elx,elz) )>= 0
              
                depx_dtheta1(ely,elx,elz) = -dfxy_dtheta1(ely,elx,elz)*cos( theta0(ely, elx, elz ) + theta(ely,elx,elz) ) + ...
                    fxy(ely,elx,elz)*sin( theta0(ely, elx, elz ) + theta(ely,elx,elz) )*( 1+dtheta0_dtheta1(ely,elx,elz) );
                
                depx_dtheta2(ely,elx,elz) = -dfxy_dtheta2(ely,elx,elz)*cos( theta0(ely, elx, elz ) + theta(ely,elx,elz) ) + ...
                    fxy(ely,elx,elz)*sin( theta0(ely, elx, elz ) + theta(ely,elx,elz) )*( dtheta0_dtheta2(ely,elx,elz) );
                
            else
                depx_dtheta1(ely,elx,elz) = dfxy_dtheta1(ely,elx,elz)*cos( theta0(ely, elx, elz ) + theta(ely,elx,elz) ) - ...
                    fxy(ely,elx,elz)*sin( theta0(ely, elx, elz ) + theta(ely,elx,elz) )*( 1+dtheta0_dtheta1(ely,elx,elz) );
                
                depx_dtheta2(ely,elx,elz) = dfxy_dtheta2(ely,elx,elz)*cos( theta0(ely, elx, elz ) + theta(ely,elx,elz) ) - ...
                    fxy(ely,elx,elz)*sin( theta0(ely, elx, elz ) + theta(ely,elx,elz) )*( dtheta0_dtheta2(ely,elx,elz) );                
            end
            
            
            if sin( theta0(ely,elx,elz) + theta(ely,elx,elz) )>=0
                depy_dtheta1(ely,elx,elz) = -dfxy_dtheta1(ely,elx,elz)*sin( theta0(ely,elx,elz) + theta(ely,elx,elz) ) - ...
                    fxy(ely,elx,elz)*cos( theta0(ely,elx,elz) + theta(ely,elx,elz) )*(1+ dtheta0_dtheta1(ely,elx,elz) );
                
                depy_dtheta2(ely,elx,elz) = -dfxy_dtheta2(ely,elx,elz)*sin( theta0(ely,elx,elz) + theta(ely,elx,elz) ) - ...
                    fxy(ely,elx,elz)*cos( theta0(ely,elx,elz) + theta(ely,elx,elz) )*dtheta0_dtheta2(ely,elx,elz);
               
            else
                depy_dtheta1(ely,elx,elz) = dfxy_dtheta1(ely,elx,elz)*sin( theta0(ely,elx,elz) + theta(ely,elx,elz) ) + ...
                    fxy( ely,elx,elz)*cos( theta0(ely,elx,elz) + theta(ely,elx,elz) )*(1+ dtheta0_dtheta1(ely,elx,elz) ); 
                depy_dtheta2(ely,elx,elz) = dfxy_dtheta2(ely,elx,elz)*sin( theta0(ely,elx,elz) + theta(ely,elx,elz) ) + ...
                    fxy(ely,elx,elz)*cos( theta0(ely,elx,elz) + theta(ely,elx,elz) )*dtheta0_dtheta2(ely,elx,elz);
            end

        end
    end
end

for ely = 1:nely
    for elx = 1:nelx
        fx(ely,elx,nelz) = -0.02;
        fy(ely,elx,nelz) = -0.002;
        fz(ely,elx,nelz) = 0.0155;
        fxy(ely, elx, nelz ) = sqrt( fx(ely,elx,nelz)^2 + fy(ely,elx,nelz)^2 );
        theta0(ely,elx,nelz) = pi + atan( fy( ely, elx, nelz )/fx( ely, elx, nelz ) );
        
        if cos( theta0(ely,elx,nelz) + theta(ely,elx, nelz) )>=0
            depx_dtheta1(ely,elx,nelz) = fxy(ely,elx,nelz)*sin( theta0(ely,elx,nelz) + theta(ely,elx, nelz) );
        else
            depx_dtheta1(ely,elx,nelz) = -fxy(ely,elx,nelz)*sin( theta0(ely,elx,nelz) + theta(ely,elx, nelz) );
        end
        
        
        if sin( theta0(ely,elx,nelz) + theta(ely,elx, nelz) )>=0
            depy_dtheta1(ely,elx,nelz) = -fxy(ely,elx,nelz)*cos( theta0(ely,elx,nelz) + theta(ely,elx, nelz) );
        else
            depy_dtheta1(ely,elx,nelz) =  fxy(ely,elx,nelz)*cos( theta0(ely,elx,nelz) + theta(ely,elx, nelz) );
        end        
        
    end

end





