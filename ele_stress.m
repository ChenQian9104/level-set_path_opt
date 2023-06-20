function [sp,vm, dvm,sigma] = ele_stress(C, epsilon,struc,p)
nele = size(epsilon,2); nely = size( struc,1); nelx = size( struc, 2); nelz = size( struc,3);
vm = zeros(nele,1);    % Von-Mises stress of each element
dvm = zeros(nele,1);   % [ dvm/dxx, dvm/dyy, ... , dvm/dxz];
sigma = zeros(nele,6); % Stress Components of each element [ xx, yy, zz, xy,yz,xz]

m =0;
for elz = 1:nelz
    for ely = 1:nely
        for elx = 1:nelx
            m = m + 1;
            sigma(m,:) = struc(ely,elx,elz)*( C*epsilon(:,m) )';
            vm(m,1) = ( sigma(m,1) - sigma(m,2) )^2 + ( sigma(m,2) - sigma(m,3) )^2 + ...
                ( sigma(m,3) - sigma(m,1) )^2 + 6*( sigma(m,4)^2 + sigma(m,5)^2 + sigma(m,6)^2 ); 
            vm(m,1) = sqrt( vm(m,1)/2 );
            if vm(m,1) == 0
                dvm(m,:) = 0;
            else
                dvm(m,1) = ( 2*sigma(m,1) - sigma(m,2) - sigma(m,3) )/(2*vm(m));
                dvm(m,2) = ( 2*sigma(m,2) - sigma(m,1) - sigma(m,3) )/(2*vm(m));
                dvm(m,3) = ( 2*sigma(m,3) - sigma(m,1) - sigma(m,2) )/(2*vm(m));
                dvm(m,4) = 3*sigma(m,4)/vm(m);
                dvm(m,5) = 3*sigma(m,5)/vm(m);
                dvm(m,6) = 3*sigma(m,6)/vm(m);
            end
        end
    end
end

sp = sum( vm(:).^p )^(1/p);

