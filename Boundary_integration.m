function [perimeter] = Boundary_integration( ScanLSgridPhi,const)
nely = size(  ScanLSgridPhi,1); 
nelx = size(  ScanLSgridPhi,2); 
dim = size(  const,2); 
perimeter = zeros(dim,1);
Dirac_S = zeros( nely,nelx);
bandwidth = 1.5;

for k = 1:dim
    for ely = 2:nely-1
        
        for elx = 2:nelx-1
            phi = ScanLSgridPhi(ely,elx) - const(1,k);
            if phi<= bandwidth & phi >= -bandwidth
                Dirac_S(ely,elx) = 1/(2*bandwidth) + cos( pi*phi/bandwidth)/(2*bandwidth);
            end
        end
        
    end
    perimeter(k) = sum( Dirac_S(:) );
end





