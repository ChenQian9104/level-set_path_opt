 function [ScanLSgridPhi] = evolve2(v,ScanLSgridPhi,stepLength)
 nely = size(v,1); nelx = size(v,2); nelz = size(v,3);
 % Extend sensitivites using a zero border
 vFull = zeros( nely + 2, nelx + 2, nelz); vFull(2:end-1,2:end-1,:) = v;
 % Choose time step for evolution based onCFL value
 %dt = 0.1/max(abs(v(:)));
 
 % Update level set function layer-by-layer
 for elz = 1:nelz
     lsf = ScanLSgridPhi(:,:,elz);
     v1 = v(:,:,elz);
     if max( abs( v1(:) ) ) > 0
         dt = 0.1/max( abs( v1(:) ) );
         for i = 1:10*stepLength
             dpx = circshift(lsf,[0,-1])-lsf;
             dmx = lsf - circshift(lsf,[0,1]);
             dpy = circshift(lsf,[-1,0]) - lsf;
             dmy = lsf - circshift(lsf,[1,0]);
             lsf = lsf - dt * min(vFull(:,:,elz) ,0).* ...
                 sqrt( min(dmx,0).^2+max(dpx,0).^2+min(dmy,0).^2+max(dpy,0).^2 ) ...
                 - dt * max(vFull(:,:,elz),0) .*...
                 sqrt( max(dmx,0).^2+min(dpx,0).^2+max(dmy,0).^2+min(dpy,0).^2 ); 
         end
         ScanLSgridPhi(:,:,elz) = lsf;
     end
    
   
 end
 



