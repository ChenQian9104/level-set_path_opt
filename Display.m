vm = struc; L = 1;
nodeIncrement = (nelx+1)*(nely+1);
nnode = (nelx+1)*(nely+1)*nelz;

figure(2)
nodeCor = zeros(nnode, 3);
m = 0;
for j = 1:nely+1
    for i = 1:nelx + 1
        m = m + 1;
        nodeCor(m,:) = [ (i-1)*L, (j-1)*L, 0];
    end
end

for k = 1:nelz
    for i = 1:nodeIncrement
        nodeCor( i + k*nodeIncrement, : ) = nodeCor( i + (k-1)*nodeIncrement,:) + [ 0 0 L];
    end
end

m=0;
for k = 1:nelz
    for j = 1:nely
        for i = 1:nelx
            m = m + 1;
            vert = [ nodeCor( edofMat(m,1),: ); nodeCor( edofMat(m,2),: );...
                nodeCor(edofMat(m,3),:); nodeCor(edofMat(m,4),:); ...
                nodeCor(edofMat(m,5),:); nodeCor(edofMat(m,6),:); ...
                nodeCor(edofMat(m,7),:); nodeCor(edofMat(m,8),:)];
            face = [ 1 2 3 4; 3 2 6 7; 7 6 5 8; 8 5 1 4; 3 4 8 7; 1 2 6 5];
            if struc( j,i,k) == 1 
%                 patch('Faces',face,'Vertices',vert,'CData',[ vm(m), vm(m),vm(m),...
%                     vm(m),vm(m),vm(m)],'FaceColor','flat');
                patch('Faces',face,'Vertices',vert,'CData',[ vm(j,i,k), vm(j,i,k),vm(j,i,k),...
                    vm(j,i,k),vm(j,i,k),vm(j,i,k)],'FaceColor','flat');                
                axis equal; axis tight;box on;view([-45 45]);
            end
        end
    end
end




