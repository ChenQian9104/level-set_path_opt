figure(5)
imagesc(Sens);colorbar('EastOutside');
axis equal; axis tight; axis on;

alpha = -pi/2:pi/180:pi/2;
x = cx(1)+0.5 + 5*cos(alpha);
y = cy(1) + 5*sin(alpha);
hold on;
plot(x,y,'b-');
fill(x,y,'w');


alpha = 0:pi/180:2*pi;
x = cx(2) + 5*cos(alpha);
y = cy(2) + 5*sin(alpha);
hold on;
plot(x,y,'b-');
fill(x,y,'w');

x = cx(3)  + 5*cos(alpha);
y = cy(3) + 5*sin(alpha);
hold on;
plot(x,y,'b-');
fill(x,y,'w');

alpha = pi/2:pi/180:3*pi/2;
x = cx(4) + 5*cos(alpha);
y = cy(4) + 5*sin(alpha);
hold on;
plot(x,y,'b-');
fill(x,y,'w');
axis equal; axis tight;