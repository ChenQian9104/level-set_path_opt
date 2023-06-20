figure(2)
contour(Q1,Q2,S2,[-9:1:8],'ShowText','off');axis equal ; axis tight

alpha = -pi/2:pi/180:pi/2;
x = cx(1) + 1.5 + 5*cos(alpha);
y = cy(1)+1 + 5*sin(alpha);
hold on;
plot(x,y,'b-');
fill(x,y,'w');


alpha = 0:pi/180:2*pi;
x = cx(2) + 1.5 + 5*cos(alpha);
y = cy(2)+1 + 5*sin(alpha);
hold on;
plot(x,y,'b-');
fill(x,y,'w');

x = cx(3) + 1.5 + 5*cos(alpha);
y = cy(3)+1 + 5*sin(alpha);
hold on;
plot(x,y,'b-');
fill(x,y,'w');

alpha = pi/2:pi/180:3*pi/2;
x = cx(4)+1.5 + 5*cos(alpha);
y = cy(4)+1 + 5*sin(alpha);
hold on;
plot(x,y,'b-');
fill(x,y,'w');
axis equal; axis tight;