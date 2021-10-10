close all; clear all; clc;
n=40;
[X,Y,Z] = sphere(n);
x = [X(:)];
y = [Y(:)];
z = [Z(:)];

% Matrix of x,y,z coordinates of points
p = [x y z];

% Viewpoint coordinates
C = [1000 1000 0];

param = 1;

% Array of indices of visible points
index = HPR(p,C,param);

scatter3(x(index),y(index),z(index))
xlabel('x')
ylabel('y')
zlabel('z')
