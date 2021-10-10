%% p_cyl.m
% arguments: a height, a radius, and a number of cell-elements
% Generates a point cloud for occlusion calculation of a cylinder of height
% h, radius r, and n points vertically, with normals pointing in the direction of 
% the cross product of the i and j vectors.
% Output quadmesh contains x-data, y-data, z-data, cell normals, and area
  
function p_cloud = p_cyl(r_h,n,varargin)

if nargin==2
    origin=[0,0,0];
elseif nargin==3
    origin=varargin{1};
elseif nargin==1
    error('missing input ''n'', number of nodes in z direction')
elseif nargin==4
    error('missing input ''z'', origin z location')
elseif nargin==5
    origin=[varargin{1},varargin{2},varargin{3}];
else
    error('unknown extra inputs')
end

x  = 1;
y  = 2;
z  = 3;
h  = 1;
r  = h*r_h;
ds = h/n;
c  = 2*pi*r;
nr = ceil(c/ds);
hp = linspace(0,2*pi,nr);
xp = r*cos(hp)';
yp = r*sin(hp)';
zp = linspace(-h/2,h/2,n);
xyr= ones(size(zp));
zr = ones(size(xp));
xp = xp.*xyr;
yp = yp.*xyr;
zp = zp.*zr;
cc(:,x) = xp(:);
cc(:,y) = yp(:);
cc(:,z) = zp(:);

cc=reshape(cc,[],3);

p_cloud=cc+origin;

end