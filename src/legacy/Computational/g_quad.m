%% g_quad.m
% arguments: a corner (xyz), two xyz vectors, and a number of cell-elements
% Generates a quadrilateral mesh, with normals pointing in the direction of 
% the cross product of the i and j vectors.
% Output quadmesh contains x-data, y-data, z-data, cell normals, and area
  
function quadmesh = g_quad(origin,ivec,jvec,varargin)
x=1;
y=2;
z=3;
if nargin==4
    ni=varargin{1};
    nj=max(1,round(norm(ivec)/norm(jvec)*ni));
elseif nargin==5
    ni=varargin{1};
    nj=varargin{2};
elseif nargin==3
    error('Missing number of cells!')
else
    error('Too many inputs!')
end

xi=linspace(0,ivec(x),ni+1);
yi=linspace(0,ivec(y),ni+1);
zi=linspace(0,ivec(z),ni+1);
xj=linspace(0,jvec(x),nj+1);
yj=linspace(0,jvec(y),nj+1);
zj=linspace(0,jvec(z),nj+1);

xn=zeros(ni+1,nj+1)+NaN;
yn=xn;
zn=xn;
for jj=1:nj+1
    xn(:,jj)=xi+xj(jj);
    yn(:,jj)=yi+yj(jj);
    zn(:,jj)=zi+zj(jj);
end

cc(:,:,x)=(xn(1:end-1,1:end-1)+xn(2:end  ,1:end-1)+...
           xn(1:end-1,2:end  )+xn(2:end  ,2:end  ))./4;
cc(:,:,y)=(yn(1:end-1,1:end-1)+yn(2:end  ,1:end-1)+...
           yn(1:end-1,2:end  )+yn(2:end  ,2:end  ))./4;
cc(:,:,z)=(zn(1:end-1,1:end-1)+zn(2:end  ,1:end-1)+...
           zn(1:end-1,2:end  )+zn(2:end  ,2:end  ))./4;

uc=zeros(ni,nj,3)+NaN;
uc(:,:,x)=xn(2:end  ,2:end  )-xn(1:end-1,1:end-1);
uc(:,:,y)=yn(2:end  ,2:end  )-yn(1:end-1,1:end-1);
uc(:,:,z)=zn(2:end  ,2:end  )-zn(1:end-1,1:end-1);

vc=zeros(ni,nj,3)+NaN;
vc(:,:,x)=xn(1:end-1,2:end  )-xn(2:end  ,1:end-1);
vc(:,:,y)=yn(1:end-1,2:end  )-yn(2:end  ,1:end-1);
vc(:,:,z)=zn(1:end-1,2:end  )-zn(2:end  ,1:end-1);

nc=cross(uc,vc,3)./2;
Ac=vecnorm(nc,2,3);
nc=nc./Ac;

cc=reshape(cc,[],3);
nc=reshape(nc,[],3);
Ac=reshape(Ac,[],1);

quadmesh.xyz=cc+origin;
quadmesh.n=nc;
quadmesh.A=Ac;
end
