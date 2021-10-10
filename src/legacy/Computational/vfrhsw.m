function vfr=vfrhsw(d_l,h_l)

n=1e3;

% Here we will handle saving/loading our view factors to a pair of external
% binary files. vfr_n.bin holds the record of complete simulations: the
% first is nn, a uint16 holding the number of records, and followed by a
% 3xnn double holding 1: the D/L of the simulation completed, 2: the H/L of
% the simulation, and 3: the number of bundles simulated.

if exist('vfr_n.bin','file')==2 % check for existence of records file
    fid=fopen('vfr_n.bin','r');
    nn =fread(fid,1,'uint16'); % read the number of records
    dhn=fread(fid,[3,nn],'double'); % read the records information
    fclose(fid); % close the records file
    % next, find the index where the D/L and H/L are equal to a record
    % currently on file. Tolerance set to 'tol'
    tol = 1e-3;
    ii = find(abs(dhn(1,:)-d_l)<tol & abs(dhn(2,:)-h_l)<tol,1);
    % if this simulation has been done before, 
    if(~isempty(ii))
        if n>dhn(3,ii) % and the requested accuracy is greater than record,
            rbool=false; % then do a new simulation
        else % otherwise
            rbool=true; % read the saved view factors
        end
    else % if this combination of parameters doesn't exist, prepare to i/o
        nn=nn+1; % number of records is end + 1
        ii=nn;   % and location of writing is at the end.
        rbool=false; % and write when done with simulatons.
    end
    nbool=false; % if the file exists, do not need to make a new file.
else
    nbool=true; % if the file does not exist, write a new file, and 
    rbool=false; % write to it when done.
    nn=1; % the number of records is now 1, and 
    ii=1; % the index of writing is also 1.
end

if rbool % if we're reading from a file, 
    fid=fopen('vfr_r.bin','r'); % open the vfr database, 
    fseek(fid,8*16*(ii-1),-1); % scan to the ii'th index,
    vfr=fread(fid,16,'double'); % and read fi-j, with i and j from 1-4,
    fclose(fid); % and then close the file. 
else % otherwise, perform monte carlo or analytic simulations. 
    if h_l==0&&d_l==0 % return planar values
        vfr=[0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0];
    elseif h_l==0 % return hex/sphere mc simulation values
        vfr=vfrhs_mc(d_l,n);
    elseif d_l==0 % return hex/wall analytic values
        vfr=vfrhw(h_l);
    else % return hex/wall/sphere mc simulation values
        vfr=vfrhsw_mc(d_l,h_l,n);
    end
    
    % then, write record to file.
    if nbool % if no records previously existed, open a new file to write
        fid=fopen('vfr_n.bin','w');
    else % otherwise, open vfr_n.bin to write into
        fid=fopen('vfr_n.bin','r+');
    end
    fwrite(fid,nn,'uint16'); % write new number of records,
    fseek(fid,8*3*(ii-1),0); % seek to record to write to,
    fwrite(fid,[d_l,h_l,n],'double'); % and write new d/l,h/l, and n
    fclose(fid); % and close
    
    % then, write vfr to file
    if nbool % if no vfr previously existed, open a new file to write
        fid=fopen('vfr_r.bin','w');
    else % otherwise, open vfr_r.bin to write into
        fid=fopen('vfr_r.bin','r+');
    end
    fseek(fid,8*16*(ii-1),-1); % seek to record to write to,
    fwrite(fid,vfr,'double'); % and write new Fi-j, for i&j 1-4
    fclose(fid); % and close
end

% and then return values for F11, F12, F13, F14, F21, ... , F44

end


% VFRHS_MC view factor relationships, hex/sphere, monte carlo simulation. 
% returns [fij] where i = [1-4], j = [1-4]
function vfr = vfrhs_mc(d_l,n)

l = 1;        % length scale of unit cell
d = d_l*l;    % diameter of central sphere
r2= d*d/4;    % radius squared
s = l/sqrt(3);% length of one side of hex
a1 = 3/2*s*l;% area of substrate
a3 = a1;    % apparent area
a2 = 4*pi*r2; % area of sphere

% monte carlo parameters

act = true(1,n); % is bundle ii active?
f11 = false(1,n);% is bundle ii incident on surface 1?
f12 = f11;       % is bundle ii incident on surface 2?
f13 = f12;       % is bundle ii incident on surface 3?
nrm = 10;          % maximum number of bounces before terminate
mcx = rand(1,n); % mc parameter for x location
mcy = rand(1,n); % mc parameter for y location
% mcz = zeroes(1,n); % not mc parameter: all at z = -d/2
mcp = rand(1,n); % mc parameter for psi direction
mcq = rand(1,n); % mc parameter for theta direction

% dbg = true;
% dbg2= dbg;

if exist('dbg','var')&&dbg
    figure
    scatter3((-s/2+s*mcx).*sqrt(mcy),l/2*sqrt(mcy),-d/2*ones(1,n))
    hold on
    plot3([s s/2 -s/2 -s -s/2 s/2 s],...
          [0 l/2 l/2 0 -l/2 -l/2 0],...
     -d/2*[1 1 1 1 1 1 1])
    plot3([s s/2 -s/2 -s -s/2 s/2 s],...
          [0 l/2 l/2 0 -l/2 -l/2 0],...
      d/2*[1 1 1 1 1 1 1])
end % if we're debugging, set up visualizer

% loop through n, determine incident surfaces for each bundle
for ii = 1:n
    ox = (-s/2+s*mcx(ii))*sqrt(mcy(ii)); % origin from x=[-s/2 s/2]*y/ymax
    oy = l/2*sqrt(mcy(ii));         % origin from y=[0 l/2]
    oz = -d/2;                      % origin at z=-d/2
    np = 2*pi*mcp(ii);              % phi from [0 2pi]
    nq = asin(sqrt(mcq(ii)));       % the from [0 pi/2]
    nx = sin(nq)*cos(np);           % pointing x
    ny = sin(nq)*sin(np);           % pointing y
    nz = cos(nq);                   % pointing z
    
    nr = 0; % reset reflection count
    
    while nr<nrm&&act(ii) % while bundle is active, until max reflection
        
        % first we will check for intersection with sphere:
        loc  = -nx*ox-ny*oy-nz*oz; % distance from o to closest approach
        lca2 = (ox^2+oy^2+oz^2)-loc^2; % closest approach, squared
        if lca2<r2 % if closest approach is less than the radius
            f12(ii) = true; % then it intersects 2
            act(ii) = false; % and terminate the bundle
            
            if exist('dbg2','var')&&dbg2
                lct2 = r2-lca2;
                s2 = loc-sqrt(lct2);
                tx = ox+nx*s2;
                ty = oy+ny*s2;
                tz = oz+nz*s2;
                plot3([ox tx],[oy,ty],[oz,tz])
            end % if we're debugging, visualize sphere intersect
            
        else
        
        % check for intersection with upper surface
        tz = d/2; % if so, we'll intersect at d/2
        s3 = (tz-oz)/nz; % scale is dz/nz;
        tx = ox+nx*s3;   % intersection in x/y is o+n*s
        ty = oy+ny*s3;
        if abs(ty)<l/2&&abs(tx)<(1-abs(ty)/(l))*s % if it's inside the hex
            f13(ii) = true; % then it intersects 3
            act(ii) = false; % and terminate
            
            if exist('dbg3','var')&&dbg3
                plot3([ox tx],[oy,ty],[oz,tz])
            end % if we're debugging, visualize hex intersect
            
        else % otherwise, it hits the wall:
            % time to determine which of the six walls we hit first:
            % we'll calculate the intersection of our ray and a circle of 
            % radius s at the origin origin. The scale from our origin to 
            % the closest approach to the bounding circle is given by:
            lx   = -ox;
            ly   = -oy;
            l2   = lx^2+ly^2; % distance from ray origin to circle center
            nxy2 = sqrt(nx^2+ny^2);
            nxx  = nx/nxy2;   % scale it for a unit normal in 2d
            nyy  = ny/nxy2;
            loc  = lx*nxx+ly*nyy; % distance from origin to closest approach
            lca2 = l2-loc^2;      % closest approach squared
            lct2 = s^2-lca2;      % distance from closest approach to circle intersection
            s4xy = loc+sqrt(lct2);% scale from origin to interception
            
            % which means the location of intersection with the circle is
            tx   = ox+s4xy*nxx;
            ty   = oy+s4xy*nyy;
            
            if exist('dbg4c','var')&&dbg4c
                tz   = oz+s4xy/nxy2*nz;
                plot3([ox tx],[oy,ty],[oz,tz])
            end % if we're debugging, visualize circle intersect
            
            % so the angle from the x axis to our circle intersection is
            qi   = sign(ty)*((tx>0)*pi+sign(tx)*atan(abs(ty/tx)));
            
            % all that to find the wall we're intersecting without guessing
            % and checking every individual wall and selecting the shortest
            % non-negative distance usw usw. We use this angle to determine
            % the wall normal, the wall intersection, and then do the shift
            % of the origin for the periodic boundary condition
            qw   = floor(qi/(pi/3))*pi/3+pi/6;
            wx   = -cos(qw);
            wy   = -sin(qw);
            s4   = (l/2-ox*wx-oy*wy)/(nx*wx+ny*wy);
            
            tx   = ox+s4*nx;
            ty   = oy+s4*ny;
            tz   = oz+s4*nz;
            
            if exist('dbg4','var')&&dbg4
                plot3([ox tx],[oy,ty],[oz,tz])
            end % if we're debugging, visualize wall intersect
            
            ox   = tx-l*wx; % periodic boundary = shift to opposite wall
            oy   = ty-l*wy;
            oz   = tz;
            
            nr   = nr+1; % if we have to period-icize, ++ reflect counter
        end
        
        end

    end
end

if exist('dbg','var')&&dbg
    axis([-s s -s s -s s])
    axis square
end

na = sum(act);
% disp([num2str(100*na/n) '% rays unterminated!'])
n = n - na;

f11 = 0;
f12 = sum(f12)/n;
f13 = sum(f13)/n;

f21 = f12*a1/a2;
f23 = f21;
f22 = 1-f21-f23;

f31 = f13*a1/a3;
f33 = 0;
f32 = 1-f31-f33;

vfr = [f11 f12 f13 0 f21 f22 f23 0 f31 f32 f33 0 0 0 0 0];

end

% VFRHW_MC view factor relationships, hex/walls, analytic relationships
% returns [fij] where i = [1-4], j = [1-4]
function vfr = vfrhw(h_l)

l = 1;
h = h_l*l;
s = l/sqrt(3);

a1 = 3/2*l*s;
a3 = a1;
a4 = 6*s*h;

f4a = C11(s,h,l);

f4b = C16(s,s,h,2*pi/3);

q   = pi/3;
fbb = C16(2*s,2*s,h,q);
fbs = C16(  s,2*s,h,q);
fsb = fbs*2/1;
fss = C16(  s,  s,h,q);
fbp = fbb-fbs;
fsp = fsb-fss;
fpb = fbp*2/1;
fps = fsp;
f4c = fpb-fps;

f44 = f4a+2*f4b+2*f4c;
f41 = (1-f44)/2;
f43 = f41;

f14 = f41*a4/a1;
f13 = 1-f14;

f31 = f13*a1/a3;
f34 = f43*a4/a3;

vfr = [0 0 f13 f14 0 0 0 0 f31 0 0 f34 f41 0 f43 f44];

end

% VFRHSW_MC view factor relationships, hex/sphere/walls
function vfr = vfrhsw_mc(d_l,h_l,n)

l = 1;        % length scale of unit cell
d = d_l*l;    % diameter of central sphere
h = h_l*l;    % height of adiabatic walls
r2= d*d/4;    % radius squared
s = l/sqrt(3);% length of one side of hex
a1 = 3/2*s*l; % area of substrate
a2 = 4*pi*r2; % area of sphere
a3 = a1;      % apparent area
a4 = 6*s*h;   % area of adbat walls

% monte carlo parameters

act = true(1,n); % is bundle ii active?
f11 = false(1,n);% is bundle ii incident on surface 1?
f12 = f11;       % is bundle ii incident on surface 2?
f13 = f12;       % is bundle ii incident on surface 3?
f14 = f13;       % is bundle ii incident on surface 4?
nrm = 100;       % maximum number of bounces before terminate
mcx = rand(1,n); % mc parameter for x location
mcy = rand(1,n); % mc parameter for y location
% mcz = zeroes(1,n); % not mc parameter: all at z = -d/2
mcp = rand(1,n); % mc parameter for psi direction
mcq = rand(1,n); % mc parameter for theta direction

% dbg = true;
% dbg2= dbg;

if exist('dbg','var')&&dbg
    figure
    scatter3((-s/2+s*mcx).*sqrt(mcy),l/2*sqrt(mcy),-d/2*ones(1,n))
    hold on
    plot3([s s/2 -s/2 -s -s/2 s/2 s],...
          [0 l/2 l/2 0 -l/2 -l/2 0],...
     -d/2*[1 1 1 1 1 1 1])
    plot3([s s/2 -s/2 -s -s/2 s/2 s],...
          [0 l/2 l/2 0 -l/2 -l/2 0],...
      d/2*[1 1 1 1 1 1 1])
end % if we're debugging, set up visualizer

% loop through n, determine incident surfaces for each bundle
for ii = 1:n
    ox = (-s/2+s*mcx(ii))*sqrt(mcy(ii)); % origin from x=[-s/2 s/2]*y/ymax
    oy = l/2*sqrt(mcy(ii));         % origin from y=[0 l/2]
    oz = -d/2;                      % origin at z=-d/2
    np = 2*pi*mcp(ii);              % phi from [0 2pi]
    nq = asin(sqrt(mcq(ii)));       % the from [0 pi/2]
    nx = sin(nq)*cos(np);           % pointing x
    ny = sin(nq)*sin(np);           % pointing y
    nz = cos(nq);                   % pointing z
    
    nr = 0; % reset reflection count
    
    while nr<nrm&&act(ii) % while bundle is active, until max reflection
        
        % first we will check for intersection with sphere:
        loc  = -nx*ox-ny*oy-nz*oz; % distance from o to closest approach
        lca2 = (ox^2+oy^2+oz^2)-loc^2; % closest approach, squared
        if lca2<r2 % if closest approach is less than the radius
            f12(ii) = true; % then it intersects 2
            act(ii) = false; % and terminate the bundle
            
            if exist('dbg2','var')&&dbg2
                lct2 = r2-lca2;
                s2 = loc-sqrt(lct2);
                tx = ox+nx*s2;
                ty = oy+ny*s2;
                tz = oz+nz*s2;
                plot3([ox tx],[oy,ty],[oz,tz])
            end % if we're debugging, visualize sphere intersect
            
        else
        
        % check for intersection with upper surface
        tz = -d/2+max([d,h]); % if so we'll intersect at top of sphere/wall
        s3 = (tz-oz)/nz; % scale is dz/nz;
        tx = ox+nx*s3;   % intersection in x/y is o+n*s
        ty = oy+ny*s3;
        if abs(ty)<l/2&&abs(tx)<(1-abs(ty)/(l))*s % if it's inside the hex
            f13(ii) = true; % then it intersects 3
            act(ii) = false; % and terminate
            
            if exist('dbg3','var')&&dbg3
                plot3([ox tx],[oy,ty],[oz,tz])
            end % if we're debugging, visualize hex intersect
            
        else % otherwise, it hits the wall:
            % time to determine which of the six walls we hit first:
            % we'll calculate the intersection of our ray and a circle of 
            % radius s at the origin origin. The scale from our origin to 
            % the closest approach to the bounding circle is given by:
            lx   = -ox;
            ly   = -oy;
            l2   = lx^2+ly^2;
            nxy2 = sqrt(nx^2+ny^2);
            nxx  = nx/nxy2;
            nyy  = ny/nxy2;
            loc  = lx*nxx+ly*nyy;
            lca2 = l2-loc^2;
            lct2 = s^2-lca2;
            s4xy = loc+sqrt(lct2);
            
            % which means the location of intersection with the circle is
            tx   = ox+s4xy*nxx;
            ty   = oy+s4xy*nyy;
            
            if exist('dbg4c','var')&&dbg4c
                tz   = oz+s4xy/nxy2*nz;
                plot3([ox tx],[oy,ty],[oz,tz])
            end % if we're debugging, visualize circle intersect
            
            % so the angle from the x axis to our circle intersection is
            qi   = sign(ty)*((tx>0)*pi+sign(tx)*atan(abs(ty/tx)));
            
            % all that to find the wall we're intersecting without guessing
            % and checking every individual wall and selecting the shortest
            % non-negative distance usw usw. We use this angle to determine
            % the wall normal, the wall intersection, and then do the shift
            % of the origin for the periodic boundary condition
            qw   = floor(qi/(pi/3))*pi/3+pi/6;
            wx   = -cos(qw);
            wy   = -sin(qw);
            s4   = (l/2-ox*wx-oy*wy)/(nx*wx+ny*wy);
            
            tx   = ox+s4*nx;
            ty   = oy+s4*ny;
            tz   = oz+s4*nz;
            
            if exist('dbg4','var')&&dbg4
                plot3([ox tx],[oy,ty],[oz,tz])
            end % if we're debugging, visualize wall intersect
            
            if tz<-d/2+h
                f14(ii)=true;
                act(ii)=false;
            else
                ox   = tx-l*wx; % periodic boundary = shift to opposite wall
                oy   = ty-l*wy;
                oz   = tz;
            
                nr   = nr+1; % ++ reflect counter
            end
        end
        
        end

    end
end

if exist('dbg','var')&&dbg
    axis([-s s -s s -s s])
    axis square
end

na = sum(act);
% disp([num2str(100*na/n) '% rays unterminated!'])
n = n - na;

f12 = sum(f12)/n;
f13 = sum(f13)/n;
f14 = sum(f14)/n;

f21 = f12*a1/a2;
f31 = f13*a1/a3;
f41 = f14*a1/a4;

% for remaining view factors, if walls are higher than the sphere, you can
% utilize analytic view factors to get the missing view factors:
if h>=d
    f22 = 0;
    f24 = 3*(C121(l/2,s,h)+C121(l/2,s,d));
    f23 = 1-f21-f22-f24;
    
    f32 = f23*a2/a3;
    f34 = 1-f31-f32;
    
    f42 = f24*a2/a4;
    f43 = f34*a3/a4;
    f44 = 1-f41-f42-f43;
    
% otherwise, we'll need to do a second monte carlo simulation
else
    
act = true(1,n); % is bundle ii active?
f41c= false(1,n);% is bundle ii incident on surface 1?
f42 = f41;       % is bundle ii incident on surface 2?
f43 = f42;       % is bundle ii incident on surface 3?
f44 = f43;       % is bundle ii incident on surface 4?
nrm = 100;       % maximum number of bounces before terminate
mcx = rand(1,n); % mc parameter for x location
% mcy = zeroes(1,n); % not mc parameter: all at y = -L/2
mcz = rand(1,n); % mc parameter for z location
mcp = rand(1,n); % mc parameter for psi direction
mcq = rand(1,n); % mc parameter for theta direction

% loop through n, determine incident surfaces for each bundle
% dbg2 = true;
for ii = 1:n
    ox = -s/2+s*mcx(ii);            % origin from x=[-s/2 s/2]
    oy = -l/2;                      % origin at y=-l/2
    oz = -d/2+h*mcz(ii);            % origin at z=[-d/2 (-d/2+h)]
    np = 2*pi*mcp(ii);              % phi from [0 2pi]
    nq = asin(sqrt(mcq(ii)));       % the from [0 pi/2]
    nx = sin(nq)*cos(np);           % pointing x
    ny = cos(nq);                   % pointing y
    nz = sin(nq)*sin(np);           % pointing z
    
    nr = 0; % reset reflection count
    
    while nr<nrm&&act(ii) % while bundle is active, until max reflection
        
        % first we will check for intersection with sphere:
        loc  = -nx*ox-ny*oy-nz*oz; % distance from o to closest approach
        lca2 = (ox^2+oy^2+oz^2)-loc^2; % closest approach, squared
        if lca2<r2 % if closest approach is less than the radius
            f42(ii) = true; % then it intersects 2
            act(ii) = false; % and terminate the bundle
            
            if exist('dbg2','var')&&dbg2
                lct2 = r2-lca2;
                s2 = loc-sqrt(lct2);
                tx = ox+nx*s2;
                ty = oy+ny*s2;
                tz = oz+nz*s2;
                plot3([ox tx],[oy,ty],[oz,tz])
            end % if we're debugging, visualize sphere intersect
            
        else
        
        % check for intersection with upper/lower surface
        tz = -d/2+(nz>0)*max([d,h]); % if so we'll intersect 
        s3 = (tz-oz)/nz; % scale is dz/nz;
        tx = ox+nx*s3;   % intersection in x/y is o+n*s
        ty = oy+ny*s3;
        if abs(ty)<l/2&&abs(tx)<(1-abs(ty)/(l))*s % if it's inside the hex
            if nz>0
                f43(ii) = true; % then it intersects 3
            else
                f41c(ii) = true; % or 1, if it's pointing down
            end
            act(ii) = false; % and terminate
            
            if exist('dbg3','var')&&dbg3
                plot3([ox tx],[oy,ty],[oz,tz])
            end % if we're debugging, visualize hex intersect
            
        else % otherwise, it hits the wall:
            % time to determine which of the six walls we hit first:
            % we'll calculate the intersection of our ray and a circle of 
            % radius s at the origin origin. The scale from our origin to 
            % the closest approach to the bounding circle is given by:
            lx   = -ox;
            ly   = -oy;
            l2   = lx^2+ly^2;
            nxy2 = sqrt(nx^2+ny^2);
            nxx  = nx/nxy2;
            nyy  = ny/nxy2;
            loc  = lx*nxx+ly*nyy;
            lca2 = l2-loc^2;
            lct2 = s^2-lca2;
            s4xy = loc+sqrt(lct2);
            
            % which means the location of intersection with the circle is
            tx   = ox+s4xy*nxx;
            ty   = oy+s4xy*nyy;
            
            if exist('dbg4c','var')&&dbg4c
                tz   = oz+s4xy/nxy2*nz;
                plot3([ox tx],[oy,ty],[oz,tz])
            end % if we're debugging, visualize circle intersect
            
            % so the angle from the x axis to our circle intersection is
            qi   = sign(ty)*((tx>0)*pi+sign(tx)*atan(abs(ty/tx)));
            
            % all that to find the wall we're intersecting without guessing
            % and checking every individual wall and selecting the shortest
            % non-negative distance usw usw. We use this angle to determine
            % the wall normal, the wall intersection, and then do the shift
            % of the origin for the periodic boundary condition
            qw   = floor(qi/(pi/3))*pi/3+pi/6;
            wx   = -cos(qw);
            wy   = -sin(qw);
            s4   = (l/2-ox*wx-oy*wy)/(nx*wx+ny*wy);
            
            tx   = ox+s4*nx;
            ty   = oy+s4*ny;
            tz   = oz+s4*nz;
            
            if exist('dbg4','var')&&dbg4
                plot3([ox tx],[oy,ty],[oz,tz])
            end % if we're debugging, visualize wall intersect
            
            if tz<-d/2+h % if we're below the adbat wall, f44 and terminate
                f44(ii)=true;
                act(ii)=false;
            else
                ox   = tx-l*wx; % periodic boundary, shift to opposite wall
                oy   = ty-l*wy;
                oz   = tz;
            
                nr   = nr+1; % ++ reflect counter
            end
        end
        
        end

    end
end

na = sum(act);
% disp([num2str(100*na/n) '% rays unterminated!'])
n = n - na;

% f41c = sum(f41c)/n;
% disp([num2str(abs(f41-f41c)/f41) '% f14/f41 error'])
f42 = sum(f42)/n;
f43 = sum(f43)/n;
f44 = sum(f44)/n;

f24 = f42*a4/a2;
f34 = f43*a4/a3;

f32 = 1-f31-f34;
f23 = f32*a3/a2;

f22 = 1-f21-f23-f24;

end

vfr = [0 f12 f13 f14 f21 f22 f23 f24 f31 f32 0 f34 f41 f42 f43 f44];

end

% C16.m returns the view factor between two adjacent rectangles, of
% shared length h, width of surface 2 a, width of surface 1 b, and angle
% between them phi. From UT Austin, C-16, thermalradiation.net
function f12=C16(a,b,h,phi)
A=a/h;
B=b/h;
C=A^2+B^2-2*A*B*cos(phi);
D=sqrt(1+A^2*sin(phi)^2);
syms xsi

% then, breaking up the giant term into several smaller ones:
t21=A*B*sin(phi)+(pi/2-phi)*(A^2+B^2)+...
    B^2*atan((A-B*cos(phi))/(B*sin(phi)))+...
    A^2*atan((B-A*cos(phi))/(A*sin(phi)));
t22=(2/sin(phi)^2-1)*log((1+A^2)*(1+B^2)/(1+C))+...
    B^2*log((B^2*(1+C))/(C*(1+B^2)))+...
    A^2*log((A^2*(1+A^2)^(cos(2*phi)))/(C*(1+C)^(cos(2*phi))));
t23=1/pi*atan(1/B)+A/pi/B*atan(1/A)-sqrt(C)/pi/B*atan(1/sqrt(C));
t24=atan(A*cos(phi)/D)+atan((B-A*cos(phi))/D);
t2h=sqrt(1+xsi^2*sin(phi)^2);
t25=t2h*(atan(xsi*cos(phi)/t2h)+atan((A-xsi*cos(phi))/t2h));
f12=-sin(2*phi)/4/pi/B*t21+...
    sin(phi)^2/4/pi/B*t22+...
    t23+...
    sin(phi)*sin(2*phi)/2/pi/B*A*D*t24+...
    cos(phi)/pi/B*int(t25,xsi,0,B);
f12=double(f12);

end

% C11.m returns the view factor between two directly opposed, parallel
% rectangles, of width a, length b, and distance between them c. 
% From UT Austin, C-11, thermalradiation.net
function f12 = C11(a,b,c)

x = a/c;
y = b/c;

f12 = 2/pi/x/y*(1/2*log((1+x^2)*(1+y^2)/(1+x^2+y^2))+...
                x*sqrt(1+y^2)*atan(x/sqrt(1+y^2))+...
                y*sqrt(1+x^2)*atan(y/sqrt(1+x^2))-...
                x*atan(x)-y*atan(y));

end

% C121.m returns the view factor between a sphere distance d away from a 
% rectangle of width w and length l, centered on and perpendicular to the 
% diameter of the sphere. From UT Austin, C-121, thermalradiation.net
function f12=C121(d,l,w)
B=((l/2)/d);
C=((w/2)/d);
n1=(2*B^2-(1-B^2)*(B^2+C^2))/((1+B^2)*(B^2+C^2));
n2=(2*C^2-(1-C^2)*(B^2+C^2))/((1+C^2)*(B^2+C^2));
f12 = 1/2/pi*(asin(n1)+asin(n2));
end