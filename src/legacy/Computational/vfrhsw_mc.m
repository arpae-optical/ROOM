function vfr = vfrhsw_mc(d_l,h_l,n)

% geometric parameters
    
l = 1;        % length scale of unit cell
d = d_l*l;    % diameter of central sphere
r2= d*d/4;    % radius squared
h = h_l*l;    % height of adiabatic walls
s = l/sqrt(3);% length of one side of hex
a1 = 3/2*s*l;% area of substrate
a3 = a1;    % apparent area
a2 = pi*d^2; % area of sphere
a4 = h*s*6;  % area of thin walls

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

% loop through n, determine incident surfaces for each bundle

for ii = 1:n
    ox = s/2*mcx(ii)*sqrt(mcy(ii)); % origin from x=[0 s/2*y/ymax]
    oy = l/2*sqrt(mcy(ii));         % origin from y=[0 l/2]
    oz = -d/2;                      % origin at z=-d/2
    np = 2*pi*mcp(ii);              % phi from [0 2pi]
    nq = asin(sqrt(mcq(ii)));       % the from [0 pi/2]
    nx = sin(nq)*cos(np);           % pointing x
    ny = sin(nq)*sin(np);           % pointing y
    nz = cos(nq);                   % pointing z
    
    nr = 0; % reset reflection count
    
    while nr<nrm&&act(ii) % while bundle is active, until max reflection
        
        % first, check for f12: does bundle intercept sphere?
        % project the distance to center of sphere onto direction of bundle
        pc = ox*nx+oy*ny+oz*nz; 
        if pc^2<=r2 % if the projected incerception distance is less than r
            f12(ii) = true;  % bundle is part of f12
            act(ii) = false; % and the bundle is no longer active
            
        % otherwise, check for f13; does bundle leave the unit cell?
        else
            s3 = (-d/2+max([h,d])-oz)/nz; % scale to reach apparent area is dz/nz
            tx = ox+nx*s3;
            ty = oy+ny*s3;
%             tz = oy*ny*s3;
            if abs(ty)<l/2&&abs(tx)<s*(abs(ty)/(l/2)) % if inside hex, f13
                f13(ii) = true;
                act(ii) = false;
                
            % otherwise check for f14; does bundle intercept adbat wall?
            else
                % time to determine which of the six walls we hit first:
                % radius of our circle intercept is s:
                aa = (nx^2+ny^2);
                bb = (ox*nx+oy*ny);
                cc = (ox^2+oy^2-s^2);
                ss = (-bb+sqrt(bb^2-4*aa*cc))/(2*aa);
                % the location on that circle, then, is given by:
                tx = ox+nx*ss; % x intercept of circle
                ty = oy+ny*ss; % y intercept of circle
                qq = (x<0)*pi-atan(ty/tx); % angle to intercept
                % now, we use this value to find the angle of the wall
                % intercepted:
                qq = pi/3*floor((qq-pi/6)/(pi/3))+pi/6; 
                wx = -sign(tx)*abs(cos(qq)); % normal vector of plane
                wy = -sign(ty)*abs(sin(qq));
                tx = cos(qq);
                ty = sin(qq);
                
                ss = ((tx-ox)*wx+(ty-oy)*wy)/(wx*nx+wy*ny);
                
                ox = ox+ss*nx;
                oy = oy+ss*ny;
                oz = oz+ss*nz;
                
                % check to see if z intercept lower than wall
                if oz<-d/2+h
                    f14(ii)=true;
                    act(ii)=false;
                    
                % otherwise, reflect    
                else
                    ox = ox+l*wx;
                    oy = oy+l*wy;
                    nr = nr+1;        % iterate number of reflections
                end
            end
        end
    end
end

na = sum(act);
disp([num2str(100*na/n) '% rays unterminated!'])
n = n - na;

f11 = 0;
f12 = sum(f12)/n;
f13 = sum(f13)/n;
f14 = sum(f14)/n;

f21 = f12*a1/a2;
%f22 = ?
%f23 = ?
%f24 = ?

f31 = f13*a1/a3;
%f32 = ?
f33 = 0;
%f34 = ?

f41 = f14*a1/a4;
%f42 = ?
%f43 = ?
%f44 = ?

% for remaining view factors, if walls are higher than the sphere, you can
% utilize analytic view factors to get the missing view factors:
if h>=d
    f22 = 0;
    f24 = 3*(vfsr(d,s,h)+vfsr(d,s,d));
    f23 = 1-f21-f22-f24;
    
    f32 = f23*a2/a3;
    f34 = 1-f31-f32-f33;
    
    f42 = f24*a2/a4;
    f43 = f34*a3/a4;
    f44 = 1-f41-f42-f43;
    
% otherwise, we'll need to do a second monte carlo simulation
else
act = true(1,n); % is bundle ii active?
f41 = false(1,n);% is bundle ii incident on surface 1?
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

for ii = 1:n
    ox = s*(-1/2+mcx(ii));          % origin from x=[-s/2 s/2]
    oy = -l/2;                      % origin at y=-L/2
    oz = -d/2+h*mcz(ii);            % origin from z=[-d/2 h-d/2]
    np = 2*pi*mcp(ii);              % phi from [0 2pi]
    nq = asin(sqrt(mcq(ii)));       % the from [0 pi/2]
    nx = sin(nq)*cos(np);           % pointing x
    ny = cos(nq);                   % pointing y
    nz = sin(nq)*sin(np);           % pointing z
    
    nr = 0; % reset reflection count
    
    while nr<nrm&&act(ii) % while bundle is active, until max reflection
        
        % first, check for f42: does bundle intercept sphere?
        % project the distance to center of sphere onto direction of bundle
        pc = ox*nx+oy*ny+oz*nz; 
        if pc^2<=r2 % if the projected incerception distance is less than r
            f42(ii) = true;  % bundle is part of f42
            act(ii) = false; % and the bundle is no longer active
            
        % otherwise, check for f41/3; does bundle leave the unit cell?
        else
            s3 = (-d/2+(nz>0)*max([h,d])-sign(nz)*oz)/nz; % scale to reach apparent area is dz/nz
            tx = ox+nx*s3;
            ty = oy+ny*s3;
%             tz = oy*ny*s3;
            if abs(ty)<l/2&&abs(tx)<s*(abs(ty)/(l/2)) % if inside hex, f41/3
                if nz<0
                    f41(ii) = true;
                else
                    f43(ii) = true;
                end
                act(ii) = false;
                
            % otherwise check for f14; does bundle intercept adbat wall?
            else
                % time to determine which of the six walls we hit first:
                % radius of our circle intercept is s:
                aa = (nx^2+ny^2);
                bb = (ox*nx+oy*ny);
                cc = (ox^2+oy^2-s^2);
                ss = (-bb+sqrt(bb^2-4*aa*cc))/(2*aa);
                % the location on that circle, then, is given by:
                tx = ox+nx*ss; % x intercept of circle
                ty = oy+ny*ss; % y intercept of circle
                qq = (x<0)*pi-atan(ty/tx); % angle to intercept
                % now, we use this value to find the angle of the wall
                % intercepted:
                qq = pi/3*floor((qq-pi/6)/(pi/3))+pi/6; 
                wx = -sign(tx)*abs(cos(qq)); % normal vector of plane
                wy = -sign(ty)*abs(sin(qq));
                tx = cos(qq);
                ty = sin(qq);
                
                ss = ((tx-ox)*wx+(ty-oy)*wy)/(wx*nx+wy*ny);
                
                ox = ox+ss*nx;
                oy = oy+ss*ny;
                oz = oz+ss*nz;
                
                % check to see if z intercept lower than wall
                if oz<-d/2+h
                    f44(ii)=true;
                    act(ii)=false;
                    
                % otherwise, reflect    
                else
                    ox = ox+l*wx;
                    oy = oy+l*wy;
                    nr = nr+1;        % iterate number of reflections
                end
            end
        end
    end
end

na = sum(act);
disp([num2str(100*na/n) '% rays unterminated!'])
n = n - na;

f41 = sum(f41)/n;
disp([num2str(abs(f41-f14*a1/a4)/f41) '% f14/f41 error!'])
f42 = sum(f42)/n;
f43 = sum(f43)/n;
f44 = sum(f44)/n;

f24 = f42*a4/a2;
f34 = f43*a4/a3;
f32 = 1-f31-f33-f34;
f23 = f32*a3/a2;
f22 = 1-f21-f23-f24;

end

vfr = [f11 f12 f13 f14 f21 f22 f23 f24 f31 f32 f33 f34 f41 f42 f43 f44];

end