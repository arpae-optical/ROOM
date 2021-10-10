%% g_hextest.m
% arguments: N, L/H
% creates quad mesh for five walls of hexagonal enclosure (2&3,L&R, and 4)
% N: number of subdivisions on each square division
% L/H: edge to edge distance of hex divided by the height
% Geometry is of form:
%    |Y             |Z
% L _4_ R          _T_
% 3/ O \3 _X     || O || _X
% 2\___/2        ||___||
%    1              B
 
function comb=g_hextest(N,L_H)
H   = 1;         % Height of walls
L   = H*L_H;     % Edge-to-edge distance of walls
s   = L/sqrt(3); % Edge length of walls

g2r = g_quad([ s/2,-L/2,-H/2],[0,0,H],[ s/2, L/2,0],N); % quadmesh for 2R
g3r = g_quad([ s  ,   0,-H/2],[0,0,H],[-s/2, L/2,0],N); % 3R
g4  = g_quad([ s/2, L/2,-H/2],[0,0,H],[-s  ,   0,0],N); % 4
g3l = g_quad([-s/2, L/2,-H/2],[0,0,H],[-s/2,-L/2,0],N); % 3L
g2l = g_quad([-s  ,   0,-H/2],[0,0,H],[ s/2,-L/2,0],N); % 2L

comb.xyz=cat(1,g2r.xyz,g3r.xyz,g4.xyz,g3l.xyz,g2l.xyz); % cat the surfaces
comb.n  =cat(1,g2r.n,  g3r.n,  g4.n,  g3l.n,  g2l.n  ); % to make the 
comb.A  =cat(1,g2r.A,  g3r.A,  g4.A,  g3l.A,  g2l.A  ); % target surface
end
