%% DNSTest.m
% arguments: N, L/H
% Just a test script for analyzing performance of DNS for a hexagonal
% enclosure. Will be updated to include a cylindrical occlusion.
% N: number of subdivisions on each square division
% L/H: edge to edge distance of hex divided by the height
% Geometry is of form:
%   |Y             |Z
%  ___            ___
% / O \  _X     || O || _X
% \___/         ||___||
%

function FWW=DNSTest(N,L_H)
H   = 1;         % Height of walls
L   = H*L_H;     % Edge-to-edge distance of walls
s   = L/sqrt(3); % Edge length of walls

gsource = g_quad([-s/2,-L/2,-H/2],[0,0,H],[ s,0,0],N);
% gtarget = g_quad([ s  ,   0,-H/2],[0,0,H],[-s/2,L/2,0],N);
gtarget = g_hextest(N,L_H);

FWW     = DNS(gsource,gtarget);
end
