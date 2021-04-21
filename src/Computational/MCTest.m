%% MCTest.m
% arguments: N, L/H
% N: number of bundles
% L/H: edge to edge distance of hex divided by the height
% Geometry is of form:
%   |Y             |Z
%  ___            ___
% / O \  _X     || O || _X
% \___/         ||___||
%
function [FWW,FWT,FWB,t]=MCTest(N,L_H)
tic
H   = 1;         % Height of walls
L   = H*L_H;     % Edge-to-edge distance of walls
s   = L/sqrt(3); % Edge length of walls

MCX = rand(1,N); % Monte-Carlo RNG for x distribution of bundles
MCY = ones(1,N); % "Distribution" of y locations: all at Y=-L/2
MCZ = rand(1,N); % Monte-Carlo RNG for z distribution of bundles
MCP = rand(1,N); % Monte-Carlo RNG for phi distribution of bundles
MCQ = rand(1,N); % Monte-Carlo RNG for theta distribution of bundles

OX  = (-1/2+MCX)*s;     % Origin of bundles
OY  = (-L/2)*MCY;       % Origin of bundles
OZ  = (-1/2+MCZ)*H;     % Origin of bundles
NP  = 2*pi*MCP;         % Phi of bundles
NQ  = pi/2*MCQ;         % Theta of bundles
NX  = sin(NQ).*cos(NP); % X direction of bundles
NY  = cos(NQ);          % Y direction of bundles
NZ  = sin(NQ).*sin(NP); % Z direction of bundles

DZ  = sign(NZ)*H/2-OZ;  % Distance to intercept Z=+/-(H/2)

S   = DZ./NZ;           % Scale of N to intercept ^

IX  = OX+S.*NX;         % X location at |Z|=H/2
IY  = OY+S.*NY;         % Y location at |Z|=H/2

% If Y@|Z|=H/2>L/2, intercepts far wall.
% If |X|>S-S/L*|Y|, intercepts side wall.
% Otherwise, if n.Z>0, intercepts top plane.
% Otherwise, if n.Z<0, intercepts bottom plane.

NW  = abs(IY)>L/2|abs(IX)>s-s/L*(abs(IY)); 
NT  = ~NW&NZ>0;
NB  = ~NW&NZ<0;

FWW = sum(NW)/N; % The view factor is the number of bundles hit divided by total
FWT = sum(NT)/N;
FWB = sum(NB)/N;

t=toc;
end
