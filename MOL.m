close all;
clc;
clear;

%%
%****************
% MOL Parameters
%****************
pmlFlag = false;
loc = './device/Free Space/';

%%
%%%%%%%%%%
%   Setup 
%%%%%%%%%%

%Grid, Source, and Device
grid = setupGrid(strcat(loc,'Grid.dat')); 
src = setupSrc('',grid);
Ex = src.esrc(1:grid.Nx);
Ey = src.esrc(grid.Nx+1:end);


device.ER(grid.Nx,grid.Ny,length(grid.L)) = 0;
device.UR = device.ER;

%PML 
percentPML = 5 * 0.01;
pmlx = round(percentPML * grid.Nx);
pmly = round(percentPML * grid.Ny);
[pml.sx, pml.sy] = calcpml2d([grid.Nx,grid.Ny],[pmlx,pmlx,pmly,pmly],pmlFlag);

if(~pmlFlag)
    pml.sx = 1;
    pml.sy = 1;
end
for layer = 1:length(grid.L)
    [device] = setupDevice(loc,grid,device,layer,pml);
end

%Derivative Matrices
[DEX,DEY,DHX,DHY] = yeeder(grid);

%%
%%%%%%%%%%%%%%%%%
%   Calculate
%%%%%%%%%%%%%%%%%
[W0,V0] = calcFreeSpace(DEX,DEY,DHX,DHY);
