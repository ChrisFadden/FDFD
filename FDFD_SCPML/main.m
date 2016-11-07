clear;
clc;

%Main Function for FDFD
%%  Simulation Parameters
fn = '../device/Free Space/';
pmlFlag = true;
srcAngle = 90;
pol = 'TE';

%%  Grid Calculation
grid = setupGrid(strcat(fn,'Grid.dat'));
NLambda = grid.lam0 / grid.dx;

%%  Device Calculation (UPML)
 percentPML = 10 * 0.01;
 pmlX = round(percentPML * grid.Nx);
 pmlY = round(percentPML * grid.Ny);

[pml.sx, pml.sy] = calcpml(grid,pmlX,pmlY);

if(~pmlFlag)
    pml.sx = 1;
    pml.sy = 1;
end

device = [];
device = setupDevice(fn,device,pml);

Sx = diag(sparse(1 ./ pml.sx(:)));
Sy = diag(sparse(1 ./ pml.sy(:)));

%%  Derivative Operator
%Derivative Matrices
[DEX,DEY,DHX,DHY] = yeeder(grid);

switch pol
    case 'TM',
        A = Sx*DHX/device.URyy*Sx*DEX + Sy*DHY/device.URxx*Sy*DEY + device.ERzz;
    case 'TE',
        A = Sx*DEX/device.ERyy*Sx*DHX + Sy*DEY/device.ERxx*Sy*DHY + device.URzz;
    otherwise,
        error('Unrecognized polarization.');
end

%%  Source Calculation
src = setupSrc(grid,A,pmlX,pmlY,srcAngle);

%%  Compute Field
Psi = A\src;

Psi = full(Psi);
Psi = reshape(Psi,grid.Nx,grid.Ny);

subplot(2,2,3)
imagesc([0:grid.Nx-1]*grid.dx / (10e-6),[0:grid.Ny-1]*grid.dy / (10e-6),real(Psi)');
xlabel('x (\mum)');
ylabel('y (\mum)');
title('E_z FIELD (TM Polarized)');
colorbar;