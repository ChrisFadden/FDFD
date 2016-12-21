clear;
clc;

%Main Function for FDFD
%%  Simulation Parameters
fn = '../device/Free Space/';
pmlFlag = true;
srcAngle = 0;

%%  Grid Calculation
grid = setupGrid(strcat(fn,'Grid.dat'));
NLambda = grid.lam0 / grid.dx;

%%  Device Calculation (UPML)
 percentPML = 10 * 0.01;
 pmlX = round(percentPML * grid.Nx);
 pmlY = round(percentPML * grid.Ny);

[pml.sx, pml.sy] = calcpml(grid,pmlX,pmlY);

if(~pmlFlag)
    pml.sz = 1;
    pml.sr = 1;
end

device = [];
tic
device = setupDevice(fn,device,pml);

Sz = diag(sparse(1 ./ pml.sx(:)));
Sr = diag(sparse(1 ./ pml.sy(:)));

%%  Derivative Operator
%Derivative Matrices
[DEZ,DER,DHZ,DHR] = yeeder(grid);

%rho = grid.dy:grid.dy:grid.Ly;
%dr = (2*pi/grid.lam0)*grid.dy;
%rho = dr:dr:grid.Ly*(2*pi/grid.lam0);
rho = (1:grid.Ny)./grid.dy;

invRho = speye(grid.Nx*grid.Ny);
jj = 1;
s = 0;
for ii = 1:grid.Nx*grid.Ny
    invRho(ii,ii) = 1/rho(jj);
    if(jj == length(rho))
        jj = 0;
    end
    jj = jj+1;
end

Az = Sz*DHZ/device.URyy*Sz*DEZ;
Ar1 = Sr*DHR/device.URxx*Sr*DER;
Ar2 = Sr*DHR/device.URxx*invRho;
A = Az + Ar1 + Ar2 + device.ERzz;
%A = Ar1 + Ar2 + device.ERzz;

%%  Source Calculation
src = setupSrc(grid,A,pmlX,pmlY,srcAngle);

%%  Compute Field
Psi = A\src;
toc
Psi = full(Psi);
Psi = reshape(Psi,grid.Nx,grid.Ny);

imagesc([0:grid.Nx-1]*grid.dx / (10e-6),[0:grid.Ny-1]*grid.dy / (10e-6),real(Psi)');
xlabel('Z (\mum)');
ylabel('\rho (\mum)');
title('E_{\phi\phi} FIELD (TM Polarized)');
colorbar;