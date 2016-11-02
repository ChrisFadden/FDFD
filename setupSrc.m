function [ src ] = setupSrc( fileName,grid)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

degrees = pi / 180;

src.theta = 0 * degrees;   %Elevation Incident Angle
src.phi = 0 * degrees;       %Azimuthal Incident Angle
src.pte = 1;               %Source Polarization
src.ptm = 0;               %Source Polarization

%Polarization unit vector
n = [0; 0; 1];
k = [sin(src.theta)*cos(src.phi); sin(src.theta)*sin(src.phi); cos(src.theta)];

if(src.theta == 0)
   src.ate = [0; 1; 0]; 
else
    %multiply by -1 to match bench...   check this
    src.ate = -1 * cross(n,k) ./ norm(cross(n,k));
end

src.atm = cross(src.ate,k) / norm(cross(src.ate,k));

Esrcx = src.ate(1)*src.pte + src.atm(1)*src.ptm;
Esrcy = src.ate(2)*src.pte + src.atm(2)*src.ptm;

x = ((0:grid.Nx-1)*grid.dx)';
y = ((0:grid.Ny-1)*grid.dy)';
src.esrc = [Esrcx * exp(-1j*(k(1) * x + k(2) * y)); Esrcy * exp(-1j*(k(1) * x + k(2)*y))];
end

