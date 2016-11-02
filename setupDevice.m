function [ device] = setupDevice(fn,grid,device,layer,pml)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    er_fn = 'ER_Layer';
    fext = '.dat';
    er_fn = strcat(fn,er_fn);
    er_fn = strcat(er_fn,num2str(layer));
    er_fn = strcat(er_fn,fext);
    
    ER(:,:,1) = importdata(er_fn);
        
    device.ERCxx(:,:,layer) = pml.sy ./ pml.sx .* ER(:,:,1);
    device.ERCyy(:,:,layer) = pml.sx ./ pml.sy .* ER(:,:,1);
    device.ERCzz(:,:,layer) = pml.sx .* pml.sy .* ER(:,:,1);
        
    ur_fn = 'UR_Layer';
    ur_fn = strcat(fn,ur_fn);
    ur_fn = strcat(ur_fn,num2str(layer));
    ur_fn = strcat(ur_fn,fext);
    
    UR(:,:,1) = importdata(ur_fn);
    device.URCxx(:,:,layer) = pml.sy ./ pml.sx .* UR(:,:,1);
    device.URCyy(:,:,layer) = pml.sx ./ pml.sy .* UR(:,:,1);
    device.URCzz(:,:,layer) = pml.sx .* pml.sy .* UR(:,:,1);
end

