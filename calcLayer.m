function Si = calcLayer(DEX,DEY,DHX,DHY,W0,V0,grid,device,layer)
  
  %Get parameters for current layer
  ERCxx = device.ERCxx(:,:,layer);
  ERCyy = device.ERCyy(:,:,layer);
  ERCzz = device.ERCzz(:,:,layer);
  
  URCxx = device.URCxx(:,:,layer);
  URCyy = device.URCyy(:,:,layer);
  URCzz = device.URCzz(:,:,layer);
  
  P00 = -DEX * inv(ERzz) * DHY;
  P01 = (URyy + DEX*inv(ERzz)*DHX);
  P10 = -(URxx + DEY * inv(ERzz)*DHY);
  P11 = DEY * inv(ERzz) * DHX;
  
  P = [P00 P01; P10 P11];
  
  Q00 = -DHX * inv(URzz) * DEY;
  Q01 = (ERyy + DHX*inv(ERzz)*DEX);
  Q10 = -(ERxx + DHY * inv(URzz)*DEY);
  Q11 = DHY * inv(ERzz) * DEX;
    
  Q = [Q00 Q01; Q10; Q11];
  
  %Compute Eigenmodes of layer
  %Eigenvalues /vectors can be out of order, but as long as they are
  %together, the end Scattering Matrices will be the same.
  OM2 = P*Q;
  [W, LAM2] = eig(OM2);

  LAM = sqrt(LAM2);
  V = Q*W*inv(LAM);

  %Compute S-matrix temp variables
  A = inv(W) * W0 + inv(V)*V0;
  B = inv(W) * W0 - inv(V)*V0;
  X = expm(-LAM * (2*pi / grid.lam0) * grid.L(layer));
      
  Si.S11 = inv(A - X*B*inv(A)*X*B)*(X*B*inv(A)*X*A - B);
  Si.S12 = inv(A - X*B*inv(A)*X*B)*X*(A - B*inv(A)*B);
  Si.S21 = Si.S12;
  Si.S22 = Si.S11;  
end
