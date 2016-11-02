function [W0, V0] = calcFreeSpace(DEX,DEY,DHX,DHY)
  
  P00 = -DEX * DHY;
  P01 = DEX * DHX;
  P01 = P01 + eye(size(P01,1));
  P10 = DEY * DHY;
  P10 = -(eye(size(P10,1)) + P10);
  P11 = DEY * DHX;
  
  P = [P00, P01; P10, P11];
  
  Q00 = -DHX * DEY;
  Q01 = DHX*DEX;
  Q01 = Q01 + eye(size(Q01,1));
  Q10 = DHY * DEY;
  Q10 = -(eye(size(Q10,1)) + Q10);
  Q11 = DHY * DEX;
    
  Q = [Q00 Q01; Q10, Q11];
    
  PQ = full(P*Q);
  
  [W0, LAM2] = eig(PQ);
  LAM = sqrt(LAM2);
  V0 = Q*W0*inv(LAM);
end
