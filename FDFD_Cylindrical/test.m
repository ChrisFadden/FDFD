clear;
clc;
%%
%   Cylindrical Test
%
    freq = 2 * 10^9;
    cc = 3*10^8;
    lam0 = cc / freq;
    
    k_rho = 2*pi / lam0;
    Nlambda = 20;
    sizeLambda = 5;
    sizeR = sizeLambda * lam0;
    numR = Nlambda * sizeLambda;
    d_rho = sizeR / numR;
    rho = d_rho:d_rho:sizeR;
    k_rho = 1;
    Psi_test = besselj(0,k_rho * rho);
    PsiDerTest = 0 * Psi_test + 1;
            
    for rr = 2:numR-1
        %Psi''
        SecondDer(rr) = (Psi_test(rr + 1) + Psi_test(rr - 1) - 2*Psi_test(rr)) / (d_rho^2);
        
        %(1/r) Psi'
        FirstDer(rr) = (1 / (rr*d_rho)) * ((Psi_test(rr + 1) - Psi_test(rr)) / (d_rho));
        
        %(k rho)^2 Psi
        constant(rr) = (rr*d_rho * k_rho)^2 * Psi_test(rr);
        
        %Psi'' + (1/r)Psi' + (k rho)^2 Psi = 0
        PsiDerTest(rr) = SecondDer(rr) + FirstDer(rr) + constant(rr);
    end
    testSol = bessely(0,rho(2:end-1));   
    
    plot(PsiDerTest(2:end-1) - testSol)  
    
    
    