function [F] = SH2017_LinearModel_DA_F(wave_number,t,modes,tildeQ1,tildeQ2,b1,b2,taur,tau1,tau2)
% Constants

% Characteristic Scales
T = 8; %hrs
L = 1500; %kms
c = 50; %m/s
p = 40000; %km circumference of the earth. 

taul = taur;
tauv = taur;

% Time and Space parameters



% Define wave number
k = wave_number*2*pi*L/p;
        
% Define matrix (model)
A_k = SWEModelAk(modes,k,tau1,tau2,b1,b2,tildeQ1,tildeQ2,taur,taul,tauv);

[S_k,D_k] = eig(A_k);
    
F = S_k*diag(exp(diag(D_k*t)))/S_k;
   
   

    
end