function Cov_Mat = Covariance_SH2017(wave_number,dt,modes,tau1,tau2,b1,b2,tildeQ1,tildeQ2,taur,taul,tauv,D)
%
%   Function for computing the covariance matrix of the Stechmann & Hottovy
%   2017 (GRL) linear stochastic model. 
%
%
%   Inputs:
%               wave_number:                  Wavenumber (integer)
%               dt:              frequency of time step
%               modes:              Number of meridional modes (odd, positive
%                                       integer, e.g. 1,3,5,...)
%               tau1,...,Dmat:      parameters of the model
%
%   Outputs: 
%               Covariance_SH2017:  Covariance Matrix 


% Define number of equations and where moisture equations start
totaleqs = 5*modes+3;
num_q0 = 3*(modes+1);

% Constants
p = 40000; %km circumference of the earth. 
beta = 2.3*10^(-11); %m^{-1}s^{-1}
theta_ref = 300; %K
g = 9.8; % m/s
H = 16; % km
cp = 1.01*10^3; % J/(kg air*K)
Lv = 2.5*10^6; % J/(kg water)


%derived constants
N = sqrt(10^(-4)); % s^{-1}
c = N*H/pi*1000; % m/s
L = sqrt(c/beta)/1000; % km
T = L/c*1000/(60^2); % hrs
alpha = H*N^2*theta_ref*1000/(pi*g);
Q = cp*alpha*1000/Lv; %g/kg

% Define non-dimensional wave number

A_k = SWEModelAk(modes,2*pi*L/p*wave_number,tau1,tau2,b1,b2,tildeQ1,tildeQ2,taur,taul,tauv);
[S_k,D_k] = eig(A_k);
        
% Define diffusion matrix
D_start = zeros(5*modes+3);
    
    
% Set diffusion matrix
for l=1:modes
    D_start(num_q0+2*(l-1)+1:num_q0+2*(l-1)+2,num_q0+2*(l-1)+1:num_q0+2*(l-1)+2) = D;
end
    
% define \tilde{D} = S_k^{-1}D
DMat = S_k\D_start;
     
% White noise vector
Noise = zeros(totaleqs,1);
        
for n=1:totaleqs
    Noise(n) = sqrt((exp(2*D_k(n,n)*dt)-1)/(2*D_k(n,n)));
end

Noise = S_k*diag(Noise)*DMat;
        
Cov_Mat = real(Noise*Noise');        

end
