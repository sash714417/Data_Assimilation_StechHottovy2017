function [qt] = RealizationFourierQMode_wF(U0,modes,t,N,F,tildeQ1,tildeQ2,b1,b2,taur,tau1,tau2,D_mat,WN)
% Constants

% Characteristic Scales
T = 8; %hrs
L = 1500; %kms
c = 50; %m/s
p = 40000; %km circumference of the earth. 

% Number of equations
totaleqs = 5*modes+3;
num_q0 = 3*(modes+1);

taul = taur;
tauv = taur;

% Time and Space parameters
omega = -WN*2*pi*L/p:2*pi*L/p:WN*2*pi*L/p;

% Create boxcar wavenumbers
Waves = [1:WN+1 N-(WN-1):N];

%Force = [zeros(num_q0,N); F*ones(2,1)*sin(2*pi*(linspace(0,1,N)-1/4)); zeros(2*(modes-1),N)];
Force = [zeros(num_q0,N); F*ones(2,1)*(exp(-(linspace(0,1,N)-1/3).^4/(2*0.003))+.75*exp(-(linspace(0,1,N)-5/6).^4/(2*0.0001))-.5); zeros(2*(modes-1),N)];
ForceFT = fft(Force,[],2);

Uhat0FullSpec = fft(U0,[],2); % Fourier transform of Q. 

% Pass through boxcar filter
Uhat0 = Uhat0FullSpec(:,Waves);
Forcehat = ForceFT(:,Waves);

Uhat0 = fftshift(Uhat0,2);
Forcehat = fftshift(Forcehat,2);

% Define U as the fourier solution
Ufinal = zeros(totaleqs,length(omega));

% Define Utilde as the diagonalized solution S^{-1}U
Utilde = zeros(totaleqs,length(omega));
Forcetilde = zeros(totaleqs,length(omega));

WdotFourierFullSpec = fft([zeros(num_q0,N);randn(2*modes,N)],[],2);
% apply boxcarfilter
WdotFourier = WdotFourierFullSpec(:,Waves);
WdotFourier = fftshift(WdotFourier,2);

for i=1:length(omega)
    
    % Define wave number
    k = omega(i);
        
    % Define matrix (model)
    A_k = SWEModelAk(modes,k,tau1,tau2,b1,b2,tildeQ1,tildeQ2,taur,taul,tauv);

    % Find the diagonal matrix of eigenvalues (D_k) and matrix which
    % diagonalizes (S_k)
    [S_k,D_k] = eig(A_k);
    
    % Define \tilde{U} = S_k^{-1} Uhat
    Uhat0tilde = S_k\Uhat0(:,i);
    Forcetilde = S_k\Forcehat(:,i);
    
    % Define diffusion matrix
    D_start = zeros(5*modes+3);
    
    
    % Set diffusion matrix
    for l=1:modes
        D_start(num_q0+2*(l-1)+1:num_q0+2*(l-1)+2,num_q0+2*(l-1)+1:num_q0+2*(l-1)+2) = D_mat;
    end
    
    % define \tilde{D} = S_k^{-1}D
    DMat = S_k\D_start;
     
    % White noise vector
    WienerVector = DMat*WdotFourier(:,i);
    
    for j=1:totaleqs
    
        Utilde(j,i) = exp(D_k(j,j)*t).*Uhat0tilde(j)+(1/(D_k(j,j)))*(exp(D_k(j,j)*t)-1).*Forcetilde(j)...
            +sqrt((exp(2*D_k(j,j)*t)-1)/(2*D_k(j,j)))*WienerVector(j);   
        
        
        
    end
   
    
    
    Ufinal(:,i) = S_k*Utilde(:,i);   
    %end

end
%[Y,I] = max(Ufinal,N,2);


Ufull_allwaves = zeros(size(U0));

Ufull_allwaves(:,Waves) = ifftshift(Ufinal,2);

qt = ifft(Ufull_allwaves,[],2,'symmetric');


    
end