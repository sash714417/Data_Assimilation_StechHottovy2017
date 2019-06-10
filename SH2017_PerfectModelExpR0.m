%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Script for perfect model experiment of SH2017 Model. %%%%%%%%%%%%%%%%%%%%%%%  
clear;
%close all;
figfontsize=14;
set(0,'defaultlinelinewidth',1.0)
close all
clc


%%

addpath ~/Dropbox/Math/Research/Stechmann/Matlab/Data_Assimilation/StechHott2017_DA/GlobalFuns/
% Windows
%addpath C:\Users\Scott\Dropbox\Math\Research\Stechmann\Matlab\NColumnCloud\CCW\BetaEffectEWasym\SWEMultMoist\GlobalFuns\


%% Constants


% Set modes
modes =3;
totaleqs = 5*modes+3;
num_q0 = 3*(modes+1);



%GRL Regime 
tau1 =1/2;
tau2 =4;
b1 = 0.1;
b2 = 0.8;
tildeQ1=.9;
tildeQ2 =.45;
taur =25;
taul =taur;
tauv =taur;
%D11 = 0;
D11 =1;
D12 =0;
D21 =0;
%D22 = 0;
D22 =1;
epsilon = 0;
b = b1;


%   Analytical Forcing Constant
F = 0;

% Define the diffusion matrix
D = [D11, D12; D21, D22];

% Set the meridional non-dimensional scale. 
ybar = 15*110/1490;

%%

% More constants Constants
p = 4*10^4;
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

% Print Dimensional Parameter values
b_low = b1*L^2/T;
b_mid = b2*L^2/T;
tau_low = T*tau1;
tau_mid = T*tau2;
D2_lowmid = D11^2*Q^2/T;
tau_u = taur*T;

dayK = 75;

% Print? if yes set to 1
printyes = 1;

% Define the contours
MJOcontours = linspace(.5,1,2);
CCKWcontours = linspace(.5,1,2);
MJOcontours_u = linspace(.5,1,2);
CCKWcontours_u = linspace(.5,1,2);
Precipcontours = linspace(.5,1,2);
precip_caxis = [.07,.2];

%% Define Kalman Filter Parameters

G = eye(totaleqs);

r0 = 1;


%% Analytical k=0 mode

minc=-100;
maxc=100;


 
% Define zonal lattice
  x = linspace(0,26.67,10^4);
  
 
 N = length(x); %1D lattice size

 % Number of wave numbers kept for the model
  WN =20;
 
  
  % Same rng number
  rng(875);

% define the total run time in days
%Tmax = 30*1/T*24; % 30 days
Tmax = 365*1/T*24; % Year
%Tmax = 3*365*1/T*24; %3 Years

% Define the time stepping
Tmin =24/T; % 1 per day

% Time vector
t = 0:Tmin:Tmax;
dt =Tmin;



% Create full vector Initial Condition?
U = zeros(length(t),totaleqs,N);
%U(1,:,:) = real(RealizationFourierQMode_wF(zeros(totaleqs,N),modes,10^4,N,0,tildeQ1,0.45,b1,b2,taur,tau1,tau2,[0.03, 0; 0 0.03],WN));

% define forecast and posterior...
forecast = zeros(size(U));
posterior = zeros(size(U));

% Define observations
v = zeros(size(U));


% Creat boxcar wavenumbers
Waves = [2:2];

% Wavenumbers
kvec = 1:1:1;

% Define R0
R0 = eye(totaleqs); % CORRECT???!!!

% Loop over times
tic
for i = 1:length(t)-1
    
    % Define Observations
    v(i,:,:) = squeeze(U(i,:,:)) + r0*randn(totaleqs, length(x)); % Single observation! V is a scalar
    
    
    % Fourier transform for obs
    v_k = zeros(totaleqs,length(kvec));
    
    % True Solution
    U_k = zeros(totaleqs,length(kvec)); % Allocate for fourier vector
    
    % Forecast Fourier
    forecast_k = zeros(size(U_k));
    
    % Posterior Fourier
    posterior_k = zeros(size(U_k));
    
    % Fourier transform of initial conditions
    Uhat0FullSpec = fft(squeeze(real(U(i,:,:))),[],2); % True solution
    forecast_kFullSpec = fft(squeeze(real(forecast(i,:,:))),[],2); % Forecast
    posterior_kFullSpec = fft(squeeze(real(posterior(i,:,:))),[],2); % Posterior
    v_kFullSpec = fft(squeeze(real(v(i,:,:))),[],2); % Observations
    
    % Pass through boxcar filter
    Uhat0 = Uhat0FullSpec(:,Waves); % True
    forecasthat0 = forecast_kFullSpec(:,Waves); % Forecast
    posteriorhat0 = posterior_kFullSpec(:,Waves); % Posterior
    vhat0 = v_kFullSpec(:,Waves); % Obs
    
    % Shift
    Uhat0 = fftshift(Uhat0,2); % True
    forecasthat0 = fftshift(forecasthat0,2); % Forecast
    posteriorhat0 = fftshift(posteriorhat0,2); % Posterior
    vhat0 = fftshift(vhat0,2);
    
    
    % Define filter for different k
    rpostk = zeros(length(kvec),totaleqs,totaleqs); % Posterior covariance matrix
    rforek = zeros(length(kvec),totaleqs,totaleqs); % forecast covariance matrix
    Kk = zeros(length(kvec),totaleqs,totaleqs); % Kalman Gain matrix
         
    %loop over wave numbers
    for j = 1:length(kvec)  
          
        k = kvec(j); % Define wave number
        
        % Call on function for the model
        F_k = SH2017_LinearModel_DA_F(k,dt,modes,tildeQ1,tildeQ2,b1,b2,taur,tau1,tau2); % Model matrix

        % Call function for the covariance matrix
        Cov_Mat = Covariance_SH2017(k,dt,modes,tau1,tau2,b1,b2,tildeQ1,tildeQ2,taur,taul,tauv,length(x)*D);
                        
        % True Solution
        U_k(:,j) = F_k*Uhat0(:,j)+ mvnrnd(zeros(totaleqs,1),Cov_Mat)'; % Run the model with Gaussian noise defined by a covariance matrix
 
        % Forecast
        forecast_k(:,j) = F_k*posteriorhat0(:,j);
        
        % define r-forecast (prior error covariance)
        rforek(j,:,:) = F_k*squeeze(rpostk(j,:,:))*F_k'+Cov_Mat;
        
        Kk(j,:,:) = squeeze(rforek(j,:,:))*G'/(G*(squeeze(rforek(j,:,:)))*G'+R0);
    
        posterior_k(:,j) = forecast_k(:,j)+squeeze(Kk(j,:,:))*(vhat0(:,j)-G*forecast_k(:,j));

        rpostk(j,:,:) = (eye(size(squeeze(Kk(j,:,:))*G))-squeeze(Kk(j,:,:))*G)*squeeze(rforek(j,:,:));
    end
     
    % Inverse transform for true solution U, posterior, and forecast 
    Ufull_allwaves = zeros(size(squeeze(U(i,:,:)))); % Allocate for solution with padded zeros
    posterior_allwaves = zeros(size(squeeze(U(i,:,:))));
    forecast_allwaves = zeros(size(squeeze(U(i,:,:))));
    
    Ufull_allwaves(:,Waves) = ifftshift(U_k,2); % Set the boxcar wavenumbers to model
    posterior_allwaves(:,Waves) = ifftshift(posterior_k,2); %Set the boxcar wavenumbers to model
    forecast_allwaves(:,Waves) = ifftshift(forecast_k,2); %Set the boxcar wavenumbers to model

    U(i+1,:,:) =ifft(Ufull_allwaves,[],2,'symmetric'); % Take the inverse transform
    posterior(i+1,:,:) = ifft(posterior_allwaves,[],2,'symmetric');
    forecast(i+1,:,:) = ifft(forecast_allwaves,[],2,'symmetric');

end
toc

Ufull = U;



%% Plots First plot in Fourier Space

PF=figure;
set(PF, 'Position', [100, 100, 1049, 895]);

% 
R0FullSpec = fft(squeeze(real(U(:,1,:))),[],2); % True solution
ObsR0FullSpec = fft(squeeze(real(v(:,1,:))),[],2); % Obs
PostR0FullSpec = fft(squeeze(posterior(:,1,:)),[],2); % Posterior
forecastR0FullSpec = fft(squeeze(forecast(:,1,:)),[],2); % Forecast

subplot(2,1,1)
%plot k=1
plot(T*t/24,real(R0FullSpec(:,2)),'--')
hold on
plot(T*t/24,real(ObsR0FullSpec(:,2)),'o')
plot(T*t/24,real(PostR0FullSpec(:,2)),'-')

%Make it look nice
title(strcat(['Time series of true filter for r_0(t) and k=1 in Fourier space']));
legend('true (\hat{r}_0)_n','obs \hat{v}_n','posterior (\hat{r}_0)_{n|n}','Location','Northeast','Orientation','Horizontal');
xlabel('time [days]');
ylabel('solution');


% Make RMS plot
subplot(2,1,2)

% Define erros
RMS_error_TruePost = abs(real(R0FullSpec(:,2))-real(PostR0FullSpec(:,2)));
RMS_error_TrueFore = abs(real(R0FullSpec(:,2))-real(forecastR0FullSpec(:,2)));
plot(T*t/24,RMS_error_TruePost,'-k')
hold on
plot(T*t/24,1*ones(size(t)),'--k')
plot(T*t/24,RMS_error_TrueFore,'-k','Marker','.')

RMS_error_TruePostTime = norm(real(R0FullSpec(:,2))-real(PostR0FullSpec(:,2)))/sqrt(length(t(300:end)));


title(strcat(['Time series for the RMS with Temporal RMS = '],num2str(RMS_error_TruePostTime)));
legend('posterior (\hat{r}_0)_{n|n}','obs error','no filter','Location','Northeast','Orientation','Horizontal');
xlabel('time');
ylabel('RMS Error');

%% Print

% Print
prtstr = strcat('~/Dropbox/Math/Research/Stechmann/Data_Assimilation/Images/SH2017_TrueModel_Fourier_R0.eps');
 
print( gcf, '-depsc2', prtstr )


