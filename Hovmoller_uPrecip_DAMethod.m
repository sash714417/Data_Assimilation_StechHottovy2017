%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Script for producing Hovmoller plots. %%%%%%%%%%%%%%%%%%%%%%%  
clear;
%close all;
figfontsize=14;
set(0,'defaultlinelinewidth',1.0)


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
U(1,:,:) = real(RealizationFourierQMode_wF(zeros(totaleqs,N),modes,10^4,N,0,tildeQ1,0.45,b1,b2,taur,tau1,tau2,[0.03, 0; 0 0.03],WN));

% Loop to solve the system

% Creat boxcar wavenumbers
Waves = [1:WN+1 N-(WN-1):N];

% Wavenumbers
kvec = -WN:1:WN;

% Loop over times


tic
for i = 1:length(t)-1
    
    U_k = zeros(totaleqs,length(kvec)); % Allocate for fourier vector
    
    WdotFourierFullSpec = fft([zeros(num_q0,N);randn(2*modes,N)],[],2);
    % apply boxcarfilter
    WdotFourier = WdotFourierFullSpec(:,Waves);
    WdotFourier = fftshift(WdotFourier,2);

    Uhat0FullSpec = fft(squeeze(real(U(i,:,:))),[],2); % Fourier transform of initial conditions. 

    % Pass through boxcar filter
    Uhat0 = Uhat0FullSpec(:,Waves);
    Uhat0 = fftshift(Uhat0,2);
    
        
    %loop over wave numbers
    for j = 1:length(kvec)  
          
        k = kvec(j); % Define wave number
        
        % Call on function for the model
        F_k = SH2017_LinearModel_DA_F(k,dt,modes,tildeQ1,tildeQ2,b1,b2,taur,tau1,tau2); % Model matrix

        % Call function for the covariance matrix
        Cov_Mat = Covariance_SH2017(k,dt,modes,tau1,tau2,b1,b2,tildeQ1,tildeQ2,taur,taul,tauv,length(x)*D);
        
        U_k(:,j) = F_k*Uhat0(:,j)+ mvnrnd(zeros(totaleqs,1),Cov_Mat)'; % Run the model with Gaussian noise defined by a covariance matrix
 
    end
        
    Ufull_allwaves = zeros(size(squeeze(U(i,:,:)))); % Allocate for solution with padded zeros

    Ufull_allwaves(:,Waves) = ifftshift(U_k,2); % Set the boxcar wavenumbers to model

    U(i+1,:,:) =ifft(Ufull_allwaves,[],2,'symmetric'); % Take the inverse transform

    
end
toc

Ufull = U;



%% Average over y

 %from -15 S to 15 N
y = linspace(-ybar,ybar,50);
dy = y(2)-y(1);

%define the parabolic cylinder functions
par_cyl_funs = ParabolicCylinderFuns(y',modes);

qlow_t = zeros(size(squeeze(U(:,1,:))));
qmid_t = zeros(size(squeeze(U(:,1,:))));
% Define u_t using the 0th and 1st modes
u_t = sqrt(2)/2*sum(par_cyl_funs(:,1))*dy*squeeze(U(:,1,:)+U(:,5,:))+...
    sqrt(2)/2*sum(par_cyl_funs(:,2))*dy*squeeze(U(:,2,:)+U(:,8,:));

for i=1:modes
    qlow_t = qlow_t + sum(par_cyl_funs(:,i))*dy*squeeze(U(:,num_q0+2*(i-1)+1,:)); 
    qmid_t = qmid_t + sum(par_cyl_funs(:,i))*dy*squeeze(U(:,num_q0+2*i,:)); 
    
    if i > 2
        u_t = u_t + sqrt(2)/2*sum(par_cyl_funs(:,i))*dy*squeeze(U(:,1+3*(i),:)+U(:,8+3*(i),:));
    end
end

q_t = qlow_t;
Precip_t = qlow_t/tau1+qmid_t/tau2;
% for i=1:length(t)
%     plot(squeeze(real(U(i,:,:)))');
%     axis([1 N -15 15])
%     
%     pause(.3)
% end
Color = [120 120 120; 0 236 236; 1 160 246; 0 0 246; 0 255 0; 0 200 0;...
    0 144 0; 255 255 0; 231 192 0; 255 144 0; 255 0 0; 214 0 0; 192 0 0;...
    255 0 255; 153 85 201; 235 235 235];

ColorZhang = [0 0 246; 0 255 0; 0 255 0;0 255 0;0 255 0;0 255 0; 0 255 0; 0 255 0; 0 255 0;0 255 0;0 255 0;0 255 0; 0 255 0;255 0 0];

Color = Color/255;
ColorZhang = ColorZhang/255;



 
        %% Plot raw data of u and precip
        figfontsize=12;
        
PF=figure;
set(PF, 'Position', [100, 100, 1049, 895]);
       subplot(3,8,[4:5,12:13])
  
       % Define time and q plots for year
       
    %contour(-N/2:1:N/2-1,freq,log10(24*(30)avext))
    h = pcolor(L*x,T*t/24,Q*(qmid_t));
    colormap(jet)
    set(h, 'EdgeColor', 'none');
    %plot(freq,log10(psdxnum),'r*')
    %Compute analytical formula

    grid off
    box on
    
   
     titlestr = strcat('q_{mid} []');
     title(titlestr,'fontsize',figfontsize);        %  put a title on the plot    %  put a title on the plot
     set(gca,'fontsize',figfontsize);
     xlabel('Longitude [deg.] ','fontsize',figfontsize);
     ylabel('day','fontsize',figfontsize)
     hcolorbar=colorbar('SouthOutside');
     %caxis([-15,15])
     %set(hcolorbar,'XTick',[-15,-10,-5,0,5,10,15])
     %caxis([-10^(-5),10^(-5)])

    
         set(hcolorbar,'Position',[0.425 0.34 0.18 0.025]);
%        set(gca,'YTick', [365,365+31,365+31+28,365+31+28+31,365+31+28+31+30,365+31+28+31+30+31,365+31+28+31+30+31+30,...
%             365+31+28+31+30+31+30+31,365+31+28+31+30+31+30+31+31,365+31+28+31+30+31+30+31+31+30,365+31+28+31+30+31+30+31+31+30+31,...
%             365+31+28+31+30+31+30+31+31+30+31+30,365+31+28+31+30+31+30+31+31+30+31+30+31]);
%         set(gca,'YTickLabel',{'Jan','Feb','March','April','May','June','July','Aug','Sept','Oct',...
%             'Nov','Dec'});
        set(gca,'XTick',[0,40000*(1/6),40000*(2/6),40000*(3/6),40000*(4/6),40000*(5/6),40000])
        set(gca,'XTickLabel',{'0','60','120','180','240','300','360'});
       
        hold on
      

        
     % Panel label 
     a = text(-5*10^3,365+364,'b)','fontsize',figfontsize,'VerticalAlignment','bottom');


    
    subplot(3,8,[1:2,9:10])
    

    %contour(-N/2:1:N/2-1,freq,log10(24*(30)avext))
    h = pcolor(L*x,T*t/24,Q*qlow_t);
    colormap(jet)
    set(h, 'EdgeColor', 'none');
    %plot(freq,log10(psdxnum),'r*')
    %Compute analytical formula
    % caxis([-10^(-5),10^(-5)])

    grid off
    box on
    
   
     titlestr = strcat('q_{low} []');
     title(titlestr,'fontsize',figfontsize);        %  put a title on the plot    %  put a title on the plot
     set(gca,'fontsize',figfontsize);
     xlabel('Longitude [deg.] ','fontsize',figfontsize);
     ylabel('day','fontsize',figfontsize)
     hcolorbar=colorbar('SouthOutside');
     %caxis([-10,10])
     %set(hcolorbar,'XTick',[-10,-5,0,5,10])

    
         set(hcolorbar,'Position',[0.125 0.34 0.18 0.025]);
          
        set(gca,'XTick',[0,40000*(1/6),40000*(2/6),40000*(3/6),40000*(4/6),40000*(5/6),40000])
        set(gca,'XTickLabel',{'0','60','120','180','240','300','360'});
        
       
        % Add panel label
         a = text(-5*10^3,365+364,'a)','fontsize',figfontsize,'VerticalAlignment','bottom');
         hold on
         
         % Precip
         
          subplot(3,8,[7:8,15:16])
  
       % Define time and q plots for year
       
    %contour(-N/2:1:N/2-1,freq,log10(24*(30)avext))
    h = pcolor(L*x,T*t/24,Q*(Precip_t)/T);
    colormap(jet)
    set(h, 'EdgeColor', 'none');
    %plot(freq,log10(psdxnum),'r*')
    %Compute analytical formula
    % caxis([-10^(-5),10^(-5)])

    grid off
    box on
    
   
     titlestr = strcat('Precip [g/kg h^{-1}]');
     title(titlestr,'fontsize',figfontsize);        %  put a title on the plot    %  put a title on the plot
     set(gca,'fontsize',figfontsize);
     xlabel('Longitude [deg.] ','fontsize',figfontsize);
     ylabel('day','fontsize',figfontsize)
     hcolorbar=colorbar('SouthOutside');
     %caxis([-15,15])
     %set(hcolorbar,'XTick',[-15,-10,-5,0,5,10,15])

    
         set(hcolorbar,'Position',[0.725 0.34 0.18 0.025]);
%        set(gca,'YTick', [365,365+31,365+31+28,365+31+28+31,365+31+28+31+30,365+31+28+31+30+31,365+31+28+31+30+31+30,...
%             365+31+28+31+30+31+30+31,365+31+28+31+30+31+30+31+31,365+31+28+31+30+31+30+31+31+30,365+31+28+31+30+31+30+31+31+30+31,...
%             365+31+28+31+30+31+30+31+31+30+31+30,365+31+28+31+30+31+30+31+31+30+31+30+31]);
%         set(gca,'YTickLabel',{'Jan','Feb','March','April','May','June','July','Aug','Sept','Oct',...
%             'Nov','Dec'});
        set(gca,'XTick',[0,40000*(1/6),40000*(2/6),40000*(3/6),40000*(4/6),40000*(5/6),40000])
        set(gca,'XTickLabel',{'0','60','120','180','240','300','360'});
       
        hold on
      

        
     % Panel label 
     a = text(-5*10^3,365+364,'b)','fontsize',figfontsize,'VerticalAlignment','bottom');

   
      



