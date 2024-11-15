%CCDF A

clc, clear
close all
format compact
addpath(genpath('Data'))
tic

% 3 RX   h = 0
%       LAT °N   LON °E
% low     7.5     81        (Sri Lanka)
% mid    45.0    -63        (Canada)
% high   64.5    -21        (Iceland)

hs = 0.5;     % height above mean sea level of the earth station [km]
lat = [7.5, 45, 64.5]; % [degN]

%% 2.2.1.1  CALCULATION OF LONG-TERM RAIN ATTENUATION STATISTICS FROM POINT RAINFALL RATE
% Carrier frequency up to 55GHz

% Input : Rainfall Rate              [mm/h]
%         Probability of exceedance  [%]
%
%         Height above mean sea level of the earth stations [km]
%         LAT, LON of earth stations [deg]
%         Elevaton angle             [deg]
%         Carrier Frequency          [GHz]
%
% Output: Attenuation exceeded       [dB]

%%  1. Mean annual rain height above mean sea level [km] (P.839)

%           h0 [km]
% low       4.616  
% mid       3.097
% high      1.120

h0 = [4.616, 3.097, 1.120];
hR = h0 + 0.36;  % [km]

%% A(θ)  p = 0.1%  f = 19GHz

%%  2. Slant-path lenght Ls [km] below rain height

th = [2 3 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90];     % [deg] elevation angle
th_ = th';
Re = 8500;  % [km] effective radius of Earth

thr = th*pi/180;

for k = 1:size(th,2)
    
if th(k) >= 5
    
    Ls(k,:) = (hR-hs)./sin(thr(k)); % [km] slant-path length
    
elseif th(k) < 5
    
    Ls(k,:) = 2*(hR-hs)./((sin(thr(k))^2+2*(hR-hs)/Re).^(1/2)+sin(thr(k))); % [km]
    
end

%%  3. Horizontal projection LG of slant-path lenght [km]

LG = Ls(k,:)*cos(thr(k)); % [km]

%%  4. Rainfall rate R001 exceeded for 0.01% of an avg year [mm/h] (P.837)
    
%           R001 [mm/h]
% low       91.209  
% mid       32.480 
% high      20.881

R001 = [91.209, 32.480, 20.881];

for i = 1:3
    
    if R001(i) == 0

        Ap(i) = 0;
        fprintf ('Ap(%d) = 0',i);
        
    end
end
        
%%  5. Specific attenuation γR [dB/km] (P.838)

% O3b mpower operating in Ka band 28Ghz (uplink) circular pol, 19GHz (downlink)

f = 19;   % carrier frequency [Ghz]
tau = 45; % polarization angle τ [deg]
          % 45° for circular polarization

% k & α as function of f = [19 28] [GHz]
% kH, kV, αH, αV constants from P.838 table 5, depending on f

kH = 0.08084;
aH = 1.0691;
kV = 0.08642;
aV = 0.9930;
taur = tau*pi/180;  % [rad]


K = (kH+kV+(kH-kV).*cos(thr(k)).^2*cos(2*taur))/2;
a = (kH.*aH+kV.*aV+(kH.*aH-kV.*aV).*cos(thr(k)).^2.*cos(2*taur))./(2*K);

gR = K.*R001.^a; % [dB/km]

%%  6. Horizontal reduction factor for 0.01% of the time

r001 = 1./(1+0.78*sqrt((LG.*gR)/f)-0.38*(1-exp(-2*LG)));

%%  7. Vertical adjustment factor for 0.01% of the time


zr = atan((hR-hs)./(LG.*r001)); % [rad]
z = zr*180/pi;                  % [deg]

dif(k,:)=z-th(k);


for i = 1:3
    if z(i) > th(k)
    
        LR(1,i) = LG(i).*r001(i)/cos(thr(k));   % [km]
    
    else

        LR(1,i) = (hR(i)-hs)/sin(thr(k));      % [km]
    
    end
end

chi = ones(1,3);

for i = 1:3
    
    if abs(lat(i)) < 36
    
        chi(i) = 36-abs(lat(i)); %[deg]
    
    else
    
        chi(i) = 0;              %[deg]
    
   end    
end

v001 = 1./(1+(sqrt(sin(thr(k))))*(31*(1-exp(-(th(k)./(1+chi)))).*((sqrt(LR.*gR))/(f^2))-0.45));

%%  8. Effective path lenght [km]

LE = LR.*v001;  % [km]

%%  9. Predicted attenuation exceeded for 0.01% of avg year [dB]

A001 = gR.*LE;  % [dB]
A001_(k,:) = A001;

%% 10. Estimated attenuation exceeded in the rage 0.001% to 5%

p = 0.1; % probability of exceedance [%]

%Ap = ones(1,3);

for i = 1:3
    for j = 1:size(p,2)

        
            if p(j) >= 1 || abs(lat(i)) >= 36

                b = 0;

            elseif p(j) < 1 && abs(lat(i)) < 36 && th(k) >= 25

                b = -0.005*(abs(lat(i))-36);

            else

                b = -0.005*(abs(lat(i))-36)+1.8-4.25*sin(thr(k));

            end

            Ap(k,i) = A001(i).*(p(j)/0.01).^(-(0.655+0.033*log(p(j))-0.045*log(A001(i))-b*(1-p(j))*sin(thr(k)))); % [dB] 
        
    end
end

end

figure
plot(th,Ap(:,1),th,Ap(:,2),th,Ap(:,3),'LineWidth',2)
title(['A(θ)  p = 0.1%   Ka band downlink frequency (19GHz)'])
xlabel('θ[deg]')
ylabel('A[dB]')

lgd = legend('low (7.5°N)','mid (45°N)','high (64.5°N)');
title(lgd,['Ground stations latitudes'])
hold on
grid on
%zoom('xon')

%% P(A) th = 3   20  45 deg   +P0    f = 19GHz 

%%  2. Slant-path lenght Ls [km] below rain height

lat = [7.5, 45, 64.5]; % [degN]

figure 
th = [5 20 45 90];     % [deg] elevation angle
Re = 8500;  % [km] effective radius of Earth

thr = th*pi/180;

for k = 1:size(th,2)
    
if th(k) >= 5
    
    Ls(k,:) = (hR-hs)./sin(thr(k)); % [km] slant-path length
    
elseif th(k) < 5
    
    Ls(k,:) = 2*(hR-hs)./((sin(thr(k))^2+2*(hR-hs)/Re).^(1/2)+sin(thr(k))); % [km]
    
end

%%  3. Horizontal projection LG of slant-path lenght [km]

LG = Ls(k,:)*cos(thr(k)); % [km]

%%  4. Rainfall rate R001 exceeded for 0.01% of an avg year [mm/h] (P.837)
    
%           R001 [mm/h]
% low       91.209  
% mid       32.480 
% high      20.881

R001 = [91.209, 32.480, 20.881];

for i = 1:3
    
    if R001(i) == 0

        Ap(i) = 0;
        fprintf ('Ap(%d) = 0',i);
        
    end
end
        
%%  5. Specific attenuation γR [dB/km] (P.838)

% O3b mpower operating in Ka band 28Ghz (uplink) circular pol, 19GHz (downlink)

f = 19;   % carrier frequency [Ghz]
tau = 45; % polarization angle τ [deg]
          % 45° for circular polarization

% k & α as function of f = [19 28] [GHz]
% kH, kV, αH, αV constants from P.838 table 5, depending on f

kH = 0.08084;
aH = 1.0691;
kV = 0.08642;
aV = 0.9930;
taur = tau*pi/180;  % [rad]


K = (kH+kV+(kH-kV).*cos(thr(k)).^2*cos(2*taur))/2;
a = (kH.*aH+kV.*aV+(kH.*aH-kV.*aV).*cos(thr(k)).^2.*cos(2*taur))./(2*K);

gR = K.*R001.^a; % [dB/km]

%%  6. Horizontal reduction factor for 0.01% of the time

r001 = (1+0.78*sqrt((LG.*gR)/f)-0.38*(1-exp(-2*LG))).^-1;

%%  7. Vertical adjustment factor for 0.01% of the time


zr = atan((hR-hs)./(LG.*r001)); % [rad]
z = zr*180/pi;                  % [deg]

dif(k,:)=z-th(k);


for i = 1:3
    if z(i) > th(k)
    
        LR(1,i) = LG(i).*r001(i)/cos(thr(k));   % [km]
    
    else

        LR(1,i) = (hR(i)-hs)/sin(thr(k));      % [km]
    
    end
end

chi = ones(1,3);

for i = 1:3
    
    if abs(lat(i)) < 36
    
        chi(i) = 36-abs(lat(i)); %[deg]
    
    else
    
        chi(i) = 0;              %[deg]
    
    end    
end

v001 = 1./(1+(sqrt(sin(thr(k))))*(31*(1-exp(-(th(k)./(1+chi)))).*((sqrt(LR.*gR))/(f^2))-0.45));

%%  8. Effective path lenght [km]

LE = LR.*v001;  % [km]

%%  9. Predicted attenuation exceeded for 0.01% of avg year [dB]

A001 = gR.*LE;  % [dB]

%% 10. Estimated attenuation exceeded in the rage 0.001% to 5%

p = [0.001 0.01 0.1 0.2 0.3 0.5 1 2 3 5]; % probability of exceedance [%]

Ap = ones(1,3);

for i = 1:3
    for j = 1:size(p,2)

        
            if p(j) >= 1 || abs(lat(i)) >= 36

                b = 0;

            elseif p(j) < 1 && abs(lat(i)) < 36 && th(k) >= 25

                b = -0.005*(abs(lat(i))-36);

            else

                b = -0.005*(abs(lat(i))-36)+1.8-4.25*sin(thr(k));

            end

            Ap(j,i) = A001(i)*(p(j)/0.01)^(-(0.655+0.033*log(p(j))-0.045*log(A001(i))-b*(1-p(j))*sin(thr(k)))); % [dB] 
            
            if k == 1
                 Ap_(j,i) = Ap(j,i);
            end
    end
end

subplot(2,2,k)
semilogy(Ap(:,1),p,Ap(:,2),p,Ap(:,3),p,'LineWidth',2)
title(['P(A)   θ = ',num2str(th(k)),'°   Ka band downlink frequency (19GHz)'])
xlabel('A[dB]')
ylabel('p[%]')

lgd = legend('low (7.5°N)','mid (45°N)','high (64.5°N)');
title(lgd,['Ground stations latitudes'])
hold on
grid on
zoom('xon')

end

%% P(A>0) for P(A) plot

%% 2.2.1.2 Probability of rain attenutation on a slant path

% Input : LAT, LON of earth stations [degN, degE]
%         Elevation angle            [deg]
%         Ls slant path length from earth station to rain height [km]
%
% Output: P(A > 0)                   [%]

%% 0. Collecting Tii and MTii data and interpolation (P.1510, P.837, P.1144)

% Bi-linear interpolation on square grid (P.1144 annex 1b) for MTii
% MT [mm] monthly mean total rainfall (P.837)
% T  [K]  monthly mean surface temperatures (P.1510)
% N       number of days in each month
%
% MTii(r,c) = (MTii(R,C)+MTii(R+1,C)+MTii(R,C+1)+MTii(R+1,C+1))*((0.5)^2);

delimiterIn = ' ';
filename = 'MT_Month01.TXT';
MT1 = importdata(filename,delimiterIn);
filename = 'MT_Month02.TXT';
MT2 = importdata(filename,delimiterIn);
filename = 'MT_Month03.TXT';
MT3 = importdata(filename,delimiterIn);
filename = 'MT_Month04.TXT';
MT4 = importdata(filename,delimiterIn);
filename = 'MT_Month05.TXT';
MT5 = importdata(filename,delimiterIn);
filename = 'MT_Month06.TXT';
MT6 = importdata(filename,delimiterIn);
filename = 'MT_Month07.TXT';
MT7 = importdata(filename,delimiterIn);
filename = 'MT_Month08.TXT';
MT8 = importdata(filename,delimiterIn);
filename = 'MT_Month09.TXT';
MT9 = importdata(filename,delimiterIn);
filename = 'MT_Month10.TXT';
MT10 = importdata(filename,delimiterIn);
filename = 'MT_Month11.TXT';
MT11 = importdata(filename,delimiterIn);
filename = 'MT_Month12.TXT';
MT12 = importdata(filename,delimiterIn);

R = [391, 541, 619]; % Rows and columns at RX Lat and Lon on file grid
C = [1045, 469, 637];

MTA = ones(12,3);
MTB = ones(12,3);
MTC = ones(12,3);
MTD = ones(12,3);

for i = 1:3

MTA(:,i) = [MT1(R(i),C(i)), MT2(R(i),C(i)), MT3(R(i),C(i)), MT4(R(i),C(i)), MT5(R(i),C(i)), MT6(R(i),C(i)), MT7(R(i),C(i)), MT8(R(i),C(i)), MT9(R(i),C(i)), MT10(R(i),C(i)), MT11(R(i),C(i)), MT12(R(i),C(i))];
MTB(:,i) = [MT1(R(i),C(i)+1), MT2(R(i),C(i)+1), MT3(R(i),C(i)+1), MT4(R(i),C(i)+1), MT5(R(i),C(i)+1), MT6(R(i),C(i)+1), MT7(R(i),C(i)+1), MT8(R(i),C(i)+1), MT9(R(i),C(i)+1), MT10(R(i),C(i)+1), MT11(R(i),C(i)+1), MT12(R(i),C(i)+1)];
MTC(:,i) = [MT1(R(i)+1,C(i)), MT2(R(i)+1,C(i)), MT3(R(i)+1,C(i)), MT4(R(i)+1,C(i)), MT5(R(i)+1,C(i)), MT6(R(i)+1,C(i)), MT7(R(i)+1,C(i)), MT8(R(i)+1,C(i)), MT9(R(i)+1,C(i)), MT10(R(i)+1,C(i)), MT11(R(i)+1,C(i)), MT12(R(i)+1,C(i))];
MTD(:,i) = [MT1(R(i)+1,C(i)+1), MT2(R(i)+1,C(i)+1), MT3(R(i)+1,C(i)+1), MT4(R(i)+1,C(i)+1), MT5(R(i)+1,C(i)+1), MT6(R(i)+1,C(i)+1), MT7(R(i)+1,C(i)+1), MT8(R(i)+1,C(i)+1), MT9(R(i)+1,C(i)+1), MT10(R(i)+1,C(i)+1), MT11(R(i)+1,C(i)+1), MT12(R(i)+1,C(i)+1)];

end

MT = ones(12,3);

for i = 1:3
    for j = 1:12

        MT(j,i) = (MTA(j,i)+MTB(j,i)+MTC(j,i)+MTD(j,i)).*((0.5).^2);

    end
end


filename = 'T_Month01.TXT';
T1 = importdata(filename,delimiterIn);
filename = 'T_Month02.TXT';
T2 = importdata(filename,delimiterIn);
filename = 'T_Month03.TXT';
T3 = importdata(filename,delimiterIn);
filename = 'T_Month04.TXT';
T4 = importdata(filename,delimiterIn);
filename = 'T_Month05.TXT';
T5 = importdata(filename,delimiterIn);
filename = 'T_Month06.TXT';
T6 = importdata(filename,delimiterIn);
filename = 'T_Month07.TXT';
T7 = importdata(filename,delimiterIn);
filename = 'T_Month08.TXT';
T8 = importdata(filename,delimiterIn);
filename = 'T_Month09.TXT';
T9 = importdata(filename,delimiterIn);
filename = 'T_Month10.TXT';
T10 = importdata(filename,delimiterIn);
filename = 'T_Month11.TXT';
T11 = importdata(filename,delimiterIn);
filename = 'T_Month12.TXT';
T12 = importdata(filename,delimiterIn);

R = [131, 181, 207];
C = [349, 157, 213];

T = ones(12,3);

for i = 1:3

    T(:,i) = [T1(R(i),C(i)), T2(R(i),C(i)), T3(R(i),C(i)), T4(R(i),C(i)), T5(R(i),C(i)), T6(R(i),C(i)), T7(R(i),C(i)), T8(R(i),C(i)), T9(R(i),C(i)), T10(R(i),C(i)), T11(R(i),C(i)), T12(R(i),C(i))];

end

N = [31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

%% 1. Annual probability of rain [%] (P.837)

t = T - 273.15; % [deg C]
r = ones(12,3);
P0 = ones(12,3);
P0a = ones(12,3);

for i = 1:3
    for j = 1:12
       
        if t(j,i) >= 0 
     
            r(j,i) = 0.5874*exp(0.0883*t(j,i)); % [mm/h] 
            
        elseif t(j,i) < 0 
        
            r(j,i) = 0.5874;
            
        end
        
        P0(j,i) = 100*MT(j,i)/(24*N(j)*r(j,i)); % [%] monthly probability of rain
        
        if P0(j,i) > 70
            
            P0(j,i) = 70;
            r(j,i) = (100/70)*(MT(j,i)/24*N(j));
            
        end
        
        P0a(j,i) = (N(j)*P0(j,i))/365.25;
    end
end

P0annual = sum(P0a); % [%] P(R>0) annual probability of rain at earth stations

%% 2. Parameter α

% Complementary cumulative normal distribution function (CCDF):
% Q(x) = (1/sqrt(2*pi))*int(e^(-(t^2)/2)dt) (x->+inf)
% Q(x) = (erfc(x/sqrt(2)))/2 
% Inverse complementary cumulative normal distribution function:
% Q(x)^-1 = sqrt(2)*erfc^-1(2*x) 
% erfc(x): complementary error function

x = P0annual/100;
alpha = sqrt(2)*erfcinv(2*x); %Q^-1

%% 3. Spatial correlation function ρ
% -1 < ρ < +1

th = [3 20 45 90]; % [deg]
thr = th*pi/180; % [rad]

for i = 1:size(th,2)
    
    d = Ls(i,:).*cos(thr(i)); % [km]
    rho = 0.59*exp(-abs(d)./31)+0.41*exp(-abs(d)/800); % spatial correlation function

%% 4. Complementary bivariate normal distribution cB
% Gaussian quadrature based on K points
% Z.Drezner and G.O.Wesolowsky method K=2,3 (approximation)

rho2(1,:) = (3-sqrt(3))*rho/6;
rho2(2,:) = (3+sqrt(3))*rho/6;

rho3(1,:) = (1-sqrt(3/5))*(rho/2);
rho3(2,:) = rho/2;
rho3(3,:) = (1+sqrt(3/5))*(rho/2);

cB2 = rho/(2*2*pi).*((1./(sqrt(1-rho2(1,:).^2)).*exp(-(alpha.^2)./(1+rho2(1,:))))+(1./(sqrt(1-rho2(2,:).^2)).*exp(-(alpha.^2)./(1+rho2(2,:)))))+x.^2;
cB3 = rho/(18*2*pi).*((5./(sqrt(1-rho3(1,:).^2)).*exp(-(alpha.^2)./(1+rho3(1,:))))+(8./(sqrt(1-rho3(2,:).^2)).*exp(-(alpha.^2)./(1+rho3(2,:))))+(5./(sqrt(1-rho3(3,:).^2)).*exp(-(alpha.^2)./(1+rho3(3,:)))))+x.^2;

%% 5. Probability of rain attenutation on a slant path P(A>0) [%]

%PA0(1,:) = 100*(1-(1-x).*((cB2-x.^2)./(x.*(1-x))).^x); % [%]
PA0(2,:) = 100*(1-(1-x).*((cB3-x.^2)./(x.*(1-x))).^x); % [%]

if th(i) == 90
    
    PA0(2,:) = P0annual; % rho = 1 is unacceptable for cB definition
    
end

PA0_(i,:) = PA0(2,:);

end

for i = 1:size(th,2)

    subplot(2,2,i)
    semilogy([0],PA0_(i,1),[0],PA0_(i,2),[0],PA0_(i,3),'Marker','o','MarkerFaceColor','auto')
   lgd = legend('low (7.5°N)','mid (45°N)','high (64.5°N)','P(A>0) low','P(A>0) mid','P(A>0) high');
    title(lgd,['Ground stations latitudes']),
    
end



%% frequency scaling for A  θ = 5°

%% 2.2.1.3 Long-term frequency and polarization scaling of rain attenuation statistics
% Prediction of propagation effects at one frequency, form knowledge of the
% propagation effect at a different frequency (typically, the predicted f is higher)

%% 2.2.1.3.2 Long-term frequency scaling of rain attenuation statistics
% Simplified method, without cumulative distribution of rain attenuation.
% Predicts the uplink rain attenuation at f2, knowing downlink
% rain attenuation at f1 at the same probability of exceedance.

f1 = 19;  % [GHz]
f2 = 28;  % [GHz]

Ap1 = Ap_;

phi1 = (f1^2)/(1+10^-4*f1^2);
phi2 = (f2^2)/(1+10^-4*f2^2);

H = 1.12*10^-3*(phi2./phi1)^0.5.*(phi1.*Ap1).^0.55;

Ap2 = Ap1.*(phi2./phi1).^(1-H); % [dB]

p = [0.001 0.01 0.1 0.2 0.3 0.5 1 2 3 5]; % [%]
figure
plot(Ap1(:,2),p, Ap2(:,2),p,'Linewidth',2)
title('P(A)      θ = 5°')
xlabel('A[dB]')
ylabel('P(A)[%]')
lgd = legend(['Ka band downlink frequency (',num2str(f1),'GHz)'],['Ka band uplink frequency (',num2str(f2),'GHz)']);

title(lgd,'Ground station at medium latitude (45°N)')
grid on

%% Worst month Pw(A)

%% 2.2.2 Seasonal variations - worst month (P.841)
% Attenuation value exceeded for specified time percentage pw
% of the worst month
% Input : pw [%] average annual worst-month time % of excess
%        (β,Q1) coefficients of Q, from table 1, for each receiver location
%
%
% Output: A [dB] attenuation exceeded for the annual time percentage P
%         

%% 1.  Annual time percentage P
% Q(Q1,β) is a function of P(average annual time percentage of excess)

beta = [0.22, 0.1, 0.16]; % coefficients of Q(P)
Q1 = [1.7, 2.7, 3.8];

pwmin = 12*(Q1./12).^(1./beta);    % pw limits
pwmax = Q1.*3.^(1-beta);

pw1(:,1) = [0.005 0.01 0.1 0.2 0.3 0.5 1 2 3 4]; % [%]
pw2(:,1) = [0.001 0.01 0.1 0.2 0.3 0.5 1 2 3 4 5 6 7];
pw3(:,1) = [0.01 0.1 0.2 0.3 0.5 1 2 3 4 5 6 7 8 9];

Q = Q1.^(1./(1+beta)).*pw1.^(-beta./(1-beta));
P1 = pw1./Q(:,1); 

Q = Q1.^(1./(1+beta)).*pw2.^(-beta./(1-beta));
P2 = pw2./Q(:,2); 

Q = Q1.^(1./(1+beta)).*pw3.^(-beta./(1-beta));
P3 = pw3./Q(:,3); 

%% 2. Attenuation exceeded for the annual time percentage P
% from 2.2.1.1 step 10.
% A is estimated attenuation also for pw% of time of the worst-month

f = 19; % [GHz]
th = 5; % [deg]
thr = th*pi/180; % [rad]
A001 = A001_(11,:);
lat = [7.5, 45, 64.5]; % [degN]

clear A1 A2 A3
for i = 1:size(P1,1)
    
    if P1(i) >= 1 || abs(lat(1)) >= 36
        
        b = 0;
        
    elseif P1(i) < 1 && abs(lat(1)) < 36 && th >= 25
       
        b = -0.005*(abs(lat(1))-36);
 
    else
        
        b = -0.005*(abs(lat(1))-36)+1.8-4.25*sin(thr);
        
    end
    
    A1(i,1) = A001(1)*(P1(i)/0.01)^(-(0.655+0.033*log(P1(i))-0.045*log(A001(1))-b*(1-P1(i))*sin(thr))); % [dB]  
    
end

for i = 1:size(P2,1)
    
    if P2(i) >= 1 || abs(lat(2)) >= 36
        
        b = 0;
        
    elseif P2(i) < 1 && abs(lat(2)) < 36 && th >= 25
       
        b = -0.005*(abs(lat(2))-36);
 
    else
        
        b = -0.005*(abs(lat(2))-36)+1.8-4.25*sin(thr);
        
    end
    
    A2(i,1) = A001(2)*(P2(i)/0.01)^(-(0.655+0.033*log(P2(i))-0.045*log(A001(2))-b*(1-P2(i))*sin(thr))); % [dB]  
    
end

for i = 1:size(P3,1)
    
    if P3(i) >= 1 || abs(lat(3)) >= 36
        
        b = 0;
        
    elseif P3(i) < 1 && abs(lat(3)) < 36 && th >= 25
       
        b = -0.005*(abs(lat(3))-36);
 
    else
        
        b = -0.005*(abs(lat(3))-36)+1.8-4.25*sin(thr);
        
    end
    
    A3(i,1) = A001(3)*(P3(i)/0.01)^(-(0.655+0.033*log(P3(i))-0.045*log(A001(3))-b*(1-P3(i))*sin(thr))); % [dB]  
    
end

figure
%subplot(1,2,1)
%semilogy(A1(:,1),pw1,A2(:,1),pw2,A3(:,1),pw3,'Linewidth',2)
%grid on
%subplot(1,2,2)
semilogy(A1(:,1),P1,A2(:,1),P2,A3(:,1),P3,'Linewidth',2)
title(['P(A) worst-month      θ = ',num2str(th),'°     Ka band downlink frequency (',num2str(f) 'GHz)'])
xlabel('A[dB]')
ylabel('p[%]')

lgd = legend('low (7.5°N)','mid (45°N)','high (64.5°N)');
title(lgd,'Earth station latitudes')
grid on

%% 2.2.4 Site diversity
% At frequencies above  20 Ghz, path impariments other than rain affect site
% diversity performance


%% 2.2.4.1 Prediction of outage probability due to rain att. with site diversity
% For unbalanced and balanced systems
% Most accurate, preferred
% Pr: joint probability that it is raining at both sites 
% Pa: conditional j.p. attenuations exceed a1 and a2, given that it
%     is raining at both sites
% R1,R2: thresholds
% Input : Pk: rain for a particular location (P.837) 2 sites [%]
%         Pi: probability the attenuation Ai is exceeded, for each Rx
%         d: separation between 2 sites [km]
%         a1, a2 [dB]
%         A001: attenuation exceeded at 0.01% for each Rx [dB]
%         thr: elevation angle [rad] 
%
% Output: PR(A1>=a1,A2>=a2): outage probability of rain[%]

d = 30; % [km]
a1 = 20; % [dB]
a2 = 35; % [dB]
th = 45; % [deg]
%th = th_(11);
thr = th*pi/180; % [rad]
A001 = A001_(11,:);
lat = [7.5, 45, 64.5]; % [degN]

%% 1. Determine Pk rain probability on the k-th path (P.837)
% Chosen from local data or from rainfall rate maps R001
% Studied for a d separation from each RX of this project,
% at constant A001 between the 2 sites 

Pk = [15 20]; % [%]

%% 2. Sets of pairs [Pi,Ai] (Pi<=Pk) 

%Pi = ones(13,3);

for i=1:3

    Pi(:,i) = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10];

end

n = size(Pi,1);

% from 2.2.1.1 step 10.

%b = ones(n,3);
%Ai = ones(n,3);

for i = 1:3
    for j = 1:n
     
        if Pi(j,i) >= 1 || abs(lat(i)) >= 36
        
            b(j,i) = 0;
        
        elseif Pi(j,i) < 1 && abs(lat(i)) < 36 && th >= 25
       
            b(j,i) = -0.005*(abs(lat(i))-36);
 
        else
        
            b(j,i) = -0.005*(abs(lat(i))-36)+1.8-4.25*sin(thr);
        
        end
    
        Ai(j,i) = A001(i)*(Pi(j,i)/0.01)^(-(0.655+0.033*log(Pi(j,i))-0.045*log(A001(i))-b(j,i)*(1-Pi(j,i))*sin(thr))); % [dB]  
    
    end
end

%% 3. Transform set of pairs
% Pi->Q^-1(Pi/Pk) inverse complementary cumulative normal distribution
% Ai->ln(Ai)

Zi_1 = sqrt(2)*(erfcinv(2*Pi/Pk(1)));
Zi_2 = sqrt(2)*(erfcinv(2*Pi/Pk(2)));

lnAi = log(Ai); 

%% 4. Determine variables mlnA and σlnA: (P.1057 annex2) Least-square
% Log-normal distribution fit of rain att. vs probability of occurrence

szln_1 = ones(n,3);
sz2_1 = ones(n,3);
szln_2 = ones(n,3);
sz2_2 = ones(n,3);

for i = 1:3
    for j = 1:n
        
        szln_1(j,i) = Zi_1(j,i)*lnAi(j,i);
        sz2_1(j,i) = Zi_1(j,i)^2;
        szln_2(j,i) = Zi_2(j,i)*lnAi(j,i);
        sz2_2(j,i) = Zi_2(j,i)^2;
        
    end
end

sszln_1 = sum(szln_1);
ssz_1 = sum(Zi_1);
ssln = sum(lnAi);
ssz2_1 = sum(sz2_1);

sszln_2 = sum(szln_2);
ssz_2 = sum(Zi_2);
ssz2_2 = sum(sz2_2);

sigma(1,:) = (n*sszln_1-ssz_1.*ssln)./(n*ssz2_1-ssz_1.^2);
sigma(2,:) = (n*sszln_2-ssz_2.*ssln)./(n*ssz2_2-ssz_2.^2);

m(1,:) = (ssln-sigma(1,:).*ssz_1)/n;
m(2,:) = (ssln-sigma(2,:).*ssz_2)/n;


%% 5.

R1 = sqrt(2)*(erfcinv(2*Pk(1)/100));
R2 = sqrt(2)*(erfcinv(2*Pk(2)/100));

rhor = 0.7*exp(-d/60)+0.3*exp(-(d/700).^2);
rhor3(1,:) = 1-(rhor/2)*sqrt(3/5);
rhor3(2,:) = rhor/2;
rhor3(3,:) = 1+(rhor/2)*sqrt(3/5);
                                                                                                                                                           
Pr = rhor/(18*2*pi).*((5./(sqrt(1-rhor3(1,:).^2)).*exp(-(R1.^2+R2.^2-2*R1.*R2.*rhor3(1,:))./(2*(1+rhor3(1,:).^2))))+(8./(sqrt(1-rhor3(2,:).^2)).*exp(-(R1.^2+R2.^2-2*R1.*R2.*rhor3(2,:))./(2*(1+rhor3(2,:).^2))))+(5./(sqrt(1-rhor3(3,:).^2)).*exp(-(R1.^2+R2.^2-2*R1.*R2.*rhor3(3,:))./(2*(1+rhor3(3,:).^2)))))+Pk(1).*Pk(2)/10000;

h1 = (log(a1)-m(1,:))./sigma(1,:);
h2 = (log(a2)-m(2,:))./sigma(2,:);

rhoa = 0.94*exp(-d/30)+0.06*exp(-(d/500).^2);
rhoa3(1,:) = 1-(rhoa/2)*sqrt(3/5);
rhoa3(2,:) = rhoa/2;
rhoa3(3,:) = 1+(rhoa/2)*sqrt(3/5);

Pa = rhoa/(18*2*pi).*((5./(sqrt(1-rhoa3(1,:).^2)).*exp(-(h1.^2+h2.^2-2*h1.*h2.*rhoa3(1,:))./(2*(1+rhoa3(1,:).^2))))+(8./(sqrt(1-rhoa3(2,:).^2)).*exp(-(h1.^2+h2.^2-2*h1.*h2.*rhoa3(2,:))./(2*(1+rhoa3(2,:).^2))))+(5./(sqrt(1-rhoa3(3,:).^2)).*exp(-(h1.^2+h2.^2-2*h1.*h2.*rhoa3(3,:))./(2*(1+rhoa3(3,:).^2)))))+erfc(h1/sqrt(2)).*erfc(h2/sqrt(2))/4;              
                                                                                                                                                                                                                                                                                                                                              
PR = 100*Pr.*Pa; % [%] PR(A1>=a1,A2>=a2) 

%% G(A)

%% 2.2.4.2 Diversity gain
% For balanced systems with short distances
% Simplified, for separation distances less than 20 km
% Less accurate
% Input : d: separation between 2 sites [km]
%         A: path rain attenuation for single site [dB]   
%         f [GHz]
%         th [deg]
%         ψ: azimuth angle of the prop path wrt baseline between sites,
%            such that ψ<=90° [deg]
%
% Output: G: net diversity gain [dB]

d = [5 15]; % [km]
A = [1 2 3 5 10 15 20 25 30 35 40 45 50]'; % [dB]
f = [19 28]; % [GHz]
th_ = [2 3 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90];     % [deg] elevation angle
th = [5 20 45]; % [deg]
psi = [0 10 30 45 60 90]; % [deg]

%fixed values:
d_f = d(2);
f_f = f(2);
th_f = th_(11);
psi_f = psi(4);

%% 1. Gain contributed by spatial separation

a = 0.78*A-1.94*(1-exp(-0.11*A));
b = 0.59*(1-exp(-0.1*A));
clear Gd
for i=1:size(A,1)
    
    Gd(i,:) = a(i).*(1-exp(-b(i).*d));
    Gd_f(i,1) = a(i).*(1-exp(-b(i).*d_f));

end

%% 2. Frequency-dependent gain 

Gf = exp(-0.025*f);
Gf_f = exp(-0.025*f_f);

%% 3. Gains term dependent on elevation angle

Gth = 1+0.006*th;
Gth_f = 1+0.006*th_f;

%% 4. Baseline-dependent term

Gpsi = 1+0.002*psi;
Gpsi_f = 1+0.002*psi_f;

%% Net diversity gain

G1 = Gd.*Gf_f.*Gth_f.*Gpsi_f; %[dB]
clear G2 G3 G4
for i = 1:size(A,1)
G2(i,:) = Gd_f(i).*Gf.*Gth_f.*Gpsi_f; %[dB]
G3(i,:) = Gd_f(i).*Gf_f.*Gth.*Gpsi_f; %[dB]
G4(i,:) = Gd_f(i).*Gf_f.*Gth_f.*Gpsi; %[dB]
end

figure

subplot(2,2,1)
plot(A,G1(:,1),A,G1(:,2),'Linewidth',2)
grid on
title('G(A)  separation distance')
xlabel('A[dB]')
ylabel('G[dB]')
legend(['d = ',num2str(d(1)),'km'],['d = ',num2str(d(2)),'km'])
    
subplot(2,2,2)
plot(A,G2(:,1),A,G2(:,2),'Linewidth',2)
grid on
title('G(A)  frequency')
xlabel('A[dB]')
ylabel('G[dB]')
legend(['f = ', num2str(f(1)),'GHz'],['f = ', num2str(f(2)),'GHz'])
hold on

subplot(2,2,3)
plot(A,G3(:,1),A,G3(:,2),A,G3(:,3),'Linewidth',2)
grid on
title('G(A)  elevation')
xlabel('A[dB]')
ylabel('G[dB]') 
legend(['θ = ', num2str(th(1)),'°'],['θ = ', num2str(th(2)),'°'],['θ = ', num2str(th(3)),'°'])
hold on
subplot(2,2,4)
plot(A,G4(:,1),A,G4(:,2),A,G4(:,3),A,G4(:,4),A,G4(:,5),A,G4(:,6),'Linewidth',2)
grid on
title('G(A)  azimuth')
xlabel('A[dB]')
ylabel('G[dB]')
legend(['ψ = ', num2str(psi(1)),'°'],['ψ = ', num2str(psi(2)),'°'],['ψ = ', num2str(psi(3)),'°'],['ψ = ', num2str(psi(4)),'°'],['ψ = ', num2str(psi(5)),'°'],['ψ = ', num2str(psi(6)),'°'])
hold on



%% 2.4 Scintillation and multipath fading

%% 2.4.1 Amplitude scintillation fading at free space, θ >= 5°
% th >= 5°
% Inputs : Nwet: (avg y) obtained from (P.453)
%                average year exceedance digital maps 
%          η: eta antenna efficiency
%          D: geometrical antenna diameter [m] 
%          th: desired apparent elevation angle [deg]
%          p(%):      0.01% < p <= 50%    
%          f: [GHz] carrier frequency
%
% Output : Afd: A(p) Fade dexceededepth  for p(%) of time 

eta = 0.5;
D = 4.5; % [m] Viasat MEOlink
th0 = [5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90]; % [deg]
th9 = th0;
th0r = th0*pi/180; % [rad]
p = 0.1; % [%] time percentage
%p = [0.1; 0.2; 0.3; 0.5; 1; 2; 3; 5; 10; 20; 30; 50]; % [%] time percentage
f = [19 28]'; % [GHz] downlink and uplink Viasat MEOLink

delimiterIn = ' ';
filename = 'NWET_Annual_01.TXT';
TNwet1 = importdata(filename,delimiterIn);
filename = 'NWET_Annual_02.TXT';
TNwet2 = importdata(filename,delimiterIn);
filename = 'NWET_Annual_03.TXT';
TNwet3 = importdata(filename,delimiterIn);
filename = 'NWET_Annual_05.TXT';
TNwet4 = importdata(filename,delimiterIn);
filename = 'NWET_Annual_1.TXT';
TNwet5 = importdata(filename,delimiterIn);
filename = 'NWET_Annual_2.TXT';
TNwet6 = importdata(filename,delimiterIn);
filename = 'NWET_Annual_3.TXT';
TNwet7 = importdata(filename,delimiterIn);
filename = 'NWET_Annual_5.TXT';
TNwet8 = importdata(filename,delimiterIn);
filename = 'NWET_Annual_10.TXT';
TNwet9 = importdata(filename,delimiterIn);
filename = 'NWET_Annual_20.TXT';
TNwet10 = importdata(filename,delimiterIn);
filename = 'NWET_Annual_30.TXT';
TNwet11 = importdata(filename,delimiterIn);
filename = 'NWET_Annual_50.TXT';
TNwet12 = importdata(filename,delimiterIn);
filename = 'NWET_Annual_60.TXT';
TNwet13 = importdata(filename,delimiterIn);
filename = 'NWET_Annual_70.TXT';
TNwet14 = importdata(filename,delimiterIn);
filename = 'NWET_Annual_80.TXT';
TNwet15 = importdata(filename,delimiterIn);
filename = 'NWET_Annual_90.TXT';
TNwet16 = importdata(filename,delimiterIn);
filename = 'NWET_Annual_95.TXT';
TNwet17 = importdata(filename,delimiterIn);
filename = 'NWET_Annual_99.TXT';
TNwet18 = importdata(filename,delimiterIn);

R = [131, 181, 207];
C = [349, 157, 213];

Nwet = ones(18,3);

for i = 1:3

    Nwet(:,i) = [TNwet1(R(i),C(i)), TNwet2(R(i),C(i)), TNwet3(R(i),C(i)), TNwet4(R(i),C(i)), TNwet5(R(i),C(i)), TNwet6(R(i),C(i)), TNwet7(R(i),C(i)), TNwet8(R(i),C(i)), TNwet9(R(i),C(i)), TNwet10(R(i),C(i)), TNwet11(R(i),C(i)), TNwet12(R(i),C(i)), TNwet13(R(i),C(i)), TNwet14(R(i),C(i)), TNwet15(R(i),C(i)), TNwet16(R(i),C(i)), TNwet17(R(i),C(i)), TNwet18(R(i),C(i))];

end

sigma_ref = 3.6*10^-3+Nwet(1:12,:)*10^-4; % std deviation of the reference signal amplitude
hL = 1000; % [m] height of turbulent layer
Deff = sqrt(eta)*D; % [m] effective antenna diameter

%% LAT MID(FREQ)
%% A(θ)  p = 0.1%  rx mid θ > 5
 
for i = 1:size(th0,2)
    
    L = 2*hL./(sqrt(sin(th0r(i)).^2+2.35*10^-4)+sin(th0r(i))); % [m] effective path lenght
    x = 1.22*Deff.^2*(f./L);
    
    for j = 1:size(f,1)
        if x(j,1) >= 7.0 % argument of the square root is negative 
    
             Afd1(j,i) = 0; 
            %fprintf('argument of the square root is negative')
    
        else

        g(j,1) = sqrt(3.86*(x(j,1).^2+1).^(11/12)*sin(11/6*atan(1./x(j,1)))-7.08*x(j,1).^(5/6)); % antenna averaging factor
        sigma(j,1) = sigma_ref(1,2).*f(j,1)^(7/12).*g(j,1)/(sin(th0r(i))^1.2); % std deviation of the signal
        
        ap = -0.061*(log10(p).^3)+0.072*(log10(p).^2)-1.71*log10(p)+3.0; % time percentage factor 

        Afd1(j,i) = ap.*sigma(j,1); % [dB]fade depth exceeded for p(%) of the time
        
        end
    end
end

% 2.4.3  θ = 5° rx mid
    Afd_ = Afd1;
%

figure
plot(th0,Afd1(1,:),th0,Afd1(2,:),'Linewidth',2) 
title(['A(θ)  θ > 5°   p = ',num2str(p),'%'])
xlabel('θ[deg]')
ylabel('A[dB]')
lgd = legend(['f = ',num2str(f(1)),'GHz'],['f = ',num2str(f(2)),'GHz']);
title(lgd,'Rx mid lat (45°N)')
hold on
grid on
zoom('xon')


%% LATs  f = 19GHz
%% A(θ)  p = 0.1%  rx mid θ > 5
 
f = 19;% [GHz]
p = 0.1; % [%]
for i = 1:size(th0,2)
    
    L = 2*hL./(sqrt(sin(th0r(i)).^2+2.35*10^-4)+sin(th0r(i))); % [m] effective path lenght
    x = 1.22*Deff.^2*(f./L);
    
    for k = 1:size(sigma_ref,2)
        if x >= 7.0 % argument of the square root is negative 
    
             Afd0(i,k) = 0; 
            %fprintf('argument of the square root is negative')
    
        else

        g = sqrt(3.86*(x.^2+1).^(11/12)*sin(11/6*atan(1./x))-7.08*x.^(5/6)); % antenna averaging factor
        
            
        
        sigma(1,k) = sigma_ref(1,k).*f^(7/12).*g/(sin(th0r(i))^1.2); % std deviation of the signal
        ap = -0.061*(log10(p).^3)+0.072*(log10(p).^2)-1.71*log10(p)+3.0; % time percentage factor 

        Afd0(i,k) = ap.*sigma(1,k); % [dB]fade depth exceeded for p(%) of the time
        
        
        end
    end
end

% 2.4.3  θ = 5° rx mid
    Afd0_ = Afd0;
%



%% P(A)  θ = 5,20,45,90°  f = 19GHz

eta = 0.5;
D = 4.5; % [m] Viasat MEOlink
th0 = [5 20 45 90]; % [deg]
th0r = th0*pi/180; % [rad]
p = [0.1; 0.2; 0.3; 0.5; 1; 2; 3; 5; 10; 20; 30; 50]; % [%] time percentage
f = 19; % [GHz] downlink and uplink Viasat MEOLink
hL = 1000; % [m] height of turbulent layer
Deff = sqrt(eta)*D; % [m] effective antenna diameter
figure

for i = 1:size(th0,2)
    for j = 1:3
        for k = 1:size(p,1)
        
            L = 2*hL./(sqrt(sin(th0r(i)).^2+2.35*10^-4)+sin(th0r(i))); % [m] effective path lenght
            x = 1.22*Deff.^2*(f./L);
    
    
            if x >= 7.0 % argument of the square root is negative 
    
                Afd(k,j) = 0; 
                %fprintf('argument of the square root is negative')
    
            else

                g = sqrt(3.86*(x.^2+1)^(11/12)*sin(11/6*atan(1./x))-7.08*x.^(5/6)); % antenna averaging factor
                sigma = sigma_ref(k,j).*f^(7/12).*g/(sin(th0r(i))^1.2); % std deviation of the signal
                ap = -0.061*(log10(p(k)).^3)+0.072*(log10(p(k)).^2)-1.71*log10(p(k))+3.0; % time percentage factor 

                Afd(k,j) = ap.*sigma; % [dB]fade depth exceeded for p(%) of the time
               
            end
           
            if i == 1
            
                Afd__1(k,j) = Afd(k,j);
                
            end
        end
    end
    
    subplot(2,2,i)
    semilogy(Afd(:,1),p,Afd(:,2),p,Afd(:,3),p,'Linewidth',2)
    title(['P(Ascintillation)   θ = ',num2str(th0(i)),'°  f = ',num2str(f),'GHz'])
    lgd=legend('low (7.5°N)','mid (45°N)','high (64.5°N)');
    title(lgd,'Latitudes')
    xlabel('A[dB]')
    ylabel('p[%]')
    grid on
    
end

%% P(A)  θ = 5,20,45,90°  f = 28GHz

eta = 0.5;
D = 4.5; % [m] Viasat MEOlink
th0 = [5 20 45 90]; % [deg]
th0r = th0*pi/180; % [rad]
p = [0.1; 0.2; 0.3; 0.5; 1; 2; 3; 5; 10; 20; 30; 50]; % [%] time percentage
f = 28; % [GHz] downlink and uplink Viasat MEOLink
hL = 1000; % [m] height of turbulent layer
Deff = sqrt(eta)*D; % [m] effective antenna diameter
figure

for i = 1:size(th0,2)
    for j = 1:3
        for k = 1:size(p,1)
        
            L = 2*hL./(sqrt(sin(th0r(i)).^2+2.35*10^-4)+sin(th0r(i))); % [m] effective path lenght
            x = 1.22*Deff.^2*(f./L);
    
    
            if x >= 7.0 % argument of the square root is negative 
    
                Afd(k,j) = 0; 
                %fprintf('argument of the square root is negative')
    
            else

                g = sqrt(3.86*(x.^2+1)^(11/12)*sin(11/6*atan(1./x))-7.08*x.^(5/6)); % antenna averaging factor
                sigma = sigma_ref(k,j).*f^(7/12).*g/(sin(th0r(i))^1.2); % std deviation of the signal
                ap = -0.061*(log10(p(k)).^3)+0.072*(log10(p(k)).^2)-1.71*log10(p(k))+3.0; % time percentage factor 

                Afd(k,j) = ap.*sigma; % [dB]fade depth exceeded for p(%) of the time
               
            end
           
            if i == 1
            
                Afd__2(k,j) = Afd(k,j);
                
            end
        end
    end
    
    subplot(2,2,i)
    semilogy(Afd(:,1),p,Afd(:,2),p,Afd(:,3),p,'Linewidth',2)
    title(['P(Ascintillation)θ = ',num2str(th0(i)),'°  f = ',num2str(f),'GHz'])
    lgd = legend('low (7.5°N)','mid (45°N)','high (64.5°N)');
    title(lgd,'Latitudes')
    xlabel('A[dB]')
    ylabel('p[%]')
    grid on
    
end

%% A(θ) θ < 5° RX MID

%% 2.4.2 Amplitude scintillation fading for fades >= 25dB  θ < 5°
% Deep fading part th <= 5°

%% MID LAT p = 0.1%

%% 1. Apparent boresight elevation angle θ (tha) [mrad] (P.834)
% τ(h,th)     0<=h<=3    thm<=th<=5°  
% given formula for τ(h,th) is reasonable approx for 10°<th<=90°
% thm: angle at which radio beam is just intercepted by surface of Earth
% Input : th: elevation angle
%         h: rx altitude
%
% Output: tha: apparent boresight elevation angle [mrad]

th0 = [0.5 0.8 0.9 1 1.1 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2 2.5 3 3.5 4 4.5 5]; % [deg] free-space elevation angle
h = 0.5; % [km]
f = [19 28]; % [GHz]
a = 0.000315;
b = 0.1361;
r = 6370; % [km] Earth's radius
n = 1+a*exp(-b*h);
%p = [0.01 0.1 0.2 0.3 0.5 1 2 3 5 10 20 30 50]; % [%]
p = 0.1; % [%]

thm = -acos(r./(r+h)*(1+a)./n); % [rad]
taus(1,:) = 1./(1.314+0.6437*thm+0.02869*thm.^2+h.*(0.2305+0.09428*thm+0.01096*thm.^2)+h.^2.*0.008583);
th0v = thm-taus(1,:); % inequality verifier

for j = 1:size(th0,2)

    if  th0v <= th0(j)
    
        taus(2,:) = 1./(1.728+0.5411*th0(j)+0.03723*th0(j).^2+h.*(0.1815+0.06272*th0(j)+0.01380*th0(j).^2)+h.^2.*(0.01727+0.008288*th0(j))); % [deg]
        tha = th0(j)+taus(2,:); % [deg]
    
    else
    
        fprintf('space station not visible, apparent elevation angle not available')
    
    end

    thar(j) = tha*pi/180*1000; % [mrad]

%% 2. Geoclimatic factor Kw for avg worst month (P.453)
% 4 seasonally representative months from maps in Figs 8 to 11 (P.453)

Clat = ones(1,3);
C0 = ones(1,3);
lat = [7.5, 45, 64.5]; % [degN]
latr = lat*pi/180; % [rad]

%for i = 1:3
    if abs(lat(2)) <= 53 % [deg]
    
        Clat = 0;
        
    elseif abs(lat(2)) > 53 && abs(lat(2)) <= 60
        
        Clat = -53+lat(i);
        
    elseif abs(lat(2)) > 60
        
        Clat = 7;
        
    end
    
% Propagation path is assumed entirely over land   
    
    if h <= 0.7 % [km]
        
       C0 = 76;
       
    elseif h > 0.7
        
       C0 = 70; 
        
    end
%end

pL = [55, 8, 3]; % (%)of time that refractivity gradient in the lowest       
                 % 100m of atm is less than 100N units/km in that month
                 % (highest pL from 4 seasons is chosen)
                 % from maps (P.453) (depending on rx locations)

Kw = pL.^(1.5).*10.^((C0+Clat)/10);

%% 3. Fade depth A(p)
% valid for      0.5° < th < 5°
% p(%) input time percentage
% for avg Year
% RX MID

% v = ones(1,3);
% AP_y = ones(8,3);

    %for i = 1:3
        if abs(lat(2)) <= 45 % [deg]
                     
            v = -1.8-5.6*log10(1.1+abs(cos(2*latr(2))).^0.7); % [dB]

        elseif abs(lat(2)) > 45
        
            v = -1.8-5.6*log10(1.1-abs(cos(2*latr(2))).^0.7); % [dB]
    
        end
        
     for k = 1:size(f,2)
        
        AP_y(j,k) = 10*log10(Kw(2))-v+9*log10(f(k))-59.5*log10(1+thar(j))-10*log10(p); % [dB]
       
     end
    %end


% for avg annual worst-month

%AP_wm = ones(1,3);

    for k = 1:size(f,2)
    
        AP_wm(j,k) = 10*log10(Kw(2))+9*log10(f(k))-55*log10(1+thar(j))-10*log10(p); % [dB]
    
    end
end

for j = 1:size(AP_y,1)
    for  k= 1:size(AP_y,2)

        if AP_y(j,k) >= 25 % [dB]
        
            AP_y_(j,k) = AP_y(j,k);
            th0y_(j) = th0(j);
            
        end
        
         if AP_wm(j,k) >= 25 % [dB]
        
            AP_wm_(j,k) = AP_wm(j,k);
            th0wm_(j) = th0(j);
         
         end
    end
end

figure

    subplot(1,2,1)
    plot(th0y_(1:7),AP_y_(1:7,1),th0y_,AP_y_(:,2),'Linewidth',2)
    hold on
    title(['Ay(θ)  θ < 5°   A >= 25dB    p = ',num2str(p),'%'])
    lgd=legend(['f = ',num2str(f(1)),'GHz'],['f = ',num2str(f(2)),'GHz']);
    title(lgd,'Rx mid lat (45°N)')
    xlabel('θ[deg]')
    ylabel('A[dB]')
    grid on
    
    subplot(1,2,2)
    plot(th0wm_(1:16),AP_wm_(1:16,1),th0wm_,AP_wm_(:,2),'Linewidth',2)
    hold on
    title(['Awm(θ)  θ < 5°   A >= 25dB   p = ',num2str(p),'%'])
    lgd=legend(['f = ',num2str(f(1)),'GHz'],['f = ',num2str(f(2)),'GHz']);
    title(lgd,'Rx mid lat(45°N)')
    xlabel('θ[deg]')
    ylabel('A[dB]')
    grid on
    
   
%% all Lats  p = 0.1%  f = 19GHz

%% 1. Apparent boresight elevation angle θ (tha) [mrad] (P.834)
% τ(h,th)     0<=h<=3    thm<=th<=5°  
% given formula for τ(h,th) is reasonable approx for 10°<th<=90°
% thm: angle at which radio beam is just intercepted by surface of Earth
% Input : th: elevation angle
%         h: rx altitude
%
% Output: tha: apparent boresight elevation angle [mrad]

th0 = [0.5 0.8 0.9 1 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2 2.5 2.6 2.65 2.7 2.8 2.9 3 3.5 4 4.5 5]; % [deg] free-space elevation angle
th99 = th0;
h = 0.5; % [km]
f = 19; % [GHz]
a = 0.000315;
b = 0.1361;
r = 6370; % [km] Earth's radius
n = 1+a*exp(-b*h);
%p = [0.01 0.1 0.2 0.3 0.5 1 2 3 5 10 20 30 50]; % [%]
p = 0.1; % [%]

thm = -acos(r./(r+h)*(1+a)./n); % [rad]
taus(1,:) = 1./(1.314+0.6437*thm+0.02869*thm.^2+h.*(0.2305+0.09428*thm+0.01096*thm.^2)+h.^2.*0.008583);
th0v = thm-taus(1,:); % inequality verifier

for j = 1:size(th0,2)

    if  th0v <= th0(j)
    
        taus(2,:) = 1./(1.728+0.5411*th0(j)+0.03723*th0(j).^2+h.*(0.1815+0.06272*th0(j)+0.01380*th0(j).^2)+h.^2.*(0.01727+0.008288*th0(j))); % [deg]
        tha = th0(j)+taus(2,:); % [deg]
    
    else
    
        fprintf('space station not visible, apparent elevation angle not available')
    
    end

    thar(j) = tha*pi/180*1000; % [mrad]

%% 2. Geoclimatic factor Kw for avg worst month (P.453)
% 4 seasonally representative months from maps in Figs 8 to 11 (P.453)

Clat = ones(1,3);
C0 = ones(1,3);
lat = [7.5, 45, 64.5]; % [degN]
latr = lat*pi/180; % [rad]

for i = 1:3
    if abs(lat(i)) <= 53 % [deg]
    
        Clat = 0;
        
    elseif abs(lat(i)) > 53 && abs(lat(i)) <= 60
        
        Clat = -53+lat(i);
        
    elseif abs(lat(i)) > 60
        
        Clat = 7;
        
    end
    
% Propagation path is assumed entirely over land   
    
    if h <= 0.7 % [km]
        
       C0 = 76;
       
    elseif h > 0.7
        
       C0 = 70; 
        
    end


pL = [55, 8, 3]; % (%)of time that refractivity gradient in the lowest       
                 % 100m of atm is less than 100N units/km in that month
                 % (highest pL from 4 seasons is chosen)
                 % from maps (P.453) (depending on rx locations)

Kw(1,i) = pL(i).^(1.5).*10.^((C0+Clat)/10);

%% 3. Fade depth A(p)
% valid for      0.5° < th < 5°
% p(%) input time percentage
% for avg Year


% v = ones(1,3);
% AP_y = ones(8,3);

    
        if abs(lat(i)) <= 45 % [deg]
                     
            v0(1,i) = -1.8-5.6*log10(1.1+abs(cos(2*latr(i))).^0.7); % [dB]

        elseif abs(lat(i)) > 45
        
            v0(1,i) = -1.8-5.6*log10(1.1-abs(cos(2*latr(i))).^0.7); % [dB]
    
        end
        
    AP_y0(j,i) = 10*log10(Kw(i))-v0(i)+9*log10(f)-59.5*log10(1+thar(j))-10*log10(p); % [dB]
       
end
end


%{
for j = 1:size(AP_y,1)
    for  k= 1:size(AP_y,2)

        if AP_y(j,k) >= 25 % [dB]
        
            AP_y_(j,k) = AP_y(j,k);
            th0y_(j) = th0(j);
            
        end
        
         if AP_wm(j,k) >= 25 % [dB]
        
            AP_wm_(j,k) = AP_wm(j,k);
            th0wm_(j) = th0(j);
         
         end
    end
end
%}    
    
%% MID LAT

%% A(θ) θ < 5° transition zone  p = 0.1% rx mid lat A < 25dB
% sequencially after 2.4.1

%% 2.4.3 Amplitude scintillation in the transition region between the above distributions
% Shallow fading part
% free space th < 5°
% fading < 25dB
% from 2.4.2: Kw(pL, lat, h), v(lat) 

%% 1. Apparent elevation angle th1

A1 = 25; % [dB]
p = 0.1; % [%] chosen arbitrarily
f = [19 28]; % [GHz]

th1r(1,:) = ((Kw(2).*f.^0.9)./(p.*10.^(A1/10))).^(1/5.5)-1; % [mrad] worst month
th1r(2,:) = ((Kw(2).*10.^(-v/10)*f.^0.9)/(p.*10.^(A1/10))).^(1/5.95)-1; % [mrad] average year

th1 = th1r*180/(1000*pi); % [deg]

%% 2. A1'

A1_(1,:) = -55./(1+th1r(1,:))*log10(exp(1)); % [dB/mrad] worth month
A1_(2,:) = -59.5./(1+th1r(2,:))*log10(exp(1)); % [dB/mrad] average year

%% 3. A2 (P.453)

A2 = Afd_(:,1)'; % [dB] fade depth exceeded for p(%) of the time

%% 4. A2'

th0 = 5; % [deg]
th0r = th0*pi/180; % [rad]
hL = 1000; %[m]

eta = 0.5;
D = 4.5; % [m] antenna diameter
Deff = sqrt(eta)*D; % [m] effective antenna diameter

L = 2*hL./(sqrt(sin(th0r).^2+2.35*10^-4)+sin(th0r)); % [m] effective path lenght
x = 1.22*Deff.^2*(f./L);

Z = 11/6.*tan(1./x);
G_ = (1770.*(x.^2+1)+2123.*x.^(1/6).*(x.^2+1)*(11/12).*(cos(Z)-x.*sin(Z)))./(12*x.^(1/6).*(x.^2+1).*(354*x.^(5/6)-193*(x.^2+1).^(11/12).*sin(Z))); % g'/g
x_th = 1.22*Deff^2*f/(2*hL)*((sin(th0r))/(sqrt(sin(th0r)^2+2.35*10^(-4))+1))*cos(th0r); % dx/dθ
A2_= A2.*(G_.*x_th-1.2/tan(th0r))/1000; % [dB/mrad]

%% 5. Apparent elevation angle th2

h = 0.5; % [km] Earth station
th0 = 5; % [deg] elevation angle

taus = 1/(1.728+0.5411*th0+0.03723*th0.^2+h.*(0.1815+0.06272*th0+0.01380*th0.^2)+h.^2.*(0.01727+0.008288*th0)); % [deg]
th2 = (th0+taus); % [deg]
th2r = th2*(pi/180)*1000; % [mrad]

%% 6. Scintillation fading A(p) at desired apparent el. angle th
% Interpolation between (θ1, A1, A1') and (θ2, A2, A2') 
% θ1 <= θ <= θ2    0 <= p <= 50%

th = [1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2 2.05 2.1 2.15 2.2 2.25 2.3 2.35 2.4 2.45 2.5 3 3.5 4 4.5 5]'; % [deg] desired apparent elevation angle
thr = th*(pi/180)*1000; % [mrad]

delta1 = th2r-th1r(1,:);
alpha_p1 = A1_(1,:)./A1;
beta_p1 = (log(A2./A1)-alpha_p1.*delta1)./(delta1.^2);
gamma_p1 = (A2_-A2.*(alpha_p1+2*beta_p1.*delta1))./(A2.*delta1.^2);
    
delta2 = th2r-th1r(2,:);
alpha_p2 = A1_(2,:)./A1;
beta_p2 = (log(A2./A1)-alpha_p2.*delta2)./(delta2.^2);
gamma_p2 = (A2_-A2.*(alpha_p2+2*beta_p2.*delta2))./(A2.*delta2.^2);
    
for i = 1:size(th,1)
    
    A_p_wm(i,:) = A1.*exp(alpha_p1.*(thr(i)-th1r(1,:))+beta_p1.*(thr(i)-th1r(1,:)).^2+gamma_p1.*(thr(i)-th1r(1,:)).^2.*(thr(i)-th2r)); % [dB]
    A_p_y(i,:) = A1.*exp(alpha_p2.*(thr(i)-th1r(2,:))+beta_p2.*(thr(i)-th1r(2,:)).^2+gamma_p2.*(thr(i)-th1r(2,:)).^2.*(thr(i)-th2r)); % [dB]

end

figure
subplot(1,2,1)
plot(th(11:size(th,1)),A_p_y(11:size(th,1),1),th(13:size(th,1)),A_p_y(13:size(th,1),2),'Linewidth',2)
title('Ay(θ) < 25dB   θ < 5°  p = 0.1%')
lgd = legend(['f = ',num2str(f(1)),'GHz'],['f = ',num2str(f(2)),'GHz']);
title(lgd,'Rx mid lat(45°N)')
xlabel('θ[deg]')
ylabel('A[dB]')
grid on

subplot(1,2,2)
plot(th(18:size(th,1)),A_p_wm(18:size(th,1),1),th(21:size(th,1)),A_p_wm(21:size(th,1),2),'Linewidth',2)
title('Awm(θ) < 25dB   θ < 5°  p = 0.1%')
lgd = legend(['f = ',num2str(f(1)),'GHz'],['f = ',num2str(f(2)),'GHz']);
title(lgd,'Rx mid lat(45°N)')
xlabel('θ[deg]')
ylabel('A[dB]')
grid on

% tot scintillation f = 19GHz
figure
plot(th9,Afd1(1,:),th0y_(1:7),AP_y_(1:7,1),th(11:size(th,1)),A_p_y(11:size(th,1),1),'Color',[0 0.4470 0.7410],'Linewidth',2)
title('A(θ)scintillation full   p = 0.1%')
grid on
xlabel('θ[deg]')
ylabel('A[dB]')
hold on
plot(th9,Afd1(2,:),th0y_,AP_y_(:,2),th(13:size(th,1)),A_p_y(13:size(th,1),2),'Color',[0.8500 0.3250 0.0980],'Linewidth',2)
lgd=legend('Ka band downlink frequency (19GHz)','','','Ka band uplink frequency (28GHz)','','');
title(lgd,'Rx mid latitude (45°N)')

%% all Lats p = 0.1%  f = 19GHz

%% 1. Apparent elevation angle th1

A1 = 25; % [dB]
p = 0.1; % [%] chosen arbitrarily
f = 19; % [GHz]

for j = 1:3
%th1r(1,:) = ((Kw(i).*f.^0.9)./(p.*10.^(A1/10))).^(1/5.5)-1; % [mrad] worst month
th1r(1,j) = ((Kw(j).*10.^(-v0(j)/10)*f.^0.9)/(p.*10.^(A1/10))).^(1/5.95)-1; % [mrad] average year

th1 = th1r*180/(1000*pi); % [deg]

%% 2. A1'

%A1_(1,:) = -55./(1+th1r(1,:))*log10(exp(1)); % [dB/mrad] worth month
A1_(1,j) = -59.5./(1+th1r(1,j))*log10(exp(1)); % [dB/mrad] average year

%% 3. A2 (P.453)

A2 = Afd0_(1,j); % [dB] fade depth exceeded for p(%) of the time

%% 4. A2'

th0 = 5; % [deg]
th0r = th0*pi/180; % [rad]
hL = 1000; %[m]

eta = 0.5;
D = 4.5; % [m] antenna diameter
Deff = sqrt(eta)*D; % [m] effective antenna diameter

L = 2*hL./(sqrt(sin(th0r).^2+2.35*10^-4)+sin(th0r)); % [m] effective path lenght
x = 1.22*Deff.^2*(f./L);

Z = 11/6.*tan(1./x);
G_ = (1770.*(x.^2+1)+2123.*x.^(1/6).*(x.^2+1)*(11/12).*(cos(Z)-x.*sin(Z)))./(12*x.^(1/6).*(x.^2+1).*(354*x.^(5/6)-193*(x.^2+1).^(11/12).*sin(Z))); % g'/g
x_th = 1.22*Deff^2*f/(2*hL)*((sin(th0r))/(sqrt(sin(th0r)^2+2.35*10^(-4))+1))*cos(th0r); % dx/dθ
A2_= A2.*(G_.*x_th-1.2/tan(th0r))/1000; % [dB/mrad]

%% 5. Apparent elevation angle th2

h = 0.5; % [km] Earth station
th0 = 5; % [deg] elevation angle

taus = 1/(1.728+0.5411*th0+0.03723*th0.^2+h.*(0.1815+0.06272*th0+0.01380*th0.^2)+h.^2.*(0.01727+0.008288*th0)); % [deg]
th2 = (th0+taus); % [deg]
th2r = th2*(pi/180)*1000; % [mrad]

%% 6. Scintillation fading A(p) at desired apparent el. angle th
% Interpolation between (θ1, A1, A1') and (θ2, A2, A2') 
% θ1 <= θ <= θ2    0 <= p <= 50%

th = [1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2 2.05 2.1 2.15 2.2 2.25 2.3 2.35 2.4 2.45 2.5 2.55 2.6 2.65 2.7 2.75 2.8 2.85 2.9 2.95 3 3.2 3.4 3.6 3.8 4 4.2 4.4 4.6 4.8 5]'; % [deg] desired apparent elevation angle
thr = th*(pi/180)*1000; % [mrad]
th999 = th;

delta1 = th2r-th1r(1,j);
alpha_p1 = A1_(1,j)./A1;
beta_p1 = (log(A2./A1)-alpha_p1.*delta1)./(delta1.^2);
gamma_p1 = (A2_-A2.*(alpha_p1+2*beta_p1.*delta1))./(A2.*delta1.^2);

%{    
delta2 = th2r-th1r(2,:);
alpha_p2 = A1_(2,:)./A1;
beta_p2 = (log(A2./A1)-alpha_p2.*delta2)./(delta2.^2);
gamma_p2 = (A2_-A2.*(alpha_p2+2*beta_p2.*delta2))./(A2.*delta2.^2);
%}    
    for i = 1:size(th,1)
    
        A_p_y0(i,j) = A1.*exp(alpha_p1.*(thr(i)-th1r(1,j))+beta_p1.*(thr(i)-th1r(1,j)).^2+gamma_p1.*(thr(i)-th1r(1,j)).^2.*(thr(i)-th2r)); % [dB]
%   A_p_y(i,:) = A1.*exp(alpha_p2.*(thr(i)-th1r(2,:))+beta_p2.*(thr(i)-th1r(2,:)).^2+gamma_p2.*(thr(i)-th1r(2,:)).^2.*(thr(i)-th2r)); % [dB]

    end
end

%% PLOT A(θ) SCINTILLATION FADING OVER ALL STATIONS  p = 0.1%   f = 19Ghz

figure
%plot(th9,Afd0(:,1),th9,Afd0(:,2),th9,Afd0(:,3),th99(1:26),AP_y0(1:26,1),th99(1:9),AP_y0(1:9,2),th99(1:6),AP_y0(1:6,3),th999(36:size(th999)),A_p_y0(36:size(th999),1),th999(11:size(th999)),A_p_y0(11:size(th999),2),th999(8:size(th999)),A_p_y0(8:size(th999),3),'Linewidth',2)
plot(th9,Afd0(:,1),th99(1:26),AP_y0(1:26,1),th999(36:size(th999)),A_p_y0(36:size(th999),1),'Color',[0 0.4470 0.7410],'Linewidth',2)

title('A(θ)scintillation full   p = 0.1%   f = 19GHz')

grid on
xlabel('θ[deg]')
ylabel('A[dB]')

hold on

plot(th9,Afd0(:,2),th99(1:9),AP_y0(1:9,2),th999(11:size(th999)),A_p_y0(11:size(th999),2),'Color',[0.8500 0.3250 0.0980],'Linewidth',2)
plot(th9,Afd0(:,3),th99(1:6),AP_y0(1:6,3),th999(8:size(th999)),A_p_y0(8:size(th999),3),'Color',[0.9290 0.6940 0.1250],'Linewidth',2)


lgd = legend('low (7.5°N)','','','mid (45°N)','','','high (64.5°N)','','');
title(lgd,'Latitudes of Earth stations')



%% P(A) CLOUDS

%% P.840 Clouds and fog
% 10 < f < 200GHz (Rayleigh approximation limit)
% gammac = [dB/km] specific attenuation within the cloud
% Kl = [(dB/km)/(g/m^3)] cloud liquid water specific att. coefficient
% M = [g/m^3] liquid water density in cloud or fog (g/m^3)
% T = [K] cloud liquid water temperature
%

%% 1.
f = [19 28]; % [GHz]

Mmf = 0.05; % [g/m^3] medium fog density 
Mtf = 0.5; % [g/m^3] thick fog density
T = 273.15; % [K]
th = 300/T; % [K^-1]
eps0 = 77.66+103.3*(th-1);
eps1 = 0.0671*eps0;
eps2 = 3.52;
fp = 20.20-146*(th-1)+316*(th-1)^2; % [GHz] principal relaxation frequency
fs = 39.8*fp; % [GHz] secondary rel. fr.

gmid = ones(size(f,2),1);
gthick = ones(size(f,2),1);
KL = ones(size(f,2),1);

for i = 1:size(f,2)
    
    eps_ = (eps0-eps1)/(1+(f(i)/fp)^2)+(eps1-eps2)/(1+(f(i)/fs)^2)+eps2;
    eps__ = f(i)*(eps0-eps1)/(fp*(1+(f(i)/fp)^2))+f(i)*(eps1-eps2)/(fs*(1+(f(i)/fs)^2));
    eta = (2+eps_)/eps__;
    Kl = 0.819*f(i)/(eps__*(1+eta^2)); % [(dB/km)/(g/m^3)]

    gammac1 = Kl*Mmf; % [dB/km] specific attenuation within cloud
    gammac2 = Kl*Mtf; % [dB/km]

    gmid(i,1) = gammac1;
    gthick(i,1) = gammac2;
    KL(i,1) = Kl; 
    
end


%% 2.
% 5° <= fi <= 90° elevation angle
% Lred = [kg/m^2] total columnar content of liquid water reduced to T
%        from (P.840) Lred Annual Maps

delimiterIn = ' ';
filename = 'Lred_01_v4.txt';
Lred1 = importdata(filename,delimiterIn);
filename = 'Lred_02_v4.txt';
Lred2 = importdata(filename,delimiterIn);
filename = 'Lred_03_v4.txt';
Lred3 = importdata(filename,delimiterIn);
filename = 'Lred_05_v4.txt';
Lred4 = importdata(filename,delimiterIn);
filename = 'Lred_1_v4.txt';
Lred5 = importdata(filename,delimiterIn);
filename = 'Lred_2_v4.txt';
Lred6 = importdata(filename,delimiterIn);
filename = 'Lred_3_v4.txt';
Lred7 = importdata(filename,delimiterIn);
filename = 'Lred_5_v4.txt';
Lred8 = importdata(filename,delimiterIn);
filename = 'Lred_10_v4.txt';
Lred9 = importdata(filename,delimiterIn);
filename = 'Lred_20_v4.txt';
Lred10 = importdata(filename,delimiterIn);
filename = 'Lred_30_v4.txt';
Lred11 = importdata(filename,delimiterIn);
filename = 'Lred_50_v4.txt';
Lred12 = importdata(filename,delimiterIn);
filename = 'Lred_60_v4.txt';
Lred13 = importdata(filename,delimiterIn);
filename = 'Lred_70_v4.txt';
Lred14 = importdata(filename,delimiterIn);
filename = 'Lred_80_v4.txt';
Lred15 = importdata(filename,delimiterIn);
filename = 'Lred_90_v4.txt';
Lred16 = importdata(filename,delimiterIn);
filename = 'Lred_95_v4.txt';
Lred17 = importdata(filename,delimiterIn);
filename = 'Lred_99_v4.txt';
Lred18 = importdata(filename,delimiterIn);

R=[74 41 24];
C=[73 265 302];

Lred = zeros(18,4);

for i = 1:3

    Lred(:,i) = [Lred1(R(i),C(i)), Lred2(R(i),C(i)), Lred3(R(i),C(i)), Lred4(R(i),C(i)), Lred5(R(i),C(i)), Lred6(R(i),C(i)), Lred7(R(i),C(i)), Lred8(R(i),C(i)), Lred9(R(i),C(i)), Lred10(R(i),C(i)), Lred11(R(i),C(i)), Lred12(R(i),C(i)), Lred13(R(i),C(i)), Lred14(R(i),C(i)), Lred15(R(i),C(i)), Lred16(R(i),C(i)), Lred17(R(i),C(i)), Lred18(R(i),C(i))];

end

fi = [5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90]; % [deg] elevation angle
fir = fi*pi/180; % [rad]

clear A1 A2 A3

for i = 1:size(KL,1)

    A1(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(1));
    A2(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(2));
    A3(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(3));
    A4(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(4));
    A5(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(5));
    A6(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(6));
    A7(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(7));
    A8(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(8));
    A9(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(9));
    A10(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(10));
    A11(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(11));
    A12(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(12));
    A13(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(13));
    A14(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(14));
    A15(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(15));
    A16(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(16));
    A17(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(17));
    A18(:,4*i-3:4*i) = Lred.*KL(i)./sin(fir(18));

end

p = [0.1 0.2 0.3 0.5 1 2 3 5 10 20 30 50 60 70 80 90 95 99]; % [%]
% f = 19GHz

figure
subplot(2,2,1)
semilogy(A1(:,1),p,A1(:,2),p,A1(:,3),p,'Linewidth',2)
title(['P(Acloud)   θ  = ',num2str(fi(1)),'°    f = 19GHz'])
grid on
xlabel('A[dB]')
ylabel('p[%]')
lgd = legend('low (7.5°N)','mid (45°N)','high (64.5°N)');
title(lgd,'Latitudes')

subplot(2,2,2)
semilogy(A4(:,1),p,A4(:,2),p,A4(:,3),p,'Linewidth',2)
title(['P(Acloud)   θ  = ',num2str(fi(4)),'°    f = 19GHz'])
grid on
xlabel('A[dB]')
ylabel('p[%]')
lgd = legend('low (7.5°N)','mid (45°N)','high (64.5°N)');
title(lgd,'Latitudes')

subplot(2,2,3)
semilogy(A9(:,1),p,A9(:,2),p,A9(:,3),p,'Linewidth',2)
title(['P(Acloud)   θ  = ',num2str(fi(9)),'°    f = 19GHz'])
grid on
xlabel('A[dB]')
ylabel('p[%]')
lgd = legend('low (7.5°N)','mid(45°N)','high(64.5°N)');
title(lgd,'Latitudes')

subplot(2,2,4)
semilogy(A18(:,1),p,A18(:,2),p,A18(:,3),p,'Linewidth',2)
title(['P(Acloud)   θ  = ',num2str(fi(18)),'°    f = 19GHz'])
grid on
xlabel('A[dB]')
ylabel('p[%]')
lgd = legend('low (7.5°N)','mid (45°N)','high (64.5°N)');
title(lgd,'Latitudes')


%
figure
semilogy(A1(:,2),p,A1(:,6),p,'Linewidth',2)
title('P(Acloud)  θ = 5°')
grid on
xlabel('A[dB]')
ylabel('p[%]')
lgd = legend(['f = ',num2str(f(1)),'GHz'],['f = ',num2str(f(2)),'GHz']);
title(lgd,'Rx mid lat (45°N)')


%% 2.5 Total attenuation due to multiple sources of simultaneously occuring atm attenuation
% f > 18GHz
% low th
% total att is combined effect of rain, gas, clouds, scintillation
%
% Inputs : AR(p) = [dB] attenuation due to rain for fixed probability 2.2.1.1 (10.)
%          AC(p) = [dB] att. due to clouds (p.840)
%          AG(p) = [dB] att. due to water vapour and oxygen (P.676)
%          AS(p) = [dB] att. due to troposhpheric scintillation 2.4.1 (9.)
%          0.001 < p < 50% probability of attenutation being exceeded
%
% Output : AT(p) = [dB] total attenuation for given probability
%
%   Combine AR, AS, AC, showing the effects at differend stations 
%   (low, mid, high latitudes) at fixed elevation angle and frequency (worst case): 
%
% th = 5 [deg]
% f = 28 [GHz]

AR_ = zeros(18,3);
AS_ = AR_;
AC_ = AR_;
p = [0.1 0.2 0.3 0.5 1 2 3 5 10 20 30 50 60 70 80 90 95 99]; % [%]

AR_(1:8,:) = Ap2(3:size(Ap2),:); % [dB] rain
AS_(1:12,:) = Afd__2; % [dB] scintillation
AC_ = A1(:,5:7); % [dB] clouds and fog
AC0_ = AC_;
for i = 1:4
    AC_(i,:) = AC_(5,:);
end
AT_RSC_ = sqrt((AR_+AC_).^2+AS_.^2); % [dB]

figure
subplot(2,1,1)
semilogy(AT_RSC_(:,1),p,AT_RSC_(:,2),p,AT_RSC_(:,3),p,'Linewidth',2)
title('P(Atot)   θ = 5°    Ka band uplink freq. (28GHz)')
grid on
xlabel('Atot[dB]')
ylabel('p[%]')
lgd = legend('low (7.5°N)','mid (45°N)','high (64.5°N)');
title(lgd,'Latitude')

% Combined effects at different carrier frequencies, at fixed elevation
% angle and latitude:
% th = 5 [deg]
% Rx Mid Lat

AR__ = zeros(18,2);
AS__ = AR__;
AC__ = AR__;

AR__(1:8,1) = Ap1(3:size(Ap2),2); % [dB] rain
AR__(1:8,2) = Ap2(3:size(Ap2),2); % [dB] rain

AS__(1:12,1) = Afd__1(:,2); % [dB] scintillation
AS__(1:12,2) = Afd__2(:,2); % [dB] scintillation

AC__(:,1) = A1(:,2); % [dB] clouds and fog
AC__(:,2) = A1(:,6); % [dB] clouds and fog

for i = 1:4
    AC__(i,1) = A1(5,2); % [dB] clouds and fog
    AC__(i,2) = A1(5,6); % [dB] clouds and fog
end
AT_RSC__ = sqrt((AR__+AC__).^2+AS__.^2); % [dB]

subplot(2,1,2)
semilogy(AT_RSC__(:,1),p,AT_RSC__(:,2),p,'Linewidth',2)
title('P(Atot)   θ = 5°')
grid on
xlabel('Atot[dB]')
ylabel('p[%]')
lgd = legend('Ka band downlink frequency (19GHz)','Ka band uplink frequency (28GHz)');
title(lgd,'Rx Mid Lat (45°N)')


%plot of all attenuations, in rx mid, same conditions
figure

semilogy(AR_(:,2),p,AS_(:,2),p,AC0_(:,2),p,'Linewidth',2)
title('Attenuation statistics   θ = 5°    Ka band uplink freq. (28GHz)')
grid on
xlabel('A[dB]')
ylabel('p[%]')
lgd = legend('Rain','Scintillations','Clouds');
title(lgd,'Rx mid lat (45°N)')

% Table for p = 0.5%  f = 19GHz
AR(1:8,:) = Ap1(3:size(Ap1),:); % [dB] rain

AR = AR(4,:); % [dB]
AS = Afd__1(4,:); % [dB]
AC = A1(5,1:3); % [dB]
A05(1,:) = sqrt((AR+AC).^2+AS.^2);
A05(2,:) = AT_RSC_(4,:);
A05 % [dB]

A05t(1,:) = AR+AC;
A05t(2,:) = AR_(4,:)+AC_(4,:);

%% 3. Noise temperature
% Inputs : A = [dB] total atm attenuation excluding scintillation fading
%          Tmr = [K] atm mean radiating temperature
%
% Output: Tsky = [K] sky noise temperature at ground station antenna

p = [0.1 0.2 0.3 0.5 1 2 3 5 10 20 30 50 60 70 80 90 95 99]; % [%]

Tmr = 275; % [K] for clear and rainy weather
clear A
A_ = AR_+AC_; % [dB]
Tsky_ = Tmr*(1-10.^(-A_/10))+2.7*10.^(-A_/10); % [K]
A__ = AR__+AC__; % [dB]
Tsky__ = Tmr*(1-10.^(-A__/10))+2.7*10.^(-A__/10); % [K]

% plot worst case θ = 5° and f = 28 GHz

figure
subplot(2,1,1)
semilogy(Tsky_(:,1),p,Tsky_(:,2),p,Tsky_(:,3),p,'Linewidth',2)
title('P(Tsky)    θ = 5°    f = 28GHz')
grid on
xlabel('T[K]')
ylabel('p[%]')
lgd = legend('low (7.5°N)','mid (45°N)','high (64.5°N)');
title(lgd,'Latitude')

subplot(2,1,2)
semilogy(Tsky__(:,1),p,Tsky__(:,2),p,'Linewidth',2)
title('P(Tsky)    θ = 5° ')
grid on
xlabel('T[K]')
ylabel('p[%]')
lgd = legend('Ka band downlink frequency (19GHz)','Ka band uplink frequency (28GHz)');
title(lgd,'Rx Mid Latitude 45°N')

% Table
T05 =  Tmr*(1-10.^(-A05t./10))+2.7*10.^(-A05t./10) % [K]

toc