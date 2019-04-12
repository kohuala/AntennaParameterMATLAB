%%

% Guan, Huihua

%%

% Given a radiation pattern pd and a fixed phi, this program 
% will generate the theta_max, maxima, and nulls of the normalized 
% radiation pattern automatically as well as the directivity, 
% max directivity, HPBW, FNBW, HPBE, and FNBE.

%% Radiation Pattern Plot
phi = pi/2;
theta = -pi+0.01:0.01:pi-0.01;
pd = abs((sin(2*pi.*sin(theta)*sin(phi))./(sin(0.5*pi.*sin(theta)*sin(phi))))).^2;

figure(1)
polar(theta, pd)
figure(2)
plot(theta, pd)
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})

%% Determining Max of Radiation Pattern
[pd_max, theta_max_idx] = max(pd);
disp('The max of the radiation pattern (pd_max) is ')
disp(pd_max)
disp('The theta_max (degrees) of the radiation pattern is ')
disp(theta(theta_max_idx).*(180/pi))
disp('The phi_max (degrees) of the radiation pattern is ')
disp(phi.*(180/pi))

%% Plotting Normalized Radiation Pattern
n_pd = pd/pd_max;
figure(3)
plot(theta, n_pd)
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
figure(4)
polar(theta, n_pd)

%% Global and Local Maximas of Normalized Radiation Pattern
TF1 = islocalmax(n_pd);
maxs = n_pd(TF1);
figure(5)
plot(theta, n_pd, theta(TF1), maxs, 'r*')
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
disp('The maximas of the radiation pattern are: ')
disp(maxs)
disp('The associated theta (degrees) at the maximas are: ')
disp(theta(TF1).*(180/pi)) 
figure(6)
polar(theta, n_pd)
hold on;
polar(theta(TF1), maxs)

%% Nulls of the Normalized Radiation Pattern
threshold = 1e-04;
TF2 = islocalmin(n_pd);
mins = n_pd(TF2);
remove = mins > threshold | mins < -threshold;
mins(remove) = [];
TF2(remove) = [];
figure(7)
plot(theta, n_pd, theta(TF2), mins, 'r*')
grid on;
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
figure(8)
polar(theta, n_pd)
hold on;
polar(theta(TF2), mins, 'r*')
disp('The nulls of the radiation pattern are located at (degrees) ')
disp(theta(TF2).*(180/pi))

%% Maximum Directivity Calculations
U = pd;
Umax = max(U);

% Prad = radiation intensity integrated over the entire steradian
func = @(theta_, phi_) sin(theta_).*abs((sin(2*pi.*sin(theta_).*sin(phi_))./(sin(0.5*pi.*sin(theta_).*sin(phi_))))).^2;
Prad = integral2(func, 0, pi, 0, 2*pi);

Dmax = (4*pi*Umax)/(Prad);
D = (4*pi*U)/(Prad);

disp('The Prad is: ')
disp(Prad)
disp('The max directivity Dmax is: ')
disp(Dmax)

%% Half Power Beamwidth Calculations
% HPBW is at the theta angle where the normalized radiation intensity is 0.5 
% Since phi is fixed, this will be the elevation HPBW
phi = pi/2;
U_ = n_pd;
Umax = pd_max;

syms theta_1;
eqn = (abs((sin(2*pi.*sin(theta_1).*sin(phi))./(sin(0.5*pi.*sin(theta_1).*sin(phi))))).^2)./Umax == 0.5;
theta_h1 = vpasolve(eqn, theta_1, [0 pi]);
theta_h = 2*theta_h1*(180/pi);

figure(9);
plot(theta,U_,'-','color',[0, 0.4470, 0.7410])
grid on;
xlim([-pi pi])
ylim([0 1])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
hold on;
plot(theta, (0.5)/(theta_h1)*theta, '-r');
hold on;
plot(theta, (-0.5)/(theta_h1)*theta, '-r');

disp('The HPBW (degrees) is: ')
disp(theta_h)
%% First Null Beamwidth Calculations
syms theta_1;
eqn = abs((sin(2*pi.*sin(theta_1).*sin(phi))./(sin(0.5*pi.*sin(theta_1).*sin(phi))))).^2 == 0;
theta_n1 = vpasolve(eqn, theta_1, [0 pi]);
theta_n = 2*theta_n1*(180/pi);

figure(10);
plot(theta,U_,'-','color',[0, 0.4470, 0.7410])
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
hold on;
xline(double(theta_n1),'-r');
hold on;
xline(-double(theta_n1),'-r');

disp('The FNBW (degrees) is: ')
disp(theta_n)
%% Half Power Beam Efficiency Calculations

% Case 1: If theta_1 is chosen as the first null
P_theta1 = integral2(func, 0, double(theta_n1), 0, 2*pi);
Ptotal = Prad;
BE1 = P_theta1/Ptotal;

disp('The first null beam efficiency BE is: ')
disp(BE1)

% Case 2: If theta_1 is chosen as the half power
P_theta2 = integral2(func, 0, double(theta_h1), 0, 2*pi);
BE2 = P_theta2/Ptotal;
disp('The half power beam efficiency BE is: ')
disp(BE2)
