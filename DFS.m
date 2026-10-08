clear
clc
close all

% ------- Importing Test data -----------------------

O = 'DFS Test data.xlsx';

% OFS test data
T2 = readtable(O, 'Sheet', 'Filtered2');
time = T2.Time;
Pg = T2.PGBH;
Tg = T2.TGBH;

%--------------Calc. of mdot-------------------------
cd=0.8; A=pi*0.006^2/4; P=Pg*10^5; T=Tg+273; 
mdot=zeros(size(P)); m = zeros(size(P)); V = 35.5*10^-3; 

for i=1:length(time)
Z = refpropm('Z','T',T(i),'P',P(i)/1000,'Hydrogen'); 
M = refpropm('M','T',T(i),'P',P(i)/1000,'Hydrogen');
g = refpropm('K','T',T(i),'P',P(i)/1000,'Hydrogen'); 
R = 8314/M; 
m(i)=P(i)*V/(Z*R*T(i));
end
mdot(2:end) = -(diff(m))/0.1;

mdot_avg = mean(mdot(1:15));
%-------------------Model---------------------

P0 = Pg(1); % Intial Pressure 
Tg0 = Tg(1)+273; % Intial Gas Temperature
Tw0 = 307; % Intial Wall Temperature
V = 35.5*10^-3; % Volume
D = 428*10^-3; % Outer Diameter
d = 14*10^-3; % Thickness
M_bot = 11; % Mass of botttle
C_bot = 565; % Specific heat of bottle
k_wall = 6.8; % Thermal conductivity of Bottle
n = 100; % No. of time steps
tf = 1.7; % Final temperature 
GasName = 'Hydrogen'; % Gas Name
mdot = mdot_avg; % Mass flow rate


[t,x] = GasBottle(P0,Tg0,Tw0,V,D,d,M_bot,C_bot,k_wall,n,tf,GasName,mdot);
%-------------------Plotting---------------------

f1 = figure();
f1.Position = [0 450 650 310];
plot(t,x(:,1)/10^5,LineWidth=1.5); 
hold on;
plot(time,Pg,LineWidth=1.5,LineStyle="--");
ylabel('Pressure (bar)')
xlabel('Time(s)')
legend('Model','Test')
title('Bottle pressure variation')


f2 = figure();
f2.Position = [0 50 650 310];
plot(t,x(:,2),LineWidth=1.5); 
hold on;
plot(time,Tg+273,LineWidth=1.5,LineStyle="--");
ylabel('Temperature (K)')
xlabel('Time(s)')
legend('Model','Test')
title('Gas Temperature variation')


% f3 = figure();
% f3.Position = [750 50 650 310];
% plot(t,x(:,3),LineWidth=1.5); 
% hold on;
% plot(t,x(:,4),LineWidth=1.5); 
% plot(t260,TSBN260+273,LineWidth=1.5,LineStyle="--");
% ylabel('Temperature (K)')
% xlabel('Time(s)')
% legend('Model-inner','Model-outer','Test')
% title('Wall Temperature variation')

f4 = figure();
f4.Position = [750 450 650 310];
plot(t,x(:,5),LineWidth=1.5); 
hold on;
plot(time,m,LineWidth=1.5,LineStyle="--");
ylabel('Mass (kg)')
xlabel('Time(s)')
legend('Model','Test')
title('Mass variation')








