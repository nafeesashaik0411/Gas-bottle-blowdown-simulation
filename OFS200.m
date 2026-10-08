clear
clc
close all

% ------- Importing Test data -----------------------

O = 'OFS Test data.xlsx';

% OFS test data
T2 = readtable(O, 'Sheet', '200 bar');
t260 = T2.TIME;
PGNR260 = T2.PGNR;
TGBN260 = T2.TGBN;
TSBN260 = T2.TSBN;
PRCMP260 = T2.PRCMP;
PRGU260 = T2.PRGU_2;
TRGU260 = T2.TRGU_2;
%--------------Calc. of mdot-------------------------
cd=0.8; A=pi*0.006^2/4; P=PGNR260*10^5; T=TGBN260+273; 
mdot=zeros(size(P)); m = zeros(size(P)); V = 117.5*10^-3; 

for i=1:length(PRGU260)
Z = refpropm('Z','T',T(i),'P',P(i)/1000,'Nitrogen'); 
M = refpropm('M','T',T(i),'P',P(i)/1000,'Nitrogen');
g = refpropm('K','T',T(i),'P',P(i)/1000,'Nitrogen'); 
R = 8314/M; 
m(i)=P(i)*V/(Z*R*T(i));
end
mdot(2:end) = -(diff(m))/0.1;

mdot_avg = mean(mdot(100:300));
%-------------------Model---------------------

P0 = PGNR260(1); % Intial Pressure 
Tg0 = TGBN260(1)+273; % Intial Gas Temperature
Tw0 = TSBN260(1)+273; % Intial Wall Temperature
V = 117*10^-3; % Volume
D = 608*10^-3; % Outer Diameter
d = 11.2*10^-3; % Thickness
M_bot = 65; % Mass of botttle
C_bot = 565; % Specific heat of bottle
k_wall = 6.8; % Thermal conductivity of Bottle
n = 100; % No. of time steps
tf = 41; % Final temperature 
GasName = 'Nitrogen'; % Gas Name
mdot = mdot_avg; % Mass flow rate

[t,x] = GasBottle(P0,Tg0,Tw0,V,D,d,M_bot,C_bot,k_wall,n,tf,GasName,mdot);
%-------------------Plotting---------------------

f1 = figure();
f1.Position = [0 450 650 310];
plot(t,x(:,1)/10^5,LineWidth=1.5); 
hold on;
plot(t260,PGNR260,LineWidth=1.5,LineStyle="--");
ylabel('Pressure (bar)')
xlabel('Time(s)')
legend('Model','Test')
title('Bottle pressure variation')


f2 = figure();
f2.Position = [0 50 650 310];
plot(t,x(:,2),LineWidth=1.5); 
hold on;
plot(t260,TGBN260+273,LineWidth=1.5,LineStyle="--");
ylabel('Temperature (K)')
xlabel('Time(s)')
legend('Model','Test')
title('Gas Temperature variation')


f3 = figure();
f3.Position = [750 50 650 310];
plot(t,x(:,3),LineWidth=1.5); 
hold on;
plot(t,x(:,4),LineWidth=1.5); 
plot(t260,TSBN260+273,LineWidth=1.5,LineStyle="--");
ylabel('Temperature (K)')
xlabel('Time(s)')
legend('Model-inner','Model-outer','Test')
title('Wall Temperature variation')

f4 = figure();
f4.Position = [750 450 650 310];
plot(t,x(:,5),LineWidth=1.5); 
hold on;
plot(t260,m,LineWidth=1.5,LineStyle="--");
ylabel('Mass (kg)')
xlabel('Time(s)')
legend('Model','Test')
title('Mass variation')








