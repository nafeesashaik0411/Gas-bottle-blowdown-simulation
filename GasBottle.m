function [t,x] = GasBottle(P0,Tg0,Tw0,V,D,d,M_bot,C_bot,k_wall,n,tf,GasName,mdot)

% Pressure input in bar
% Temperature input in K
% Volume of bottle in m^3
% Dia of bottle in m
% mdot in Kg/s
% thickness (d) in m
% M_bot - Kg
% C_bot -SI

%------------------------------Inputs--------------------------
Tg = Tg0; % gas temp
%-----Spliting the bottle in to two nodes for conduction------
Twi = Tw0; % internal wall temp
Two = Tw0; % External wall temp
Mi = M_bot/2; % internal wall Mass
Mo = M_bot/2; % external wall Mass
%-------------------------------------------------------------
P = P0*10^5; % Pressure, KPa
A = pi*D^2; % Area of Spherical gas bottle (inside)
Ao = pi*(D+2*d)^2; % Area of Spherical gas bottle (outside)
Aor = pi*0.012^2/4; % Area of bottle outlet
u = refpropm('U','T',Tg,'P',P/1000,GasName); % Internal energy
he = refpropm('H','T',Tg,'P',P/1000,GasName); % Exit Enthalpy
Z = refpropm('Z','T',Tg,'P',P/1000,GasName); % Compressibility factor
M = refpropm('M','T',Tg,'P',P/1000,GasName); % Molar mass
R = 8314/M; % Gas Const
m = P*V/(Z*R*Tg); 

tstep = tf/(n-1);
P_x = zeros(1,n+1);
T_gx = zeros(1,n+1);
T_wxi = zeros(1,n+1);
T_wxo = zeros(1,n+1);
m_x = zeros(1,n+1);

P_x(1) = P;
T_gx(1) = Tg;
T_wxi(1) = Twi;
T_wxo(1) = Two;
m_x(1) = m;

for i = 1:1:n 
%------------------------------Gas properties------------------
g = 9.81; % Acc. due to gravity
Tf = (Twi+303)/2; % Internal film temperature: wall/gas
gamma = refpropm('K','T',Tg,'P',P/1000,GasName); % Specific Heat
Z = refpropm('Z','T',Tg,'P',P/1000,GasName); % Compressibility factor
nu = refpropm('$','T',Tf,'P',P/1000,GasName)*10^-4; % Kinematic viscosity
beta = refpropm('B','T',Tf,'P',P/1000,GasName); % Volumetric expansivity
k = refpropm('L','T',Tf,'P',P/1000,GasName); %Thermal conductivity
Pr = refpropm('^','T',Tf,'P',P/1000,GasName); % Prandl No.
Gr = g*beta*abs(Twi-Tg)*D^3/nu^2; % Grashoff No.
Ra = Gr*Pr; % Rayleigh No.
rho = refpropm('D','T',Tf,'P',P/1000,GasName); % Density.

%------------------------------Heat Transfer-----------------------------

% Conduction Heat transfer
ri = D/2; % Internal radius of bottle
ro = (D+2*d)/2; % External radius of bottle
Q_cond = 4*pi*k_wall*ri*ro/(ro-ri)*(Two-Twi); 

% Heat transfer Coeff: Assuming Gas botlle is spherical for gas -
% Convection Inside bottle
vb = sqrt(g*beta*abs(Twi-Tg)*D); % Buoyancy velocity
vj = mdot/(rho*Aor); % Jet velocity
alpha = 0.1; % eff. circulation factor varies from 0.05-0.15
Uc = alpha*vj; % effective circulation velocity
v = sqrt(Uc^2+vb^2); % Combination of Jet and buoyancy
Re = v*D/nu; % Reynolds No.
Ri = Gr/Re^2; % Richardson No.

if Ri<0.1 % Forced Convection
    Nu = 0.023*Re^0.8*Pr^0.4;
elseif Ri>=0.1 && Ri<=10 % Mixed Convection
    Nu_n = 2+0.589*Ra^(1/4)/(1+(0.469/Pr)^(9/16))^(4/9); % Natural 
    Nu_f = 0.023*Re^0.8*Pr^0.4; % Forced
    a = 3;
    Nu = (Nu_n^a+Nu_f^a)^(1/a); 
elseif Ri>10 % Natural Convection
    Nu = 2+0.589*Ra^(1/4)/(1+(0.469/Pr)^(9/16))^(4/9);
end

h = k/D*Nu;   
Q = h*A*(Twi-Tg);
%----------------------------Mass flow rate calc.(Orifice)----------------
% Cd = 0.8; Aor = pi*(6E-3)^2/4; 
% mdot = Cd*P*Aor*sqrt(gamma/(Z*R*Tg))*(2/(gamma+1))^((gamma+1)/(2*(gamma-1)));

%------------------------------Energy Eqn--------------------------------
m_new = m - mdot*tstep;
U_new = m*u + (Q - mdot*he)*tstep;
u = U_new/m_new;
m = m_new;
Tg = refpropm('T','P',P/1000,'U',u,GasName);

%------------------------------Wall Temp---------------------------------
% External Heat transfer
Tf0 = (Two+303)/2; % External film temperature: wall/ambient
rho0 = refpropm('D','T',Tf0,'P',1.103*10^2,'Nitrogen'); % Density of air
ko = refpropm('L','T',Tf0,'P',1.103*10^2,'Nitrogen'); % Thermal conductivity of air
Pro = refpropm('^','T',Tf0,'P',1.103*10^2,'Nitrogen'); % Prandl No. of air
nuo = refpropm('$','T',Tf0,'P',1.103*10^2,'Nitrogen')*10^-4; % Kinematic viscosity
betao = refpropm('B','T',Tf0,'P',1.103*10^2,'Nitrogen'); % Volumetric expansivity
mu_inf = refpropm('V','T',303,'P',1.103*10^2,'Nitrogen'); % Dynamic viscosity
Gro = g*betao*abs(Two-303)*D^3/nuo^2; % Grashoff No. of air
Rao = Gro*Pro; % Rayleigh No. of air

% Heat transfer Coeff: Assuming Gas botlle is spherical for air
% External convection 

v = 1; % External Air velocity
Reo = v*D/nuo; % Reynolds No.
Rio = Gro/Reo^2; % Richardson No.
mu_s = refpropm('V','T',Two,'P',1.103*10^2,'Nitrogen'); % Dynamic viscosity

if Rio<0.1 % Forced Convection
    Nuo = 2+(0.4*Reo^(1/2)+0.06*Reo^(2/3))*Pro^0.4*(mu_inf/mu_s)^(1/4);
elseif Rio>=0.1 && Rio<=10 % Mixed Convection
    Nu_n = 2+0.589*Rao^(1/4)/(1+(0.469/Pro)^(9/16))^(4/9); % Natural 
    Nu_f = 2+(0.4*Reo^(1/2)+0.06*Reo^(2/3))*Pro^0.4*(mu_inf/mu_s)^(1/4); % Forced
    a = 3;
    Nuo = (Nu_n^a+Nu_f^a)^(1/a); 
elseif Rio>10 % Natural Convection
    Nuo = 2+0.589*Rao^(1/4)/(1+(0.469/Pro)^(9/16))^(4/9);
end

ho = ko/(D+2*d)*Nuo;
Cp = C_bot+0.5*(Twi-293); % Cp dependency on temp.

dTwidt = (Q_cond-h*A*(Twi-Tg))/(Mi*Cp);
Twi = Twi+dTwidt*tstep;

dTwodt = (ho*Ao*(303-Two)-Q_cond)/(Mo*Cp);
Two = Two+dTwodt*tstep;

%--------------------------Mass Updation--------------------------------

P = Z*m*R*Tg/V;

Pguess = P;
for j = 1:20
    Z = refpropm('Z','T',Tg,'P',Pguess/1000,GasName);
    Pnew = Z*m*R*Tg/V;
    if abs(Pnew-Pguess) < 1      % 1 Pa tolerance
        break
    end
    Pguess = Pnew;
end
P = Pnew;


he = refpropm('H','T',Tg,'P',P/1000,GasName);
%--------------------------Assignments----------------------------------

P_x(i+1) = P;
T_gx(i+1) = Tg;
T_wxi(i+1) = Twi;
T_wxo(i+1) = Two;
m_x(i+1) = m;

end

x = [P_x' T_gx' T_wxi' T_wxo' m_x'];
t = linspace(0,tf,n+1);
end