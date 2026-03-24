function term_vel = terminalVelocity(max_area, max_circ_area, mass)
% This function calculates the terminal velocity of hydrometeors using a
% particle by particle method. Inputs include: max area of hydrometeor, the
% corresponding ellipse area, and its mass. 

% Constants
rho_air = 0.9; % density of air [kg/m^3]
g = 9.8; % gravitationa constant [m/s^2]
eta = 1.81*10^-5; % dynamic viscosity of air [kg/m*s] 
nu = eta/rho_air; % kinematic viscosity of air [m^2/s] 

% Derived values

area_ratio = max_circ_area./max_area;
x1 = 8*mass.*g*rho_air;
x2 = pi*eta^2;
x3 = area_ratio.^0.25;
X = x1.*x3./x2;
p1 = (1+0.1519.*X.^0.5).^0.5;
Re = 8.5*(p1-1).^2;
v = Re.*eta.*(pi./max_circ_area).^0.5;

term_vel = v./(2.*rho_air);