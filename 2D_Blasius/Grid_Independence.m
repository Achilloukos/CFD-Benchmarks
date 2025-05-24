function [CFD,NoGP] = Grid_Independence(dx,dy)

%% Data %%
len = 11; % [m]
height = 0.1; % [m]
visc = 1.57 * 10^-5; % [m^2/s], kinematic viscocity of air
rho = 1.2754; % [kg/m^3],  density of air at sea level 
U_inf = 1; % [m/s], air speed away from the plane
%% Discretization %%
% --> Uniform Grid
x = [0:dx:11];
y = [0:dy:height];
Noy = length(y);
Nox = length(x);

% Matrix initialization %
u = zeros(Noy,Nox);
u(:,1) = U_inf;
u(1,:) = U_inf;
u(end,:) = 0; % No-slip condition

v = zeros(Noy,Nox);

%% Solver (Explicit) %%
% --> Explicit Method


for i = 1:size(u,2)-1
    for j = size(u,1)-1:-1:2
        % Calculation of u at (i+1,j) via Momentum Equation
        c1 = visc * ((u(j-1,i) - 2*u(j,i) + u(j+1,i))/u(j,i)) * (dx/(dy)^2);
        c2 = 0.5*(u(j-1,i) - u(j+1,i)) * (v(j,i)/u(j,i)) * (dx/dy);
        u(j,i+1) = u(j,i) + c1 -c2;
        % Calculation of v at (i+1,j) via Continuity Equation
        c3 = (u(j,i+1) - u(j,i))/dx;
        c4 = (u(j+1,i+1) - u(j+1,i))/dx;
        v(j,i+1) = v(j+1,i+1) - 0.5*(c3 + c4)*dy;
    end
end
%% Boundary Layer metrics %%

% Boundary Layer thickness, displacement thickness, momentum thickness %

% --> Blasius analytical solution
delta99_Blasius = zeros(1,Nox);

for k = 1:Nox
    delta99_Blasius(k) = 5 * sqrt(visc * (k*dx)/U_inf);
end




% --> CFD Results
% --> explicit
delta_CFD_exp = zeros(1,Nox);

for i = 1:length(u)

    yind = find(u(:,i) >= 0.99*U_inf,1,'last');
    delta_CFD_exp(i) = dy*(Noy-(yind-1));
            
           
        
    
end


%% Visualization %%

% --> Explicit Method Results

% --> Boundary Layer thickness
CFD = delta_CFD_exp(round(Nox/2));
NoGP = Nox*Noy;

