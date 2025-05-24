%% CFD Assignment 1 Bonus %%
% --> Calculation of the boundary layer for a flat plane (non-uniform grid)

%% Data %%
len = 11; % [m]
height = 0.01; % [m]
visc = 1.57 * 10^-5; % [m^2/s], kinematic viscocity of air
rho = 1.2754; % [kg/m^3],  density of air at sea level 
U_inf = 1; % [m/s], air speed away from the plane

%% Discretization %%
% --> NON-Uniform Grid

dksi = 0.001; % [m]
deta = 0.001; % [m]

%Noy = height/deta + 1;
Noy = 101;
Nox = len/dksi + 1;

% Matrix initialization %
u = zeros(Noy,Nox);
u(:,1) = U_inf;
u(1,:) = U_inf;
u(end,:) = 0; % No-slip condition

v = zeros(Noy,Nox);

for k = 1:Noy
    y(k) = exp(deta*(k-1)) - 1;
end

eta = log(y+1);
eta = eta(end:-1:1);

%% Solver (Explicit) %%
% --> Explicit Method


for i = 1:size(u,2)-1
    for j = size(u,1)-1:-1:2
        % Calculation of u at (i+1,j) via Momentum Equation
        c1 = 0.5*(u(j-1,i) - u(j+1,i)) * (v(j,i)/u(j,i)) * (dksi/deta)* exp(-eta(j));
        c2 = exp(-2*eta(j)) * visc * dksi/u(j,i);
        c3 = ((u(j-1,i) - 2*u(j,i) + u(j+1,i)) * (1/(deta)^2) - (u(j-1,i) - u(j+1,i))/deta);
        u(j,i+1) = u(j,i) - c1 + c2*c3;
        % Calculation of v at (i+1,j) via Continuity Equation
        c4 = (u(j,i+1) - u(j,i))/dksi;
        c5 = (u(j+1,i+1) - u(j+1,i))/dksi;
        v(j,i+1) = v(j+1,i+1) - 0.5*(c4 + c5)*deta*exp(eta(j));
    end
end

%% Boundary Layer metrics %%

% Boundary Layer thickness %

% --> Blasius analytical solution
delta_Blasius = zeros(1,Nox);
for k = 1:Nox
    delta_Blasius(k) = 5 * sqrt(visc * (k*dksi)/U_inf);
end

% --> CFD Results
% --> explicit
delta_CFD_NU = zeros(1,Nox);

for i = 1:length(u)
    for j = size(u,1)-1:-1:2
        if u(j,i)/U_inf >= 0.99
            delta_CFD_NU(i) = y(end+1-j);
            
            break
        end
    end
end

%% Visualization %%

% --> Explicit Method Results

% --> Boundary Layer thickness
figure(1)
plot([0:dksi:len],delta_CFD_NU)
hold on
plot([0:dksi:len],delta_Blasius)
xlabel('Distance from Leading Edge [m]')
ylabel('Boundary Layer thickness [m]')
title('Boundary Layer thickness, non-uniform grid (Explicit Method)')
legend('CFD Results (Explicit)','Blasius Analytical Results')