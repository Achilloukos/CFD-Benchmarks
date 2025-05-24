%% CFD Assignment 1 %%
% --> Calculation of the boundary layer for a flat plane

%% Data %%
len = 11; % [m]
height = 0.1; % [m]
visc = 1.57 * 10^-5; % [m^2/s], kinematic viscocity of air
rho = 1.2754; % [kg/m^3],  density of air at sea level 
U_inf = 1; % [m/s], air speed away from the plane
%% Discretization %%
% --> Uniform Grid

dx = 0.001; % [m]
dy = 0.001; % [m]

Noy = height/dy + 1;
Nox = len/dx + 1;

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

%% Solver (Implicit) %%

% --> Implicit Method (Crank-Nicolson)

U = zeros(Noy,Nox);
U(:,1) = U_inf;
U(1,:) = U_inf;
U(end,:) = 0; % No-slip condition

V = zeros(Noy,Nox);




for i = 1:size(u,2)-1
    counter = 1;
    A = zeros(Noy-2,Noy-2);
    b = zeros(Noy-2,1);
    for j = size(U,1)-1:-1:2

        C1 = -0.5*visc/(dy^2);
        C2 = U(j,i)/dx + visc/(dy^2);
        C3 = -0.5*visc/(dy^2);
        C4 = 0.5*visc*((U(j-1,i) - 2*U(j,i) + U(j+1,i))/(dy^2)) - V(j,i)*(U(j-1,i) - U(j+1,i))*0.5/dy + ((U(j,i))^2)/dx;
        if counter > 1 && counter < Noy-2
            A(counter,counter-1) = C1;
            A(counter,counter) = C2;
            A(counter,counter+1) = C3;
            b(counter) = C4;
        elseif counter == 1
            A(1,1) = C2;
            A(1,2) = C3;
            b(1) = C4; % = C4 - 0, bcz u(end,i+1) = 0
        
        else 
            A(end,end-1) = C1;
            A(end,end) = C2;
            b(end) = C4 - C3*U_inf;
            counter = counter + 1;
        end
        counter = counter + 1;
    end
    
    % --> Solving the algebraic system for the u-velocities at x_i
    U(end-1:-1:2,i+1) = A\b;
    % Calculation of v at (i+1,j) via Continuity Equation (Explicitly)
    for j = size(U,1)-1:-1:2
        c3 = (U(j,i+1) - U(j,i))/dx;
        c4 = (U(j+1,i+1) - U(j+1,i))/dx;
        V(j,i+1) = V(j+1,i+1) - 0.5*(c3 + c4)*dy;
    end
    
end


%% Boundary Layer metrics %%

% Boundary Layer thickness, displacement thickness, momentum thickness %

% --> Blasius analytical solution
delta99_Blasius = zeros(1,Nox);

for k = 1:Nox
    delta99_Blasius(k) = 5 * sqrt(visc * (k*dx)/U_inf);
end

delta1_Blasius = 0.34*delta99_Blasius;
delta2_Blasius = 0.13*delta99_Blasius;


% --> CFD Results
% --> explicit
delta_CFD_exp = zeros(1,Nox);

for i = 1:length(u)
    for j = size(u,1)-1:-1:2
        if u(j,i)/U_inf >= 0.99
            delta_CFD_exp(i) = (Noy - j)*dy;
            
            break
        end
    end
end

delta1_CFD_exp = zeros(1,Nox);
for k = 1:Nox
    delta1_x = 0;
    for m = 1:Noy
        delta1_x = delta1_x + (1 - u(end+1-m,k)/U_inf)*dy;
    end
    delta1_CFD_exp(k) = delta1_x;
end

delta2_CFD_exp = zeros(1,Nox);
for k = 1:Nox
    delta2_x = 0;
    for m = 1:Noy
        delta2_x = delta2_x + (u(end+1-m,k)/U_inf)*(1 - u(end+1-m,k)/U_inf)*dy;
    end
    delta2_CFD_exp(k) = delta2_x;
end

% --> implicit
delta_CFD_imp = zeros(1,Nox);

for i = 1:length(U)
    for j = size(U,1)-1:-1:2
        if U(j,i)/U_inf >= 0.99
            delta_CFD_imp(i) = (Noy - j)*dy;
            
            break
        end
    end
end

delta1_CFD_imp = zeros(1,Nox);
for k = 1:Nox
    delta1_x = 0;
    for m = 1:Noy
        delta1_x = delta1_x + (1 - U(end+1-m,k)/U_inf)*dy;
    end
    delta1_CFD_imp(k) = delta1_x;
end

delta2_CFD_imp = zeros(1,Nox);
for k = 1:Nox
    delta2_x = 0;
    for m = 1:Noy
        delta2_x = delta2_x + (U(end+1-m,k)/U_inf)*(1 - U(end+1-m,k)/U_inf)*dy;
    end
    delta2_CFD_imp(k) = delta2_x;
end
     
        
% Skin Friction Coefficient %

% --> Blasius analytical solution
Cf_Blasius = zeros(1,Nox);
for k = 1:Nox
    Cf_Blasius(k) = 0.664 * sqrt(visc/((k*dx)*U_inf));
end
        
% --> CFD Results
Cf_CFD_exp = zeros(1,Nox);
Cf_CFD_imp = zeros(1,Nox);
for k = 1:Nox
    Cf_CFD_exp(k) = visc * (u(end-1,k)/dy)/(0.5*U_inf);
    Cf_CFD_imp(k) = visc * (U(end-1,k)/dy)/(0.5*U_inf);
end



%% Visualization %%

% --> Explicit Method Results

% --> Boundary Layer thickness
figure(1)
plot([0:dx:len],delta_CFD_exp)
hold on
plot([0:dx:len],delta99_Blasius)
xlabel('Distance from Leading Edge [m]')
ylabel('Boundary Layer thickness [m]')
title('Boundary Layer thickness for flow over a flat plane (Explicit Method)')
legend('CFD Results (Explicit)','Blasius Analytical Results')

% --> Displacement thickness 
figure(2)
plot([0:dx:len],delta1_CFD_exp)
hold on
plot([0:dx:len],delta1_Blasius)
xlabel('Distance from Leading Edge [m]')
ylabel('Boundary Layer thickness [m]')
title('Displacement thickness for flow over a flat plane (Explicit Method)')
legend('CFD Results (Explicit)','Blasius Analytical Results')

% --> Momentum thickness 
figure(3)
plot([0:dx:len],delta2_CFD_exp)
hold on
plot([0:dx:len],delta2_Blasius)
xlabel('Distance from Leading Edge [m]')
ylabel('Boundary Layer thickness [m]')
title('Momentum thickness for flow over a flat plane (Explicit Method)')
legend('CFD Results (Explicit)','Blasius Analytical Results')

% --> Skin Friction Coefficient
figure(4)
plot([0:dx:len],Cf_CFD_exp)
hold on
plot([0:dx:len],Cf_Blasius)
xlabel('Distance from Leading Edge [m]')
ylabel('Skin Friction Coefficient [-]')
title('Skin Friction Coefficient for flow over a flat plane (Explicit Method)')
legend('CFD Results (Explicit)','Blasius Analytical Results')

% --> Implicit Method Results

% --> Boundary Layer thickness
figure(5)
plot([0:dx:len],delta_CFD_imp)
hold on
plot([0:dx:len],delta99_Blasius)
xlabel('Distance from Leading Edge [m]')
ylabel('Boundary Layer thickness [m]')
title('Boundary Layer thickness for flow over a flat plane (Implicit Method)')
legend('CFD Results (Implicit)','Blasius Analytical Results')

% --> Displacement thickness 
figure(6)
plot([0:dx:len],delta1_CFD_imp)
hold on
plot([0:dx:len],delta1_Blasius)
xlabel('Distance from Leading Edge [m]')
ylabel('Boundary Layer thickness [m]')
title('Displacement thickness for flow over a flat plane (implicit Method)')
legend('CFD Results (implicit)','Blasius Analytical Results')

% --> Momentum thickness 
figure(7)
plot([0:dx:len],delta2_CFD_imp)
hold on
plot([0:dx:len],delta2_Blasius)
xlabel('Distance from Leading Edge [m]')
ylabel('Boundary Layer thickness [m]')
title('Momentum thickness for flow over a flat plane (implicit Method)')
legend('CFD Results (implicit)','Blasius Analytical Results')

% --> Skin Friction Coefficient
figure(8)
plot([0:dx:len],Cf_CFD_imp)
hold on
plot([0:dx:len],Cf_Blasius)
xlabel('Distance from Leading Edge [m]')
ylabel('Skin Friction Coefficient [-]')
title('Skin Friction Coefficient for flow over a flat plane (Implicit Method)')
legend('CFD Results (implicit)','Blasius Analytical Results')