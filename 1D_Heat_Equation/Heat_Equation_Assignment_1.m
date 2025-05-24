%CFD Pre-assignment
%This is a script in order to solve the 1-D heat equation for a linear beam

%% Data Import
%du/dt = a*d2u/d2x (1D heat equation)
a = 0.0001; %Thermal Diffusivity Constant
L = 1; % [m]

%% Spacetime Discretization
dx = 0.01;
dt = (1/2)*1;
Sol_Time = 2000;
x_disc = L/dx+1;
t_disc = Sol_Time/dt+1;

Temp_Field = zeros(Sol_Time/dt,L/dx);

%% Initial and Boundary Conditions
Temp_Field(1:end,1) = 20;
Temp_Field(1:end,end) = 100;
Temp_Field(1,2:end-1) = 20;
%% Solver
% for i = 2:Sol_Time/dt
%     for j = 2:L/dx-1
%         Temp_Field(i,j) = Temp_Field(i-1,j) + a*(dt/(dx)^2)*(Temp_Field(i-1,j+1) - 2*Temp_Field(i-1,j) + Temp_Field(i-1,j-1));
%         Perc_diff = abs(Temp_Field(i,:) - Temp_Field(i-1,:))./Temp_Field(i-1,:);
%         if max(Perc_diff) <= 0.1
%             break
%             break
%         end
%     end
% end
% 
crit = 1;
i = 2;

while i < 50000
    for j = 2:L/dx-1
        Temp_Field(i,j) = Temp_Field(i-1,j) + a*(dt/(dx)^2)*(Temp_Field(i-1,j+1) - 2*Temp_Field(i-1,j) + Temp_Field(i-1,j-1));
    end
    crit = max(abs(Temp_Field(i,:) - Temp_Field(i-1,:))./Temp_Field(i-1,:));
    if crit <= 0.001
        break
    end
    i = i+1;
end



%% Visualization
% figure(1)
% colormap hot
% imagesc(Temp_Field)
% colorbar
% xlabel('# x-discrete point [m]')
% ylabel('# Discrete time step')

figure(1)
cmap = colormap(jet(x_disc));
taf = 0;
for k = 1:i
    c = 0;
  
    for j = 1:L/dx-1
        color = cmap(round(Temp_Field(k,j)),:);
        plot([c c+j*dx],[3 3],'Color',color,'LineWidth',15)
        c = c + j*dx;
        hold on
        
    end
    taf = taf + dt;
    title(num2str(taf),'[s]')
    drawnow
    hold off
end

