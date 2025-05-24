
DxDy = [0.01:-0.0002:0.0008];
CFD_vec = zeros(1,length(DxDy));
NoGP_vec = CFD_vec;
for i = 1:length(DxDy)
    [CFD_vec(i),NoGP_vec(i)] = Grid_Independence(DxDy(i),DxDy(i));
end

delta99_Blasius = 5 * sqrt(1.57 * 10^-5 * (5.5)/1);

figure(1)
plot(NoGP_vec,CFD_vec,'-*')
hold on
plot(NoGP_vec,ones(size(CFD_vec))*delta99_Blasius)

ylabel('Boundary Layer Thickness [m]')
xlabel('Number of Grid Points [#]')
title('Boundary Layer thickness at x = 5.5 [m], for different grid sizes')
legend('CFD Results (Explicit)','Blasius Analytical Results')
