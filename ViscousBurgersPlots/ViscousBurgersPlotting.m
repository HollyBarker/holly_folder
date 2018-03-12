hold on
%initial solution in black: stars are analytical, line is oomph; t=0
plot(x0_valueana,sin(x0_valueana),'k*')
plot(x_oomph,initial_oomph,'k')

%300th timestep in red; t=0.3
plot(x0_valueana, solution_ana(:,3),'r*')
plot(x_oomph,solution_oomph(:,3),'r')

%800th timestep in blue; t=0.8
plot(x0_valueana, solution_ana(:,8),'b*')
plot(x_oomph,solution_oomph(:,8),'b')

%1000th timestep in green; t=1.0
plot(x0_valueana, solution_ana(:,10),'g*')
plot(x_oomph,solution_oomph(:,10),'g')

xlim([0,2*pi])

title({' ','1D Viscous Burgers Equation: ',' oomph-lib and Cole-Hopf Transformation Solutions',' '})

xlabel('x')
ylabel('u')