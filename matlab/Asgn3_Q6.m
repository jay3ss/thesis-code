function Asgn3_Q6
%Asgn3_Q6
%    Sample code for optimal control problem Asgn3_Q6.
%    This problem is from D.E.Kirk's Optimal control theory: an introduction
%    example 5.1-1 on page 198 - 202

dt = 0.001;
sim_time = 2;

time = 0:dt:sim_time;

its = round(sim_time/dt);

X1 = [];
X2 = [];
U = [];

for t=0:its
x1 = -6.103+7.289*t+6.696*exp(-t)-0.593*exp(t);
x2 = 7.289-6.696*exp(-t)-0.593*exp(t);
u = 7.289*(1-exp(t))-6.103*exp(t);

X1 = [X1, x1];
X1 = [X1, x1];
U = [U, u];
end

% State equations
syms x1 x2 p1 p2 u;
Dx1 = x2;
Dx2 = -x2 + u;

% Cost function inside the integral
syms g;
g = 0.5*u^2;

% Hamiltonian
syms p1 p2 H;
H = g + p1*Dx1 + p2*Dx2;

% Costate equations
Dp1 = -diff(H,x1);
Dp2 = -diff(H,x2);

% solve for control u
du = diff(H,u);
sol_u = solve(du, u);

% Substitute u to state equations
Dx2 = subs(Dx2, u, sol_u);

% convert symbolic objects to strings for using 'dsolve'
eq1 = strcat('Dx1=',char(Dx1));
eq2 = strcat('Dx2=',char(Dx2));
eq3 = strcat('Dp1=',char(Dp1));
eq4 = strcat('Dp2=',char(Dp2));

sol_h = dsolve(eq1,eq2,eq3,eq4); %dsolve gives Symbolic solution of ordinary differential equations.

%% use boundary conditions to determine the coefficients
%    case a: (a) x1(0)=x2(0)=0; x1(2) = 5; x2(2) = 2;
conA1 = 'x1(0) = 0';
conA2 = 'x2(0) = 0';
conA3 = 'x1(2) = 5';
conA4 = 'x2(2) = 2';
sol_a = dsolve(eq1,eq2,eq3,eq4,conA1,conA2,conA3,conA4);

% Compare the Analytical solution from Kirk book and Matlab code
sol_book = {@(t)(7.289*t-6.103+6.696*exp(-t)-0.593*exp(t)),...
            @(t)(7.289-6.696*exp(-t)-0.593*exp(t))};
time = linspace(0,2,20);
s_book = [sol_book{1}(time); sol_book{2}(time)];

% plot both solutions
figure(1);
ezplot(sol_a.x1,[0 2]); hold on;
ezplot(sol_a.x2,[0 2]);
ezplot(-sol_a.p2,[0 2]);    % plot the control: u=-p2
plot(time, s_book,'*');
axis([0 2 -1.6 7]);
text(0.6,0.5,'x_1(t)');
text(0.4,2.5,'x_2(t)');
text(1.6,0.5,'u(t)');
xlabel('Time');
ylabel('State Trajectories');
title('Solution Comparison');
hold off;
print -depsc -r300 eg1a.eps
