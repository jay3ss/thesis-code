%Author: Shaurya Agarwal
%EE5630 - Assignment2, Q4
%This is the main file which cals other functions

clc
clear all
close all

%Declaring the System Matrices after comparison with Riccatti Equation
A = [0 1; -2 1];
B = [0; 1];
Q = [2 3; 3 5];
R=1;
K0 = [20; 0; 0; 0];

%%%%%---------------------------
%This will calculate the differential equation of K(t) in symbolic form
%%%%%---------------------------

%Declaring K(t) matrix
%k12 and k21 are essentially equal as K(t) is symmetric
syms k11 k12 k21 k22
K=[k11 k12; k21 k22];
K_dot = -A.'*K - K*A + K*B*(R^-1)*B.'*K;


%%%%%---------------------------
%Using ODE45 to solve numerically
[T, K] = ode45(@(t,K)mRiccati(t, K, A, B, Q,R), [5 0], K0);

plot(T, K)
xlabel('Time') % x-axis label
ylabel('K(t)') % y-axis label
legend('K_{11}(t)','K_{12}(t)','K_{21}(t)','K_{22}(t)')
figure
%%%%%---------------------------


%%%%%---------------------------
%Now as we have K(t), lets simulate another ODE system to solve for the states
%Solving for states x1 and x2
X0=[1,1];
[T2, X] = ode45(@(t,X)MySys(t, K, X,R,B, T), [0 5], X0);
plot(T2, X)
xlabel('Time') % x-axis label
ylabel('States X(t)') % y-axis label
legend('X_{1}(t)','X_{2}(t)')
%

%%%%%---------------------------
% Once we have K(t) and X(t) we can get the u(t).
% Interesting thisng to note is that values of K(t) and X(t) are at different time intervals.
% Hence we use interpolation

for j = 1:1:501
    
    K1 = K(:,1);
    K2 = K(:,2);
    K3 = K(:,3);
    K4 = K(:,4);
    t_new=(j-1)/100;
    K1=interp1(T,K1,t_new);
    K2=interp1(T,K2,t_new);
    K3=interp1(T,K3,t_new);
    K4=interp1(T,K4,t_new);
    new_K=[K1 K2; K3 K4];
    
    X1 = X(:,1);
    X2 = X(:,2);
    X1=interp1(T2,X1,t_new);
    X2=interp1(T2,X2,t_new);
    new_X=[X1; X2];
    
    u(j) =-(R^-1)*B.'*new_K*new_X;
    T3(j)=t_new;
end

figure
plot(T3,u)
xlabel('Time') % x-axis label
ylabel('Optimal Control Action u*(t)') % y-axis label
legend('u*(t)')

function dKdt = mRiccati(t, K, A, B, Q, R)

K = reshape(K, size(A)); %Converting K from "n^2"-by-1 to "n"-by-"n"
dKdt = -Q -A.'*K - K*A + K*B*(R^-1)*B.'*K; %Differential Equation
dKdt = dKdt(:); %Convert from "n"-by-"n" back to "n^2"-by-1
end

function dXdt = MySys(t, K, X,R,B, T)


%%%Interpolation is required for K(t) as time steps are different for
%%%diferent ODE solvers

    K1 = K(:,1);
    K2 = K(:,2);
    K3 = K(:,3);
    K4 = K(:,4);

    K1=interp1(T,K1,t);
    K2=interp1(T,K2,t);
    K3=interp1(T,K3,t);
    K4=interp1(T,K4,t);

    new_K=[K1 K2; K3 K4];

    dXdt = zeros(2,1);
    dXdt(1) = X(2);
    dXdt(2) = -X(1)-2*X(2)+(-(R^-1)*B.'*new_K*X);
end