%// ---------- V2X diffusion study

% Logic of first car is always different

clear all;
close all;
clc;

dt = 0.1;                                            % time step
N = 1500;                                            % number of time steps
n = 3;                                               % number of cars

% IDM parameters
% See https://github.com/movsim/traffic-simulation-de/blob/master/doc/TrafficSimulationDe.pdf
v0 = 80; % m/s                  desired speed of all vehicles
T = 2; % s                     reaction time
s0 = 15; % m                      minimum bumper-to-bumper distance
a = 0.2; % m/s^2                 acceleration in everyday traffic
b = 4.0; % m/s^2                 comfortable braking in everyday traffic

% Initial conditions
Vi = 7.5; % m/s

% vf = 80;
% rhom = 86;
% rhos = 43;
% sm = 1/rhom;

L= 70; % m                             Radius of Ring

d =  5; % m                             length of car
edges = [50:.5:80];                                            % define bin edges for histogram

% dx = [60, 23, 35];       % velocities
dx = [Vi, Vi, Vi];       % velocities
x = [0, 10, 55];    % positions
s(1) = x(2) - x(1) - d; % headways
s(2) = x(3) - x(2) - d;
s(3) = L - x(3) - d;

%x = [  1     2     3     4     5     6     7     8     9    10];

for k = 1:1:1                                %%%%%%%%%%%%%--->>>loop for random positions
    k
   % r = rand(1,10).*L;         % randomized initial positions
    
   % x = sort(r);
%    x = [0 .25 .5 .75 1 1.25 1.3 1.55 1.64 1.75];
    
%     R = randi([10,10],1)  % randomized position of the AV vehicle    
    
    for j = 1:1:N                                             % Time Loop
        
        
        for i  = 1:1:n                                              % Loop for number of cars
            % find the relative velocity
            if i == n  % if the car is lead car
%                 vl = dx(1) - dx(i);
                vl = dx(i) - dx(1);
            else
%                 vl = dx(i+1) - dx(i);
                vl = dx(i) - dx(i+1);
            end
            % Updates in this order:
            % desired gap
            % acceleration
            % velocity
            % position
            % headway
            % NOTE: Work out the indexing correctly!!!
            s_star = s0 + max(0, dx(i)*T + 0.5*dx(i)*vl/sqrt(a*b));
            dvdt(i) = a*(1- (dx(i)/v0))^4 - a*(s_star/s(i))^2;
            dx(i) = dx(i) + dvdt(i)*dt;
            if dx(i) < 0
                x(i) = x(i) - 0.5*dx(i)*dx(i)/dvdt(i);
                dx(i) = 0;
            else
                x(i) = x(i) + dx(i)*dt + 0.5*dvdt(i)*dt*dt;
            end
            % s(i) = s(i) - vl*dt; % this looks like it may be part of the issue
            % disp(i); disp(s(i))
            
            theta(i) = x(i)*2*pi/L; % x converted to angle
            if theta(i) > 2*pi
               theta(i) = theta(i) - 2*pi; 
            end
            
%             if (x(i) >= L) % Circular Loop Logic
%                 x(i) = x(i) - L;
%             end
        end
        s(1) = x(2) - x(1);
        s(2) = x(3) - x(2);
        s(3) = x(1) - x(3);
        if s(1) < 0
            s(1) = s(1) + L;
        end
        
        if s(2) < 0
            s(2) = s(2) + L;
        end
        
        if s(3) < 0
            s(3) = s(3) + L;
        end
        Theta(j,:) = theta(:);  
        X(j,:) = x(:);         % Location of each car at this time step
        V(j,:) = dx(:);        % Speed of each car at this time step
        S(j,:) = s(:);         % Headways
        t(j) = dt*(j-1);       % actual running time
 
    end
end


% xlabel('Bins','FontSize',16)
% %set(gca,'XTick',[0,.25,.5])
% ylabel('Probability','FontSize',16)
% title('Normalized Histogram','FontSize',16)
% figure
figure
for i = 1:1:n
    st(:,i) = sin(Theta(:,i))*L/(2*pi);
    ct(:,i) = cos(Theta(:,i))*L/(2*pi);
    plot3(st(:,i),ct(:,i),t)                              % Plotting Angle
    hold on
end
%plot(t, x_1, t, x_2,t, x_3,t, x_4,t, x_5);
xlabel('x','FontSize',16)
%set(gca,'XTick',[0,.25,.5])
ylabel('y','FontSize',16)
zlabel('Time (h)','FontSize',16)
title('Location of Vehicles','FontSize',16)
legend('First Car','Second Car','Third Car','Location','southeast')


% for i = 1:1:n                                           % Plotting Locations
%     plot(t, X(:,i))
%     hold on
% end
% %plot(t, x_1, t, x_2,t, x_3,t, x_4,t, x_5);
% xlabel('Time (h)','FontSize',16)
% %set(gca,'XTick',[0,.25,.5])
% ylabel('Location','FontSize',16)
% title('Location of Vehicles','FontSize',16)
% %set(gca,'YTick',[0, 5, 10, 15, 20, 25])
% %set(gca,'YLim',[5 55])
% %legend('\rho_1','\rho_2','\rho_3','\rho_4','\rho_5','Location','southeast')
% 
figure

for i = 1:1:n                                           % Plotting Speed
    plot(t, V(:,i))
    hold on
end
%plot(t, x_1, t, x_2,t, x_3,t, x_4,t, x_5);
xlabel('Time (s)','FontSize',16)
%set(gca,'XTick',[0,.25,.5])
ylabel('Speed (m/s)','FontSize',16)
title('Speed of Vehicles','FontSize',16)
legend('First Car','Second Car','Third Car','Location','southeast')
%set(gca,'YTick',[0, 5, 10, 15, 20, 25])
%set(gca,'YLim',[52 67])
% legend('Tenth Car','Nineth Car','Eight Car','Seventh Car','Sixth Car','Fifth Car','Fourth Car','Third Car','Second Car','First Car','Location','southeast')


% figure
% 
% plot(t, Entropy(:))                                     % Plotting Entropy

%plot(t, x_1, t, x_2,t, x_3,t, x_4,t, x_5);
% xlabel('Time (h)','FontSize',16)
%set(gca,'XTick',[0,.25,.5])
% ylabel('Entropy','FontSize',16)
% title('Entropy Analysis - Last Car is AV','FontSize',16)
%set(gca,'YTick',[0, 5, 10, 15, 20, 25])
%set(gca,'YLim',[-.1 16])
%legend('First Car','Second Car','Third Car','Fourth Car','Fifth Car','Sixth Car','Seventh Car','Eight Car','Nineth Car','Tenth Car','Location','southeast')

% Plot headways
figure
for i = 1:1:n
    plot(t, S(:, i))
    hold on
end
title('Headway of vehicles','FontSize',16)
xlabel('Time (s)', 'FontSize', 16)
ylabel('Headway (m)', 'FontSize', 16)
legend('First Car','Second Car','Third Car','Location','southeast')

% % First car
% figure
% 
% st1 = sin(Theta(:,1))*L/(2*pi);
% ct1 = cos(Theta(:,1))*L/(2*pi);
% plot3(st1,ct1,t)
% legend('First Car', 'Location','southeast')
% 
% % Second car
% figure
% 
% st2 = sin(Theta(:,2))*L/(2*pi);
% ct2 = cos(Theta(:,2))*L/(2*pi);
% plot3(st2,ct2,t)
% legend('Second Car','Location','southeast')
% 
% % Third car
% figure
% 
% st3 = sin(Theta(:,3))*L/(2*pi);
% ct3 = cos(Theta(:,3))*L/(2*pi);
% plot3(st3,ct3,t)
% legend('Third Car','Location','southeast')