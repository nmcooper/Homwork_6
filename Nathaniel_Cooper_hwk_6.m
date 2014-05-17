%% Nathaniel Cooper Hwk 6
% Qustion 1

% Parameters
time=1:25; % 25 year simulation
a=0.001; % Attact rate
e=0.1; % prey to preditor coverstion efficieny
d=0.1; %predator birth rate
b=0.5; % prey birth rate
h=0.0015; % prey handiling time
k=15000; % carring capacity (1st value) 

n0=[8000;2000]; % Vi=8000 and Pi=2000 inital values

% ODE solver 
[T, Y] = ode45(@(t,y) LV_Pred_RM(y,b,a,e,d,K,h),time,n0); % solves de

figure;
subplot(1,2,1);
plot(T,Y); % plots P and V vs time
xlabel('Time'); ylabel('Abundance, Predator and Prey Populations');


