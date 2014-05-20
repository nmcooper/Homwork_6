%% Nathaniel Cooper Hwk 6
% Qustion 1

% Parameters
time=1:25; % 25 year simulation
a=0.001; % Attact rate
e=0.1; % prey to preditor coverstion efficieny
d=0.1; %predator birth rate
b=0.5; % prey birth rate
h=0.0015; % prey handiling time
K=15000; % carring capacity (1st value) 

n0=[8000;2000]; % Vi=8000 and Pi=2000 inital values

% ODE solver 
[T, Y] = ode45(@(t,y) LV_Pred_RM(y,b,a,e,d,K,h),time,n0); % SEE ATTACHED FUNCTION: LV_Pred_RM

figure;
subplot(1,2,1);
plot(T,Y); % plots P and V vs time
xlabel('Time'); ylabel('Abundance, Predator and Prey Populations');
legend({'V, Prey abundance','P, Predator abundance'}); 
subplot(1,2,2);
hold on
isoV = @(x) b/a + (b.*x.*(a*h -1/K - (a/K).*h.*x))./a  ; % Prey isocline
plot(1:8000,feval(isoV,1:8000))
isoP = d/(e*a - a*d*h); % predator isocline
plot(Y(:,1),Y(:,2),'r-','LineWidth',2)
xlabel('V, prey abundance')
ylabel('P, predator abundance')
line([isoP isoP],[0 10000],'Color','g');
legend({'Prey isocline','Population trajectory','Predator isocline'});  
ylim([0 10000]); xlim([0 8000])
hold off


%% Qusetion 2
% set matrix paramators for Jacobian

V=d/(a*e-a*d*h); % V at equalibrium 
P=(b/a)+(b*V(a*h-(1/K)-(a/K)*h*V)/a); % P at equlibrium

a11=b-(P*a)/(V*a*h+1)^2-(2*V*b)/K; % p(dv/dv)

a12=-(V*a)/(V*a*h+1); % P(dv/dp)

a21=(P*a*e)/(V*a*h+1)^2; %p(dp/dv)

a22=(V*a*e)/(V*a*h+1)-d; % P(dp/dv)

J=[a11  a12; a21  a22]; % Jacobian Matrix 

J_eig=eig(J); % egien valuse of Jacobian matrix












