%% Nathaniel Cooper Hwk 6
% Qustion 1

% Parameters
time=1:25; % 25 year simulation
a=0.0001; % Attack rate
e=0.1; % prey to preditor coverstion efficieny
d=0.1; %predator birth rate
b=0.5; % prey birth rate
h=0.0015; % prey handiling time
K=1000000; % carrying capacity (1st value) 

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
plot(1:15000,feval(isoV,1:15000))
isoP = d/(e*a - a*d*h); % predator isocline
plot(Y(:,1),Y(:,2),'r-','LineWidth',2)
xlabel('V, prey abundance')
ylabel('P, predator abundance')
line([isoP isoP],[0 10000],'Color','g');
legend({'Prey isocline','Population trajectory','Predator isocline'});  
ylim([0 10000]); xlim([0 15000])
hold off
disp(' as carrying capacity increases the population becomes less stable, stable limit cycles approch the axis ')
disp(' if K= 100,000 it is likly that the system will go extinct givin time');
disp (' if K=10,000,000 then it is more likly that the system will go extinct givin time');  

%% Qusetion 2
% set matrix paramators for Jacobian solve using simbolic math 
syms a e d b h K V P % 
aa = 0.0001;      % attack rate
ee = 0.1;       % conversion efficiency rate
dd = 0.1;       % predator death rate
bb = 0.5;       % prey birth rate
hh = .5;         % prey handling time
KK = 50000;       % prey carrying capacity

% L-V growth equations for Lotka Volterra model
dV_dt = b*V*(1-V/K) - (a*V*P)/(1+a*h*V);
dP_dt = e*((a*V*P)/(1+a*h*V)) -d*P;

V0 = simplify(solve(dP_dt,V)); % value of V when P at equilibtrium (P isocline)
P0 = simplify(solve(dV_dt,P)); % value of P when V at equilibrium (V isocline)

V0 = subs(V0,[a e d h],[aa ee dd hh]);
P0 = subs(P0,[a b h K V],[aa bb hh KK V0]);

a11a = b-(P*a)/(V*a*h+1)^2-(2*V*b)/K; % P(dv_dv)
a11 = subs(a11a,[a e d b h K V P],[aa ee dd bb hh KK V0 P0]);

a12a = -(V*a)/(V*a*h+1); % P(dv_dp)
a12 = subs(a12a,[a e d b h K V P],[aa ee dd bb hh KK V0 P0]);

a21a =(P*a*e)/(V*a*h+1)^2; % P(dp_dV)
a21 = subs(a21a,[a e d b h K V P],[aa ee dd bb hh KK V0 P0]);

a22a =(V*a*e)/(V*a*h+1)-d; % P(dP_dP)
a22 = subs(a22a,[a e d b h K V P],[aa ee dd bb hh KK V0 P0]);

J = [a11 a12; a21 a22] % Jacobian matrix 
J_eig = eig(J) % eigan values of the jacobian matrix 

disp(' for K=15,000 the eigen values are -0.7610, and 0.0110');
disp(' for K=50,000 the eigen valuea are -0.0250-0.1199i and -0.0250-0.1199i');
disp(' for K=100,000 the eigen values are 0.0500-0.1323i and 0.0500+0.1323i');
disp(' for K=1,000,000 the eigen values are 0.1175-0.1034i and 0.1175+0.1034i');
disp(' for K=10,000,000 the eigen values are 0.1242-0.0975i and 0.1242+0.0975i');
disp(' the eigen analyses shows a progressivly unstable system with more variable stable limit cycles like the numarerical simulations');











