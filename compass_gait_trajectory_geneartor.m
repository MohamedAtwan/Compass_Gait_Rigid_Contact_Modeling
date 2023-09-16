clear; clc;close all;

N = 100;
k = 10;
m = 1; %kg
mH = 5; %kg
g = 9.81; %m/s^2
T = 1;
gamma = 3*pi/180;
Q = eye(2,2);
% Q(1,1) = 1e-6;
% Q(2,2) = 1e-6;

fun = @(x)trajectoryCostFun(x);
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations',5e6, 'MaxIterations',5e6);
A = [];
b = [];
Aeq = [];
beq = [];
lb = ones(2,k+1)*-0.6;
ub = ones(2,k+1)*0.6;
nonlcon = @dribbelConstraints;
% xb0 = zeros(1,2*(k+1))*1e-3;
xb0 = [ -0.128745118459292        -0.171165036493644         0.123433540956095         0.193546813107829 0.390226938963612       -0.0200657208747879        0.0046241028487897         -0.19248104952802 -0.6         0.500625529479417                     -0.05,...
      -0.282077607059959         0.307745090395092         0.598803817739429        -0.471501041471335 -0.338914116930161        -0.197392460596769         0.578140065402312        -0.105168179934522 -0.290363144092493      0.000727576548405127                       0.1];
xb = fmincon(fun,xb0,A,b,Aeq,beq,lb,ub,nonlcon,options);
xb = reshape(xb,[k+1,2]);
xb = xb';
q = zeros(2,N);
qd = zeros(2,N);
qdd = zeros(2,N);

qr = zeros(2,N);
qrd = zeros(2,N);
qrdd = zeros(2,N);

tau = zeros(2,N);
t = linspace(0,T,N);
c2 = zeros(2,N);
for tm = 1:N
    ts = ones(k+1,1);
    for i = k:-1:1
        ts(k-i+1) = t(tm)^i;
    end

    q(:,tm) = xb*ts;
    bd = [modifyPolynomial(polyder(xb(1,:)),k); modifyPolynomial(polyder(xb(2,:)),k)];
    qd(:,tm) = bd*ts;
    bdd = [modifyPolynomial(polyder(bd(1,:)),k); modifyPolynomial(polyder(bd(2,:)),k)];
    qdd(:,tm) = bdd*ts;
    tau(:,tm) = inverseDribbelDynamics(q(:,tm),qd(:,tm),qdd(:,tm));

    c2_temp = forwardKinematicsLeft([0;0;q(:,tm)]);
    c2(:,tm) = c2_temp(2:end);
end

tau = tau';
e = abs(2*q(1,:)+q(2,:));
q_t = timeseries(q',linspace(0,T,N));
qd_t = timeseries(qd',linspace(0,T,N));
qdd_t = timeseries(qdd',linspace(0,T,N));
e_t = timeseries(e',linspace(0,T,N));
c2_t = timeseries([zeros(N,1) c2'],linspace(0,T,N));

NumCycles = 2;
tau_t = timeseries(repmat(tau,NumCycles,1),linspace(0,T*NumCycles,N*NumCycles));
% q_t = timeseries(repmat(q',NumCycles,1),linspace(0,T*NumCycles,N*NumCycles));
% qd_t = timeseries(repmat(qd',NumCycles,1),linspace(0,T*NumCycles,N*NumCycles));
% qdd_t = timeseries(repmat(qdd',NumCycles,1),linspace(0,T*NumCycles,N*NumCycles));

disp(['q1: ', num2str(q(:,1)')]);
disp(['q2: ', num2str(q(:,end)')]);

disp(['qd1: ', num2str(qd(:,1)')]);
disp(['qd2: ', num2str(qd(:,end)')]);
MPM = [(-m-2*m*cos(2*q(1,end))+4*(m+mH)*cos(4*q(1,end)))/(3*m+4*mH-2*m*cos(4*q(1,end))),                 (-2*m*cos(2*q(1,end)))/(3*m+4*mH-2*m*cos(4*q(1,end)));
       (8*(m+mH)*(1+2*cos(2*q(1,end)))*sin(q(1,end))^2)/(3*m+4*mH-2*m*cos(4*q(1,end))),                  (-m+2*m*cos(2*q(1,end)))/(3*m+4*mH-2*m*cos(4*q(1,end)))];

disp(['qd2: ', num2str((MPM*qd(:,end))')]);




figure()
plot(t,q(1,:))
hold on
plot(t,q(2,:))
hold off
legend('q3','q4')
xlabel('Time')
ylabel('q')
title('Angles');

figure()
plot(c2(1,:),c2(2,:))
xlabel('c2y')
ylabel('c2z')
title('Workspace');

figure()
plot(t,tau(:,1))
hold on
plot(t,tau(:,2))
hold off
legend('tau1','tau2')
xlabel('Time')
ylabel('Tau')
title('Torques')
