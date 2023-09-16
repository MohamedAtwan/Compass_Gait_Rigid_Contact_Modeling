function f = trajectoryCostFun(b)
global N k m mH g T gamma Q
f = 0;
t = linspace(0,T,N);
% lambda = [1; 1]*1;
bm = b;
b = zeros(2,k+1);
b(1,:) = modifyPolynomial(bm(1,1:k+1),k);
b(2,:) = modifyPolynomial(bm(1,k+2:end),k);
for tm = 1:N
    ts = ones(k+1,1);
    for i = k:-1:0
        ts(k-i+1) = t(tm)^i;
    end
    q = b*ts;
    bd = [modifyPolynomial(polyder(b(1,:)),k); modifyPolynomial(polyder(b(2,:)),k)];
    qd = bd*ts;
    qdd = [modifyPolynomial(polyder(bd(1,:)),k); modifyPolynomial(polyder(bd(2,:)),k)]*ts;
%     c = forwardKinematicsLeft([0;0;q]);
%     if c(3) < 0.0
% 
%         dvdq = zeros(2,1);
%         % dvdq(1) = -m*g*0.5*sin(q(1)-gamma) - mH*g*sin(q(1)-gamma) - m*g*sin(q(1)-gamma) + 0.5*m*g*sin(q(1)+q(2)-gamma);
%         dvdq(1) = -((3/2)*m*g+mH*g)*sin(q(1)-gamma)+(1/2)*m*g*sin(q(1)+q(2)-gamma);
%         dvdq(2) = m*g*0.5*sin(q(1)+q(2)-gamma);
%         A = [cos(q(1)+q(2))-cos(q(1)) cos(q(1)+q(2));
%              sin(q(1)+q(2))-sin(q(1)) sin(q(1)+q(2))];
%         tau = -A*lambda+dvdq;
%     else
    tau = inverseDribbelDynamics(q,qd,qdd);
%     end

    f = f + (tau'*Q*tau);
f = (T/N)*f;
end
