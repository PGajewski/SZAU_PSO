function [J] = PIDGoalFun(x,y_zad, T,sim_time)
%PIDGOALFUN Summary of this function goes here
%   Detailed explanation goes here
global alfa1;
global alfa2;

global beta1;
global beta2;

global g1;
global g2;
global u_min;
global u_max;
global u_delay;

K = x(1);
Ti = x(2);
Td = x(3);

%Initialization.
y_vector = zeros(1,sim_time);
x_vector = zeros(2,sim_time);
u_vector = zeros(1,sim_time);

%Init state
x_vector(:,1) = [0;0];

%Init PID
e_vector = zeros(1,sim_time);
r2 = K*Td/T;
r1 = K*(T/(2*Ti)-2*Td/T-1);
r0 = K*(1+T/(2*Ti)+Td/T);

%Main simulation loop.
for k=u_delay+1:sim_time
    %Object simulation.
    x_vector(1,k+1) = -alfa1 * x_vector(1,k) + x_vector(2,k) + beta1*g1(u_vector(1,k-u_delay));
    x_vector(2,k+1) = -alfa2 * x_vector(1,k) + beta2*g1(u_vector(1,k-u_delay));
    y_vector(1,k) = g2(x_vector(1,k)); %without noise
    %y_vector(1,k) = g2(x_vector(1,k)) + 0.02*(rand(1,1)-0.5); % with noise
    %y_vector(1,k) = w20 + w2*tanh(w10 + w1*[flip(u_vector(k-nb:k-tau)) flip(y_vector(k-na:k-1))]');
    tic

    %Count control signal.
    e_vector(k) = y_zad(k) - y_vector(k);
    u_vector(k) = r2*e_vector(k-2)+r1*e_vector(k-1)+r0*e_vector(k) + u_vector(k-1);

    if u_vector(k) < u_min
        u_vector(k) = u_min;
    elseif u_vector(k) > u_max
        u_vector(k) = u_max;
    end
end

J = sum((y_zad - y_vector').^2);
end

