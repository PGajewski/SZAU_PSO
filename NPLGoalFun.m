function [J] = NPLGoalFun(x,y_zad,sim_time, w10, w1, w20, w2, na, nb, tau)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Process params.
global alfa1;
global alfa2;

global beta1;
global beta2;

global g1;
global g2;
global u_min;
global u_max;
global u_delay;

%neural network params.

lambda_limit = 0.01;

lambda = x(1);
N = round(x(2));
Nu = round(x(3));

% Coefficient limits.
if lambda < lambda_limit
    lambda = lambda_limit;
end    
if N < 1
    J = Inf;
    return;
end

if Nu < 1
    J = Inf;
    return;
end

if Nu > N
   J = Inf;
   return;
end

%Initialization.
y_vector = zeros(1,sim_time);
x_vector = zeros(2,sim_time);
u_vector = zeros(1,sim_time);

%Init state
x_vector(:,1) = [0;0];

lin_type = 'Anal'; %Num, Anal
%Main simulation loop.
for k=nb+1:sim_time
    %Object simulation.
    x_vector(1,k+1) = -alfa1 * x_vector(1,k) + x_vector(2,k) + beta1*g1(u_vector(1,k-u_delay));
    x_vector(2,k+1) = -alfa2 * x_vector(1,k) + beta2*g1(u_vector(1,k-u_delay));
    y_vector(1,k) = g2(x_vector(1,k)); %without noise
    %y_vector(1,k) = g2(x_vector(1,k)) + 0.02*(rand(1,1)-0.5); % with noise
    %y_vector(1,k) = w20 + w2*tanh(w10 + w1*[flip(u_vector(k-nb:k-tau)) flip(y_vector(k-na:k-1))]');


    %Create linearization point vector.
    x = [flip(u_vector(k-nb:k-tau)) flip(y_vector(k-na:k-1))]';

    %Count linearization coefficients.
    switch lin_type
        case 'Num'
            [a,b] = networkLinearizationNum(w10, w1, w20, w2, na, nb, tau, x, 'tanh');
        case 'Anal'
            [a,b] = networkLinearization(w10, w1, w20, w2, na, nb, tau, x, 'tanh');
    end
    %Count step response.
    s = zeros(N,1);
    for j=1:N
       for i=1:min(j,nb)
           s(j) = s(j) + b(i);
       end
       for i=1:min(j-1,na)
           s(j) = s(j) - a(i)*s(j-i);
       end
    end

    %Construct dynamic matrix.
    M = zeros(N,Nu);
    for j=1: Nu
       M(:,j) = [zeros(j-1,1); s(1:N-j+1)]';
    end

    K = (M' * M + lambda * eye(Nu))\(M');       
    
  
    %Construct free trajectory.
    y = w20 + w2*tanh(w10 + w1*[flip(u_vector(k-nb:k-tau)) flip(y_vector(k-na:k-1))]');
    dk = y_vector(1,k) - y;

    y0 = [y_vector(1:k-1) y ones(1,sim_time+N-k)*dk];
    u0 = [u_vector(1:k-1) ones(1,sim_time+N-k+1)*u_vector(k-1)];
    for j=1:N
        y0(k+j) = y0(k+j) + w20 + w2*tanh(w10 + w1*[flip(u0(k-nb+j:k-tau+j)) flip(y0(k-na+j:k-1+j))]');
    end

    y0 = y0(k+1:k+N)';
    
    %Count control signal.
    
    du = K*(y_zad(k)*ones(N,1) - y0);
    u_vector(k) = du(1) + u_vector(k-1);

    if u_vector(k) < u_min
        u_vector(k) = u_min;
    elseif u_vector(k) > u_max
        u_vector(k) = u_max;
    end
   
end

% Count error.
J = sum((y_zad - y_vector').^2);
end

