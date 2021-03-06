%% PSO in academic problem.
clear all;
% Himmelblau function
him_fun = @(x) (x(1)^2+x(2)-11)^2 + (x(1)+x(2)^2-7)^2;

% Function plot.
x1_vector = -5:0.1:5;
x2_vector = -5:0.1:5;
y_vector = zeros(length(x1_vector),length(x2_vector));

for i=1:length(x1_vector)
    for j=1:length(x2_vector)
        y_vector(i,j) = him_fun([x1_vector(i) x2_vector(j)]);
    end
end

figure;
surf(x1_vector, x2_vector, y_vector);
title('Himmelblau function');
xlabel('x_1');
ylabel('x_2');
zlabel('y');
minimal = fmincon(him_fun, [0, 0], [], [], [], [], [-5; -5], [5, 5]);
% PSO made by Wael Korani.
n = 150;          % Size of the swarm " no of birds "
bird_setp  = 100; % Maximum number of "birds steps"

c1 =1.6;          % PSO parameter C1 
c2 = 1.1;        % PSO parameter C2 
w =0.1;           % pso momentum or inertia
figure;
[best_position,Jbest_min] = PS0Function(him_fun,2,n,bird_setp, c1,c2,w);
title(sprintf('PSO for Himmelblau, c1 = %g, c2 = %g, w = %g', c1, c2, w));
xlabel('x_1');
ylabel('x_2');
%% Process
global alfa1;
global alfa2;
global beta1;
global beta2;

global g1;
global g2;
global u_min;
global u_max;
global u_delay;

alfa1 = -1.489028;
alfa2 = 0.535261;
beta1 = 0.012757;
beta2 = 0.010360;
u_min = -1;
u_max = 1;
u_delay = 5;

g1 = @(u)((exp(7*u)-1)/(exp(7*u)+1));
g2 = @(x)(0.25*(1-exp(-2.5*x)));

%% Goal y trajectory.

sim_time = 1000;
%Generate changing set values vector.
steps = [0,  0;
        100, 0.1;
        200, 0;
        300, 0.05;
        400, 0.07;
        500, -0.01;
        600, 0.02;
        700, 0.1;
        800, -0.03;
        900, 0.01
        1000, 0];
    
y_zad = zeros(sim_time,1);
for i=1:sim_time
    for j=2:length(steps)
        if steps(j,1) > i
            y_zad(i) = steps(j-1,2);
           break; 
        end
    end
end


%% PID regulator.
T = 0.1;

pid_fun = @(x) PIDGoalFun(x,y_zad,T,sim_time);

%PSO params.
n = 150;          % Size of the swarm " no of birds "
bird_setp  = 100; % Maximum number of "birds steps"

c1 =1.5;          % PSO parameter C1 
c2 = 0.5;        % PSO parameter C2 
w =0.5;           % pso momentum or inertia

[best_coef,Jbest_min] = PS0Function(pid_fun,3,n,bird_setp, c1,c2,w);

%Check best regulator.
K = best_coef(1);
Ti = best_coef(2);
Td = best_coef(3);

% Method engeneer
%K = 2.5;
%Ti = 1.3;
%Td = 0.15;

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

figure;
stairs(1:sim_time, y_vector);
hold on;
stairs(1:sim_time, y_zad);
hold off;
xlabel('k');
ylabel('y');
title('Przebiegi najlepszego regulatora PID');
legend('y','y_{zad}');

figure;
stairs(1:sim_time, u_vector);
title('Przebiegi najlepszego regulatora PID');
xlabel('k');
ylabel('u');
legend('u');

%% NPL regulator
[tau, nb, na, K, max_iter, error, algorithm] = readConfig();

load('bestmodel','w10','w1','w20','w2');

npl_fun = @(x) NPLGoalFun(x,y_zad,sim_time, w10, w1, w20, w2, na, nb, tau);

%PSO params.
n = 150;          % Size of the swarm " no of birds "
bird_setp  = 100; % Maximum number of "birds steps"

c2 =1;          % PSO parameter C1 
c1 = 0.7;        % PSO parameter C2 
w =0.6;           % pso momentum or inertia

[best_coef,Jbest_min] = PS0Function(npl_fun,3,n,bird_setp, c1,c2,w);

%Check best regulator.
lambda_limit = 0.01;
x = best_coef;
lambda = x(1);
N = round(x(2));
Nu = round(x(3));

%Initialization.
y_vector = zeros(1,sim_time);
x_vector = zeros(2,sim_time);
u_vector = zeros(1,sim_time);

%Init state
x_vector(:,1) = [0;0];

lin_type = 'Anal'; %Num, Anal
%Main simulation loop.
for k=nb+1:sim_time
    k
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

figure;
stairs(1:sim_time, y_vector);
hold on;
stairs(1:sim_time, y_zad);
hold off;
xlabel('k');
ylabel('y');
title('Przebiegi najlepszego regulatora NPL');
legend('y','y_{zad}');

figure;
stairs(1:sim_time, u_vector);
title('Przebiegi najlepszego regulatora NPL');
xlabel('k');
ylabel('u');
legend('u');

%% Neural Network
[tau, nb, na, K, max_iter, error, algorithm] = readConfig();

[u_learning, y_learning] = readData('dane.txt');

nn_fun = @(x) NNGoalFun(x, K, na, nb, tau, u_learning, y_learning);

%PSO params.
n = 150;          % Size of the swarm " no of birds "
bird_setp  = 100; % Maximum number of "birds steps"

c1 =1.4;          % PSO parameter C1 
c2 = 0.12;        % PSO parameter C2 
w =0.9;           % pso momentum or inertia

[best_coef,Jbest_min] = PS0Function(nn_fun,K*(nb-tau+1+na)+2+K,n,bird_setp, c1,c2,w);

%Check best network.
%Parse data from x vector into NN wages.
w10 = best_coef(1);
w1 = zeros(K,nb-tau+1+na);
for i=1:K
    %2+(i-1)*(nb-tau+1+na)
    %2+i*(nb-tau+1+na)-1
    w1(i,:) = best_coef(2+(i-1)*(nb-tau+1+na):(2+i*(nb-tau+1+na)-1))';
end

w20 = best_coef(K*(nb-tau+1+na)+2);
w2 = best_coef(K*(nb-tau+1+na)+3:K*(nb-tau+1+na)+2+K)';

%Count response of network for learning data.
y_vector = zeros(1,length(u_learning));
y_vector(1:nb) = y_learning(1:nb);

for i=(nb)+1:length(u_learning)
    y_vector(1,i) = w20 + w2*tanh(w10 + w1*[flip(u_learning(i-nb:i-tau)) flip(y_vector(i-na:i-1))]');
end

%Compare model for learning data.
figure;
hold on;
plot(1:length(y_vector),y_vector);
plot(1:length(y_learning),y_learning);
plot(1:length(u_learning),u_learning);
hold off;
title('Symulacja modelu neuronowego dla danych ucz.');
legend('y_{mod}','y_{learning}','u');
xlabel('t');
ylabel('y,u');

figure;
h = scatter(y_learning,y_vector,0.1);
h.Marker='.';
title('Relacja wyj�cia procesu i modelu, dane ucz.');
xlabel('y_{learning}');
ylabel('y_{mod}');

%Validate model.
[u_val, y_val] = readData('dane_wer.txt');

%Count response of network for validating data.
y_vector = zeros(1,length(u_val));
y_vector(1:nb) = y_val(1:nb);

for i=(nb)+1:length(u_val)
    y_vector(1,i) = w20 + w2*tanh(w10 + w1*[flip(u_val(i-nb:i-tau)) flip(y_vector(i-na:i-1))]');
end

%Compare model for validating data.
figure;
hold on;
plot(1:length(y_vector),y_vector);
plot(1:length(y_val),y_val);
plot(1:length(u_val),u_val);
hold off;
title('Symulacja modelu neuronowego dla danych ver.');
legend('y_{mod}','y_{val}','u');
xlabel('t');
ylabel('y,u');
errorr = (y_val-y_vector)*(y_val-y_vector)'
figure;
h = scatter(y_val,y_vector,0.1);
h.Marker='.';
title('Relacja wyj�cia procesu i modelu, dane ver.');
xlabel('y_{val}');
ylabel('y_{mod}');



%% Compare PSO model with experimental method.
com_choose = 'PID'; %PID, NPL, NN

switch com_choose
    case 'PID'
        load('best_project_2','K','Ti','Td');
        
        %Initialization.
        y_check_vector = zeros(1,sim_time);
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
            y_check_vector(1,k) = g2(x_vector(1,k)); %without noise
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
    case 'NPL'
        load('best_project_2','lambda','N','Nu');
        
        %Initialization.
        y_vector = zeros(1,sim_time);
        x_vector = zeros(2,sim_time);
        u_vector = zeros(1,sim_time);

        %Init state
        x_vector(:,1) = [0;0];

        lin_type = 'Anal'; %Num, Anal
        %Main simulation loop.
        for k=nb+1:sim_time
            k
            %Object simulation.
            x_vector(1,k+1) = -alfa1 * x_vector(1,k) + x_vector(2,k) + beta1*g1(u_vector(1,k-u_delay));
            x_vector(2,k+1) = -alfa2 * x_vector(1,k) + beta2*g1(u_vector(1,k-u_delay));
            y_check_vector(1,k) = g2(x_vector(1,k)); %without noise
            %y_vector(1,k) = g2(x_vector(1,k)) + 0.02*(rand(1,1)-0.5); % with noise
            %y_vector(1,k) = w20 + w2*tanh(w10 + w1*[flip(u_vector(k-nb:k-tau)) flip(y_vector(k-na:k-1))]');


            %Create linearization point vector.
            x = [flip(u_vector(k-nb:k-tau)) flip(y_check_vector(k-na:k-1))]';

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
            y = w20 + w2*tanh(w10 + w1*[flip(u_vector(k-nb:k-tau)) flip(y_check_vector(k-na:k-1))]');
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
        
    case 'NN'
        load('best_project_2','w1','w10','w2','w20');
        
        %Count response of network for validating data.
        y_check_vector = zeros(1,length(u_val));
        y_check_vector(1:nb) = y_val(1:nb);

        for i=(nb)+1:length(u_val)
            y_check_vector(1,i) = w20 + w2*tanh(w10 + w1*[flip(u_val(i-nb:i-tau)) flip(y_check_vector(i-na:i-1))]');
        end
end

%Plot results.
figure;
hold on;
plot(1:length(y_vector),y_vector);
plot(1:length(y_check_vector),y_check_vector);
hold off;
title(sprintf('PSO mode vs. %s',com_choose));
legend('y_{PSO}','y_{exp}');
xlabel('t');
ylabel('y');

