%% PSO in academic problem.

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

% PSO made by Wael Korani.
n = 150;          % Size of the swarm " no of birds "
bird_setp  = 100; % Maximum number of "birds steps"

c2 =1.4;          % PSO parameter C1 
c1 = 0.12;        % PSO parameter C2 
w =0.9;           % pso momentum or inertia

[best_position,Jbest_min] = PS0Function(him_fun,n,bird_setp, c1,c2,w);

%% Process
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

