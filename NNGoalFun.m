function [J] =NNGoalFun(x, K, na, nb, tau, u_learning, y_learning)
%NNGoalFun Summary of this function goes here
%   Detailed explanation goes here

%Parse data from x vector into NN wages.
w10 = x(1);
w1 = zeros(K,nb-tau+1+na);
for i=1:K
    %2+(i-1)*(nb-tau+1+na)
    %2+i*(nb-tau+1+na)-1
    w1(i,:) = x(2+(i-1)*(nb-tau+1+na):(2+i*(nb-tau+1+na)-1))';
end

w20 = x(K*(nb-tau+1+na)+2);
w2 = x(K*(nb-tau+1+na)+3:K*(nb-tau+1+na)+2+K)';

%Count response of network for learning data.
y_vector = zeros(1,length(u_learning));
y_vector(1:nb) = y_learning(1:nb);

for i=(nb)+1:length(u_learning)
    y_vector(1,i) = w20 + w2*tanh(w10 + w1*[flip(u_learning(i-nb:i-tau)) flip(y_vector(i-na:i-1))]');
end

%Count error.
J = sum((y_learning - y_vector).^2);
