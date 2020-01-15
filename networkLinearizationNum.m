function [a,b] = networkLinearizationNum(w10, w1, w20, w2, na, nb, tau, x0, fun)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
a = zeros(1,na);
b = zeros(1,nb);
delta = 10^-5;

val = w20 + w2*tanh(w10 + w1*x0);
for i=tau:nb
    x = x0;
    x(i-tau+1) = x(i-tau+1) + delta;
    b(i) = (w20 + w2*tanh(w10 + w1*x)-val)/delta;
end
for i=1:na
    x = x0;
    x(nb-tau+1+i) = x(nb-tau+1+i) + delta;
    a(i) = -(w20 + w2*tanh(w10 + w1*x)-val)/delta;
end
end

