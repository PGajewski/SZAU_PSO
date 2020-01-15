function [a,b] = networkLinearization(w10, w1, w20, w2, na, nb, tau, x0, fun)
%UNTITLED Summary of this function goes here
tan_deriv = @(x) (1 - tanh(x).^2);
a = zeros(1,na);
b = zeros(1,nb);
if strcmp(fun, 'tanh')
    for l=1:na
        a(l) =-sum(w2*(tan_deriv(w10 + w1*x0).*w1(:,nb-tau+1+l)));
    end
    for l=1:nb
        if l < tau
           b(l) = 0; 
        else
           b(l) = sum(w2*(tan_deriv(w10 + w1*x0).*w1(:,l-tau+1)));
       end
    end
end
end

