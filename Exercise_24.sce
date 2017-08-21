clc; clear; //clears the console and all previously stored variables

function [V0, epsilon] = BS_EuOption_MC_AV (S0, r, sigma, T, K, M, g)
    // Computing the 95% confidence interval for the european option price 
    // via the MC approach with antithetic variables:
    V_hat(:) = g(S0*exp((r-sigma^2/2)*T+sigma*grand(M, 1, "nor", 0, sqrt(T))))*exp(-r*T);
    V_(:) = g(S0*exp((r-sigma^2/2)*T+sigma*(-grand(M, 1, "nor", 0, sqrt(T)))))*exp(-r*T);
    V0 = 0.5*mean(V_hat+V_);
    // Sample variance formula is used, instead is it also possible to use
    // the build in "variance" function: variance((V_hat+V_)/2)
    var_AV = M/(M-2)*mean(((V_hat+V_)/2)^2-V0^2)
    epsilon = 1.96*sqrt(var_AV/M);
endfunction

// Function from ex. 21 is modified, so it is possible to compare the 
// radius epsolon of the 95% CI of both approaches
function [V0, epsilon] = EuOption_BS_MC (S0, r, sigma, T, K, M, g)
    // Computing the 95% confidence intervals for the european option price 
    // via the MC approach
    V_(:) = g(S0*exp((r-sigma^2/2)*T+sigma*grand(M, 1, "nor", 0, sqrt(T))))*exp(-r*T);
    V0 = mean(V_);
    var_hat = M/(M-1)*mean(V_^2-V0^2)
    epsilon =1.96*sqrt(var_hat/M);
endfunction

function V0 = BS_EuOption_Put(S0, r, sigma, T, K)
    // Integration formula
    function p = Phi(x)
       p = cdfnor("PQ", x, zeros(x), ones(x));
    endfunction
    d1 = ( log(S0/K) + (r+1/2*sigma^2)*T ) / ( sigma*sqrt(T) );
    d2 = d1 - sigma*sqrt(T);
    V0 = K*exp(-r*T)*Phi(-d2) - S0*Phi(-d1); 
endfunction

// Defining the european put function
function y = g(x)
    y = max(K - x, 0)
endfunction

// Assigning given values to variables
K=100; S0=100; r=0.05; sigma=0.2; T=1; M=100000;

// Calling the first function and displaying results
[V0, epsilon] = BS_EuOption_MC_AV (S0, r, sigma, T, K, M, g)
disp("Antithetic variables: V0 = "+string(V0)+ " Epsilon of the 95% CI: "+string(epsilon))

// Calling the second function and displaying results
[V0, epsilon] = EuOption_BS_MC (S0, r, sigma, T, K, M, g)
disp("Plain MC simulation:  V0 = "+string(V0)+ " Epsilon of the 95% CI: "+string(epsilon))

// Calling the third function and displaying results
V0 = BS_EuOption_Put(S0, r, sigma, T, K)
disp("BS EuPut exact price: V0 = "+string(V0))
