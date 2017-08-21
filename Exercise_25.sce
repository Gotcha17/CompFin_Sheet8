clc; clear; //clears the console and all previously stored variables

function V0 = BS_EuOption_MC_CV (S0, r, sigma, T, K, M, g)
    // First MC simulation to estimate optimal value for Beta
    V_T(:) = g(S0*exp((r-sigma^2/2)*T+sigma*grand(M, 1, "nor", 0, sqrt(T))));
    S_T(:) = S0*exp((r-sigma^2/2)*T+sigma*grand(M, 1, "nor", 0, sqrt(T)));
    // E(S_T) is known, as it is the equivalent of investing the value of the
    // stock S0 for period of time T at interest rate r
    Cov = mean((V_T - mean(V_T)) .* (S_T - mean(S_T)));
    Beta = Cov / variance(S_T);
    // Now we just plug in the values in the formula for control variables
    // S_T is multiplied with previously simulated optimal Beta to reduce variance
    // Finally value of the option is discounted to the present value
    V0 = exp(-r*T) * (mean(V_T) - mean(Beta * S_T) + Beta * S0)
endfunction

// Calculation of exacte option price with Black Scholes formula for the
// European Call
function V0 = BS_EuOption_Call(S0, r, sigma, T, K)
    // Integration formula
    function p = Phi(x)
       p = cdfnor("PQ", x, zeros(x), ones(x));
    endfunction
    d1 = ( log(S0/K) + (r+1/2*sigma^2)*T ) / ( sigma*sqrt(T) );
    d2 = d1 - sigma*sqrt(T);
    V0 = S0*Phi(d1) - K*exp(-r*T)*Phi(d2); 
endfunction

// Defining the european call function
function y = g(x)
    y = max(x - K, 0)
endfunction

// Assigning given values to variables
K=100; S0=120; r=0.05; sigma=0.2; T=1; M=100000;

// Calling the first function and displaying results
V0 = BS_EuOption_MC_CV (S0, r, sigma, T, K, M, g)
disp("Control variates: V0 = "+string(V0))

// Calling the second function and displaying results
V0 = BS_EuOption_Call(S0, r, sigma, T, K)
disp("Exact BS Eu_Call : V0 = "+string(V0))
