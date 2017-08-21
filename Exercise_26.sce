clc; clear; //clears the console and all previously stored variables

function [V0, CIl, CIr] = BS_EuCall_MC_IS (S0, r, sigma, K, T, mu, N, alpha)
    function V_ = f (Y)
        V_ = exp(-r*T-Y.*mu+mu^2/2).*max(S0*exp((r-sigma^2/2)*T+sigma*sqrt(T).*Y)-K, 0);
    endfunction
    V1(:) = f(grand(N, 1, "nor", mu, 1));
    V0 = mean(V1);
    var_V_is = variance(V1);
    // function to determine the value of the standard deviation of the mean
    // for a given alpha (probability), since std. nor. is symetric: q=(1-alpha)/2
    function std_dev = CI_std (alpha)
        std_dev = cdfnor("X", 0, 1, 1-(1-alpha)/2, (1-alpha)/2);
    endfunction    
    epsilon = CI_std(alpha)*sqrt(var_V_is/N);
    CIl=V0-epsilon;
    CIr=V0+epsilon;
endfunction

// Assigning given values to variables
S0=100; r=0.05; sigma=0.2; K=200; T=1; N=10000; alpha=0.95; mu=[0:0.1:6]

// For-loop over all mu
j=0 // auxillary index variable
CI_matrix = zeros(2,length(mu)) // Initializing matrix for the CI values
for i=mu
    [V0, CIl, CIr] = BS_EuCall_MC_IS (S0, r, sigma, K, T, i, N, alpha)
    j=j+1 // Increasing index variable by for each iteration of the for-loop
    V_plot(j) = V0 // Storing values for V0
    CI_matrix(1,j) = CIl // Storing values for left CI
    CI_matrix(2,j) = CIr // Storing values for right CI
end

// Plotting the V0 in dependence on mu
scf(0)
clf()
plot(mu, V_plot')
// It can be seen that mu should probably be choosen in the range of [2,5] as
// the variation of the price is lower than before or after this range

// Plotting the CI in dependence on mu to see if variance is really the lowest
// in the range of [2,5]
plot(mu,CI_matrix(1,:),"r")
plot(mu,CI_matrix(2,:),"r")
// Labeling the axis and adding a legend for better overview
legend(["Option Value"; string(alpha*100)+"% Confidence Intervals"])
ylabel("Value", "fontsize", 4, "color", "blue")
xlabel("mu", "fontsize", 4, "color", "blue")
// It is obviously the case, so mu should be chosen in the range of [3,4] as 
// the CI are the nearest in this region.
