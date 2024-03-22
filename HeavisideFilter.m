function [rho_h,drho_h] = HeavisideFilter(beta,eta,rho)
% Volume preserving nonlinear density filter based on heaviside functions
% Shengli Xu Yuanwu Cai Gengdong Cheng
% 
% INPUT:
% beta      Parameter to smooth the Heaviside function
% eta       Parameter to control volume of filtered densities
% rho       Input pseudo-density 
% 
% OUTPUT:
% rho_h     Filtered density with the modified Heaviside function
% drho_h    Derivative of the modified density w.r.t. the filtered density

rho_h = zeros(length(rho),1);     % Pre-allocation
drho_h = zeros(length(rho),1);   % Pre-allocation

for e = 1:length(rho)
    rh = [];    % Clear temporary variables
    drh = [];   % Clear temporary variables
    ro = rho(e);
    
    if ro <= eta
        rh = eta*( exp(-beta*(1-ro/eta)) -(1-ro/eta)*exp(-beta) );
        drh = beta*exp( -beta*(1-ro/eta) ) + exp(-beta);
    elseif ro > eta
        rh = (1-eta)*( 1 - exp( -beta*(ro-eta)/(1-eta) ) ...
            + ( (ro-eta)*exp(-beta) )/(1-eta) ) + eta;
        drh = beta*exp( -beta*(ro-eta)/(1-eta) ) + exp(-beta);
    end
    
    rho_h(e) = rh;
    drho_h(e) = drh;
end

% Set density lower limit
limit = 1e-6;
rho_h(rho_h<limit) = limit;

end
