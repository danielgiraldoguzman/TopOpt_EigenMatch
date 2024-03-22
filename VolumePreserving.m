function eta = VolumePreserving(rho_d,beta,ve,MaxVol,MinVol,Vol)
% Find eta parameter using the bi-section method

% Inital variables
a = 0;
b = 1;
TOL = 1e-3;
N = 1;
Before = rho_d'*ve; % Volume before the filter
Vol_d = Before/sum(ve);
if nargin == 5
    if Vol_d > MaxVol
        Before = sum(ve)*MaxVol;
    elseif Vol_d < MinVol
        Before = sum(ve)*MinVol;
    end
elseif nargin == 6
    switch Vol
        case 'MaxVol'
            Before = sum(ve)*MaxVol;
        case 'MinVol'
            Before = sum(ve)*MinVol;
    end
end

while N <= 100
    eta = (a + b)/2; % Mid point
    rho_h = HeavisideFilter(beta,eta,rho_d);
    
    % Before = rho_d'*ve; % Volume before the filter
    After = rho_h'*ve;  % Volume after the filter
    
    % Check if the problem converges
    if abs(After - Before)/Before < TOL && (b - a)/2 < TOL
       break
    elseif After > Before % More volume after the filter
        a = eta; % Move eta towards 1        
    elseif After < Before % Less volume after the filter
        b = eta; % Move eta towards 0    
    end
    N = N+1;
end

% End of Function
end