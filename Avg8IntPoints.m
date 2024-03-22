function [E] = Avg8IntPoints(E,ini)

% Average the integration point strain and stress output data
switch ini.ElementType
    case 8 % This is for 8-node elements with reduced integration
        fprintf('No average of integration points neccesary.\n')
        
    case 20 % This is for 20-node elements with reduced integration
        % check the number of elements, data size, and integration points
        IntP = 8; % Number of integration points used in ABAQUS
        m = size(E,1)/IntP; % This should be equal to number of elements
        EE = zeros(m,6);    % Pre-allocation
        parfor n = 1:m
            EE(n,:) = mean( E((n-1)*IntP+1:n*IntP,:) ); %#ok
        end
        E = EE;
        fprintf('Integration points averaged for 20-node element.\n')
    otherwise
        error('Incorrect ini.ElementType definition')
end