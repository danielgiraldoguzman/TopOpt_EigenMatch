function [Strain,eps] = StrainBrick20(xyz,u)
% Stress-strain calculations using the quadrature integration points
%
% INPUT:
%   xyz: nodal element coordinates
%   u: nodal element displacement

%% Gaussian Quadrature Points and weigths for numerical integration
% Gaussian Quadrature with n = 2
e = [-1/sqrt(3); 1/sqrt(3)];
w = [1; 1];
n = 2;

% Gaussian Quadrature with n = 3
% e = [sqrt(0.6); -sqrt(0.6); 0];
% w = [5/9; 5/9; 5/9];
% n = 3;

% Gaussian Quadrature with n = 4
% e = [sqrt(3 + sqrt(4.8))/7;
%     -sqrt(3 + sqrt(4.8))/7;
%     sqrt(3 - sqrt(4.8))/7;
%     -sqrt(3 - sqrt(4.8))/7];
% w = [1/2 - 1/(3*sqrt(4.8));
%     1/2 - 1/(3*sqrt(4.8));
%     1/2 + 1/(3*sqrt(4.8));
%     1/2 + 1/(3*sqrt(4.8))];
% n = 4;

%% Numerical evaluation
dof = 20*3;
B  = zeros(6,dof);   % Strain-displacement matrix
eps = zeros(6,n^3); % Strain vectors arrangement
pp = 0; % Pointer

% Extrapolation matrix for strain and stress to nodes
% a = (5 + sqrt(3))/4;
% b = -(sqrt(3) + 1)/4;
% c = (sqrt(3) - 1)/4;
% d = (5 - sqrt(3))/4;
% L = [a b c b b c d c;
%      b a b c c b c d;
%      c b a b d c b c;
%      b c b a c d c b;
%      b c d c a b c d;
%      c b c d b a b c;
%      d c b c c b a b;
%      c d c b d c b a];

for i=1:n
    for j=1:n
        for k=1:n
            % Integrations points
            ri = e(k);
            sj = e(j);
            tk = e(i);
            
            % Weights
            wx = w(k);
            wy = w(j);
            wz = w(i);
            
            % Shape Derivates Function with respect to r,s,t
            dNdr = [ ((sj - 1)*(tk - 1)*(2*ri + sj + tk + 1))/8;
                    -((sj - 1)*(tk - 1)*(sj - 2*ri + tk + 1))/8;
                    -((sj + 1)*(tk - 1)*(2*ri + sj - tk - 1))/8;
                    -((sj + 1)*(tk - 1)*(2*ri - sj + tk + 1))/8;
                    -((sj - 1)*(tk + 1)*(2*ri + sj - tk + 1))/8;
                    -((sj - 1)*(tk + 1)*(2*ri - sj + tk - 1))/8;
                     ((sj + 1)*(tk + 1)*(2*ri + sj + tk - 1))/8;
                     ((sj + 1)*(tk + 1)*(2*ri - sj - tk + 1))/8;
                                      -(ri*(sj - 1)*(tk - 1))/2;
                                 ((sj - 1)*(sj + 1)*(tk - 1))/4;
                                       (ri*(sj + 1)*(tk - 1))/2;
                                -((sj - 1)*(sj + 1)*(tk - 1))/4;
                                       (ri*(sj - 1)*(tk + 1))/2;
                                -((sj - 1)*(sj + 1)*(tk + 1))/4;
                                      -(ri*(sj + 1)*(tk + 1))/2;
                                 ((sj - 1)*(sj + 1)*(tk + 1))/4;
                                -((sj - 1)*(tk - 1)*(tk + 1))/4;
                                 ((sj - 1)*(tk - 1)*(tk + 1))/4;
                                -((sj + 1)*(tk - 1)*(tk + 1))/4;
                                 ((sj + 1)*(tk - 1)*(tk + 1))/4];
                             
            dNds = [ ((ri - 1)*(tk - 1)*(ri + 2*sj + tk + 1))/8;
                    -((ri + 1)*(tk - 1)*(2*sj - ri + tk + 1))/8;
                    -((ri + 1)*(tk - 1)*(ri + 2*sj - tk - 1))/8;
                    -((ri - 1)*(tk - 1)*(ri - 2*sj + tk + 1))/8;
                    -((ri - 1)*(tk + 1)*(ri + 2*sj - tk + 1))/8;
                    -((ri + 1)*(tk + 1)*(ri - 2*sj + tk - 1))/8;
                     ((ri + 1)*(tk + 1)*(ri + 2*sj + tk - 1))/8;
                     ((ri - 1)*(tk + 1)*(ri - 2*sj - tk + 1))/8;
                                -((ri - 1)*(ri + 1)*(tk - 1))/4;
                                       (sj*(ri + 1)*(tk - 1))/2;
                                 ((ri - 1)*(ri + 1)*(tk - 1))/4;
                                      -(sj*(ri - 1)*(tk - 1))/2;
                                 ((ri - 1)*(ri + 1)*(tk + 1))/4;
                                      -(sj*(ri + 1)*(tk + 1))/2;
                                -((ri - 1)*(ri + 1)*(tk + 1))/4;
                                       (sj*(ri - 1)*(tk + 1))/2;
                                -((ri - 1)*(tk - 1)*(tk + 1))/4;
                                 ((ri + 1)*(tk - 1)*(tk + 1))/4;
                                -((ri + 1)*(tk - 1)*(tk + 1))/4;
                                 ((ri - 1)*(tk - 1)*(tk + 1))/4];
                             
            dNdt = [ ((ri - 1)*(sj - 1)*(ri + sj + 2*tk + 1))/8;
                    -((ri + 1)*(sj - 1)*(sj - ri + 2*tk + 1))/8;
                    -((ri + 1)*(sj + 1)*(ri + sj - 2*tk - 1))/8;
                    -((ri - 1)*(sj + 1)*(ri - sj + 2*tk + 1))/8;
                    -((ri - 1)*(sj - 1)*(ri + sj - 2*tk + 1))/8;
                    -((ri + 1)*(sj - 1)*(ri - sj + 2*tk - 1))/8;
                     ((ri + 1)*(sj + 1)*(ri + sj + 2*tk - 1))/8;
                     ((ri - 1)*(sj + 1)*(ri - sj - 2*tk + 1))/8;
                                -((ri - 1)*(ri + 1)*(sj - 1))/4;
                                 ((ri + 1)*(sj - 1)*(sj + 1))/4;
                                 ((ri - 1)*(ri + 1)*(sj + 1))/4;
                                -((ri - 1)*(sj - 1)*(sj + 1))/4;
                                 ((ri - 1)*(ri + 1)*(sj - 1))/4;
                                -((ri + 1)*(sj - 1)*(sj + 1))/4;
                                -((ri - 1)*(ri + 1)*(sj + 1))/4;
                                 ((ri - 1)*(sj - 1)*(sj + 1))/4;
                                      -(tk*(ri - 1)*(sj - 1))/2;
                                       (tk*(ri + 1)*(sj - 1))/2;
                                      -(tk*(ri + 1)*(sj + 1))/2;
                                       (tk*(ri - 1)*(sj + 1))/2];


            % Jacobian
            Jacobian = [dNdr'; dNds'; dNdt']*xyz;
            
            %  Shape Function derivatives with respect to x,y,z
            dN = Jacobian\([dNdr'; dNds'; dNdt']);
            dNdx = dN(1,:);
            dNdy = dN(2,:);
            dNdz = dN(3,:);
            
            % Mechanical strain-displacement matrix
            B(1,1:3:dof) = dNdx;
            B(2,2:3:dof) = dNdy;
            B(3,3:3:dof) = dNdz;
            B(4,1:3:dof) = dNdy; B(4,2:3:dof) = dNdx;
            B(6,2:3:dof) = dNdz; B(6,3:3:dof) = dNdy; % Following ABAQUS
            B(5,1:3:dof) = dNdz; B(5,3:3:dof) = dNdx; % order for shear E.
            % NOTE: ABAQUS uses the following order:
            % xx, yy, zz, xy, xz, yz
            
            % Strain {eps}=[B]{u} for integartion point pp (i,j,k)
            pp = pp + 1; % Pointer
            eps(:,pp) = ( B*u(1:dof) )*wx*wy*wz;
        end
    end
end

% Average strain
Strain = mean(eps,2);