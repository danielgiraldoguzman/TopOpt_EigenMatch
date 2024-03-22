function [U,S,E,EigF] = ...
    ReadFILResults(Inp_file,recordKeys,ElementType,NumElements)

filFile=[ Inp_file '.fil'];
% ini.recordKeys = [101 11 21 1980]; % 101-U, 11-S, 21-E, 1980->EigF
out = readFil(filFile,recordKeys);

% out{1,1} - displacements (key 101) data order
% Column 1   | Col 2 | Col 3 | Col 4 | Col 5 | Col 6 | Col 7  |
% Node index | X, Y, Z real part     | X, Y, Z imaginary part |

% out{1,2} - stress (key 11) data order
% Col 1 | Col 2 | Col 3 | Col 4 | Col 5 | Col 6 |
% ------------------ real part ------------------
% S11   | S22   | S33   | S12   | S13   | S23   |
%
% Col 7 | Col 8 | Col 9 | Col 10 | Col 11 | Col 12 |
% ----------------- imaginary part -----------------
% S11   | S22   | S33   | S12    | S13    | S23    |

% out{1,3} - strain (key 21) data order
% Col 1 | Col 2 | Col 3 | Col 4 | Col 5 | Col 6 |
% ------------------ real part ------------------
% E11   | E22   | E33   | E12   | E13   | E23   |
%
% Col 7 | Col 8 | Col 9 | Col 10 | Col 11 | Col 12 |
% ----------------- imaginary part -----------------
% E11   | E22   | E33   | E12    | E13    | E23    |

% out{1,4} - Modal solution (key 1980) data order
% Col 1             | Col 2      | Col 3            | Col 4             |
% Eigenvalue number | Eigenvalue | Generalized mass | Composite damping |
% -----------------------------------------------------------------------
% Participation factor per Mode #
% Col 5       | Col 7   | Col 9   | Col 11     | Col 13 | Col 15
% x-component | y-comp. | z-comp. | x-rotation | y-rot. | z-rot.
% -----------------------------------------------------------------------
% Effective mass per Mode #
% Col 6       | Col 8   | Col 10  | Col 12     | Col 14 | Col 16
% x-component | y-comp. | z-comp. | x-rotation | y-rot. | z-rot.

%% Find record keys index
idxDisp =   recordKeys == 101;
idxStress = recordKeys == 11;
idxStrain = recordKeys == 21;
idxModal =  recordKeys == 1980;

%% Check for displacement data (including eigenvectors)
if any(idxDisp) && any(out{1,idxDisp},'all')
% Pre-allocation
NumFreq = length(out{1,idxDisp})/max(out{1,idxDisp}(:,1));
Set1_size = size(out{1,idxDisp},1)/NumFreq; % output length
U = zeros(Set1_size,3,NumFreq); % 3 is for Ux,Uy,Uz
for n = 1:NumFreq
    idx0 = Set1_size*(n-1)+1 : Set1_size*n;
    U(:,1,n) = out{1,idxDisp}(idx0,2) + 1i*out{1,idxDisp}(idx0,5);
    U(:,2,n) = out{1,idxDisp}(idx0,3) + 1i*out{1,idxDisp}(idx0,6);
    U(:,3,n) = out{1,idxDisp}(idx0,4) + 1i*out{1,idxDisp}(idx0,7);
end

else
    U = [];
end

%% Check for Stress data
if any(idxStress) && any(out{1,idxStress},'all')
if ElementType == 20; IntP = 8; elseif ElementType == 8; IntP = 1; end
% Pre-allocation
% NumFreq = length(out{1,idxStress})/max(out{1,idxStress}(:,1));
Set2_size = size(out{1,idxStress},1); % length for stress output
NumFreq = (Set2_size/IntP)/NumElements;
S = zeros(Set2_size,6,NumFreq); % 6 is for Sxx,Syy,Szz,Sxy,Sxz,Syz
for n = 1:NumFreq
    idx1 = Set2_size*(n-1)+1 : Set2_size*n;
    S(:,1,n) = out{1,idxStress}(idx1,1) + 1i*out{1,idxStress}(idx1,7);
    S(:,2,n) = out{1,idxStress}(idx1,2) + 1i*out{1,idxStress}(idx1,8);
    S(:,3,n) = out{1,idxStress}(idx1,3) + 1i*out{1,idxStress}(idx1,9);
    S(:,4,n) = out{1,idxStress}(idx1,4) + 1i*out{1,idxStress}(idx1,10);
    S(:,5,n) = out{1,idxStress}(idx1,5) + 1i*out{1,idxStress}(idx1,11);
    S(:,6,n) = out{1,idxStress}(idx1,6) + 1i*out{1,idxStress}(idx1,12);
end

else
    S = [];
end

%% Check for Strain data
if any(idxStrain) && any(out{1,idxStrain},'all')
if ElementType == 20; IntP = 8; elseif ElementType == 8; IntP = 1; end
% Pre-allocation
% NumFreq = length(out{1,idxStrain})/max(out{1,idxStrain}(:,1));
Set3_size = size(out{1,idxStrain},1); % length for stress output
NumFreq = (Set3_size/IntP)/NumElements;
E = zeros(Set3_size,6,NumFreq); % 6 is for Exx,Eyy,Ezz,Exy,Exz,Eyz
for n = 1:NumFreq
    E(:,1,n) = out{1,idxStrain}(idx1,1) + 1i*out{1,idxStrain}(idx1,7);
    E(:,2,n) = out{1,idxStrain}(idx1,2) + 1i*out{1,idxStrain}(idx1,8);
    E(:,3,n) = out{1,idxStrain}(idx1,3) + 1i*out{1,idxStrain}(idx1,9);
    E(:,4,n) = out{1,idxStrain}(idx1,4) + 1i*out{1,idxStrain}(idx1,10);
    E(:,5,n) = out{1,idxStrain}(idx1,5) + 1i*out{1,idxStrain}(idx1,11);
    E(:,6,n) = out{1,idxStrain}(idx1,6) + 1i*out{1,idxStrain}(idx1,12);
end

else
    E = [];
end

%% Check for eigenfrequency data
if any(idxModal) && any( cell2mat(out{1,idxModal}),'all' )
    lambda = cell2mat( out{1,idxModal}(:,2) );
    lambda(lambda<0) = 0;
    EigF = sqrt(lambda)/(2*pi);
    
else
    EigF = [];
end

% End of function
end