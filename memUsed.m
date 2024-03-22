function memUsed
usr = memory;
y = usr.MemUsedMATLAB/1e9; % 1e9-> GigaBytes
fprintf('Memory used: %.5g GigaBytes \n', y)
% only works on WINDOWS