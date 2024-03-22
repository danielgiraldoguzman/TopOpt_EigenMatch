%% Starts the TOPslicer GUI
fprintf('Starting TOPslicer v0.8.4 -- Tomas Zegard\n')

% This call starts TOPslicer with an empty input file field and the default
% output file prefix.
TOPslicer

% Use the following call instead to start TOPslicer and set the input
% filename and output file prefix to some user defined strings.
% TOPslicer('Input','example2.mat','outputprefix','example2_output')