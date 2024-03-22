function Rec=fil2Str(resultsFileName)
%
% Convert the data in the ABAQUS results file into a single row string
%
% Syntax
%     #Rec#=fil2Str(#resultsFileName#);
%
% Description
%     Convert the data contained in an ABAQUS results file in ASCII format
%     into a string that has one row.
%     The following option with parameter has to be specified in the ABAQUS
%     input file for the results file to be created:
%         *FILE FORMAT, ASCII
%
% Input parameters
%     #resultsFileName# (string) is a string containing the name of the
%         ABAQUS results file, along with its extension. The results file
%         is generated by Abaqus after the analysis has been completed.
%
% Output parameters
%     #Rec# ([1 x #m#]) is a string containing the data of the Abaqus
%         results file assembled in one row.
%
% _________________________________________________________________________
% Abaqus2Matlab - www.abaqus2matlab.com
% Copyright (c) 2017 by George Papazafeiropoulos
%
% If using this toolbox for research or industrial purposes, please cite:
% G. Papazafeiropoulos, M. Muniz-Calvente, E. Martinez-Paneda.
% Abaqus2Matlab: a suitable tool for finite element post-processing.
% Advances in Engineering Software. Vol 105. March 2017. Pages 9-16. (2017) 
% DOI:10.1016/j.advengsoft.2017.01.006



% Open the results file for reading
[fileID,errmsg]=fopen(resultsFileName,'r');
if fileID < 0
    error(errmsg)
end
% Read data from results file as a string and assign them to a cell array
% Concatenate each line without specifying delimiter, whitespace or end of
% line characters
try
    C=textscan (fileID, '%s', 'CollectOutput', '1', 'delimiter', ...
        '','whitespace','','endofline','');
catch
    C=textscan (fileID, '%s', 'CollectOutput', 1, 'delimiter', ...
        '','whitespace','','endofline','');
end
% Close the results file
fclose(fileID);
% Assign A
A=C{1}{1};
% Remove newline characters
A1=strrep(A,sprintf('\n'),'');
% Remove carriage return characters
Rec=strrep(A1,sprintf('\r'),'');


end

