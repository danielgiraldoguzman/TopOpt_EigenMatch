function Freq = ReadFreq(Inp_file)

% Read all Histoy Output data
% Pyton scripts are currently set up for 'TopOpt' and 'Step-2'.
% To modify these parameters, go to the Python file readAllHistoryOdb.py
odb_name= [Inp_file '.odb'];
[odbOut,odbDat,~] = readAllHistoryOdb(odb_name);

% fprintf('\nOutput database created.\n')

% Get size of indexing data
Sets = max(odbDat(:,2));
SubSets = max(odbDat(:,3));
Steps = max(odbDat(:,4));

% Retrieve index to distribute data in array format
r = odbDat(:,2);
n = odbDat(:,3);
s = odbDat(:,4);
q = repmat([1;2],Steps*Sets*SubSets,1);
SIZ = [Steps 2 Sets SubSets];
IND = sub2ind(SIZ,s,q,r,n); % Generate array-index based on odbData idexing

% Populate data array
Dat = zeros(Steps,2,Sets,SubSets);
Dat(IND) = odbOut;

% Extract frequency vector
Freq = Dat(:,1,1,1);

end