# -*- coding: mbcs -*-
# (taken from abaqus *.rec and *.rpy files)
# or:
# -*- coding: utf-8 -*-
# Necessary import statements
import os
from odbAccess import *
# In addition, if your script refers to any of the Symbolic Constants defined
# in the Abaqus Scripting Interface, your script must contain the following
# statement:
from abaqusConstants import *
#==========================================================================
# INPUT VARIABLES
# odbName is a string containing the name of the odb file which is to be postprocessed
# along with its extension
odbName='FreqResZ_1372.odb'
# Create a lock file to notify that the function is running
lckFile=open(odbName[0:-4]+'.lck','w')
try:
    # stepName can be either a string containing the name of the step to be postprocessed
    # or (None) in which case all the steps contained in the odb file will be postprocessed.
    # Otherwise, no postprocessing takes place and an error message is issued.
    stepName=None
    # outputVar can be either a string containing the name of the history output identifier
    # to be postprocessed, or (None) in which case all the history output data contained
    # in the odb file will be postprocessed. Otherwise, postprocessing for the specific step
    # and history region is skipped and continues to the next, while a warning message
    # is issued.
    outputVar=None
    # indOut is a boolean which defines if indexing data of the Abaqus results will be
    # printed in a separate ind file
    indOut=True
    #==========================================================================
    # MAIN FUNCTION DEFINITION
    def getHistoryOutput(odbName,stepName,outputVar,indOut):
        # Open the Abaqus odb file only with reading permission
        odb1=openOdb(path=odbName, readOnly=True)
        # Open the output file in which the results will be printed
        # Delete any existing file with the name (NameOfFile) without asking for
        # permission (option 'w' in function open)
        NameOfFile=odbName[0:-4]
        outFile=open(NameOfFile+'.out','w')
        #==========================================================================
        # Check if the odb contains any steps
        try:
            # If the odb contains any steps find their number
            outSteps=odb1.steps.values()
            numSteps=len(odb1.steps.keys())
        except:
            # If the odb does not contain any steps, issue an error in the output
            # file and terminate the postprocessing
            out='Error: No step exists in the output database %s' \
                % (odbName)
            print(out)
            outFile.close()
            return
        # If the number of steps are zero, issue an error in the output file and
        # terminate the postprocessing
        if numSteps==0:
            out='Error: No step exists in the output database %s' \
                % (odbName)
            print(out)
            outFile.close()
            return
        # Check to see if the step (stepName) exists in the odb 
        if stepName:
            try:
                # If the step (stepName) exists in the odb only one step (specStep) will
                # be postprocessed
                specStep=odb1.steps[stepName]
                outSteps=[];
                outSteps.append(specStep)
                numSteps=1
            except:
                # If the step (stepName) does not exist in the odb, issue a warning
                # in the output file and postprocess all steps
                out='Warning: The step %s does not exist in the output database %s. ' \
                    'All steps will be postprocessed.' \
                    % (stepName, odbName)
                print(out)
        #==========================================================================
        # If a valid outputVar is given, post process only the corresponding history
        # output result. If not, issue a warning and postprocess all history output
        # results.
        if outputVar:
            allVars=False
        else:
            allVars=True
            out='Warning: All history output variables will be printed.'
            print(out)
        #==========================================================================
        # Print Abaqus results to the *.out file
        # Look in '34.8 HistoryOutput object' of Abaqus Documentation - Scripting
        # reference manual
        # Loop over steps:
        for step in outSteps:
            allHistoryRegions=step.historyRegions
            # Loop over history regions:
            for key in allHistoryRegions.keys():
                outputAll=allHistoryRegions[key]
                # Check first if all history output results are to be written to
                # the output file (*.out)
                if allVars:
                    outputSet=outputAll
                    # Loop over history outputs:
                    for values in outputSet.historyOutputs.values():
                        # Loop over history output data:
                        for data in values.data:
                            for item in data:
                                outFile.write(str(item)+'\n')
                # Else if the specified output identifier (outputVar) exists in
                # historyOutputs, print only the desired results.
                elif outputAll.historyOutputs.has_key(outputVar):
                    outputSet=outputAll.getSubset(variableName=outputVar)
                    # Loop over history outputs:
                    for values in outputSet.historyOutputs.values():
                        # Loop over history output data:
                        for data in values.data:
                            for item in data:
                                outFile.write(str(item)+'\n')
                # Else issue a warning and do not postprocess anything
                else:
                    out='Warning: In step %s, history region %s, there is no ' \
                        'history output with identifier %s.' \
                        % (step.name, key, outputVar)
                    print(out)
        # Flush the internal buffer
        outFile.flush()
        # Close the output file
        outFile.close()
        #==========================================================================
        if indOut:
            # Print indexing data for Abaqus results to the *.ind file
            # Open the data file in which the indexing data for the results will be
            # printed
            # Delete any existing file with the name (NameOfFile) without asking for
            # permission (option 'w' in function open)
            datFile=open(NameOfFile+'.ind','w')
            i1=0
            # Loop over steps:
            for step in outSteps:
                allHistoryRegions=step.historyRegions
                i1=i1+1
                i2=0
                # Loop over history regions:
                for key in allHistoryRegions.keys():
                    outputAll=allHistoryRegions[key]
                    # Check first if all history output results are to be written
                    # to the output file (*.out)
                    if allVars:
                        outputSet=outputAll
                        i2=i2+1
                        i3=0
                        # Loop over history outputs:
                        for values in outputSet.historyOutputs.values():
                            i3=i3+1
                            i4=0
                            # Loop over history output data:
                            for data in values.data:
                                i4=i4+1
                                for item in data:
                                    datFile.write(str(i1)+','+str(i2)+','+str(i3)+','+str(i4)+'\n')
                    # Else if the specified output identifier (outputVar) exists in
                    # historyOutputs, print only the desired results.
                    elif outputAll.historyOutputs.has_key(outputVar):
                        outputSet=outputAll.getSubset(variableName=outputVar)
                        i2=i2+1
                        i3=0
                        # Loop over history outputs:
                        for values in outputSet.historyOutputs.values():
                            i3=i3+1
                            i4=0
                            # Loop over history output data:
                            for data in values.data:
                                i4=i4+1
                                for item in data:
                                    datFile.write(str(i1)+','+str(i2)+','+str(i3)+','+str(i4)+'\n')
            # Flush the internal buffer
            datFile.flush()
            # Close the output file
            datFile.close()
        #==========================================================================
        # Close the output database before exiting the program
        odb1.close()
    #==========================================================================
    # POST PROCESSING OF RESULTS
    getHistoryOutput(odbName,stepName,outputVar,indOut)
except:
    print('A2MOut '+'Error in execution of Python script')
# Delete the lock file after the function has terminated
lckFile.close()
os.remove(odbName[0:-4]+'.lck')
# _________________________________________________________________________
# Abaqus2Matlab - www.abaqus2matlab.com
# Copyright (c) 2017 by George Papazafeiropoulos
#
# If using this toolbox for research or industrial purposes, please cite:
# G. Papazafeiropoulos, M. Muniz-Calvente, E. Martinez-Paneda.
# Abaqus2Matlab: a suitable tool for finite element post-processing.
# Advances in Engineering Software. Vol 105. March 2017. Pages 9-16. (2017) 
# DOI:10.1016/j.advengsoft.2017.01.006
