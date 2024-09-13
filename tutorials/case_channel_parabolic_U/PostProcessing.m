%% "Read the files"
case1BinaryLog = readtable("/Volumes/OpenFOAM/multiformFoam/Case1/case1BinaryLog");
case1MembraneLog = readtable("/Volumes/OpenFOAM/multiformFoam/Case1/case1MembraneLog");
case1FluxLog = readtable("/Volumes/OpenFOAM/multiformFoam/Case1/case1FluxLog");

%%  % remove the following columns with this headings
case1MembraneLog = removevars(case1MembraneLog,{'Var1', 'Var2', 'Var3'});
case1FluxLog = removevars(case1FluxLog,{'Var1', 'Var2', 'Var3'});
case1BinaryLog = removevars(case1BinaryLog,{'Var1', 'Var2', 'Var3', 'Var4','Var6', 'Var7', 'Var9','Var10', 'Var12', 'Var13'});


%% Replace the existing table headings with this

case1BinaryLog.Properties.VariableNames = {'Time' 'Permeability' 'Porosity' 'DeltaP'};
case1MembraneLog.Properties.VariableNames = {'Conc Mem' 'Conc in' 'Average'};
case1FluxLog.Properties.VariableNames = {'Flux Mem' 'Flux in' 'Average'};


%% Find the unique rows with unique values in the first variable "Time"

C = unique(case1BinaryLog);
[C, ia] = unique(case1BinaryLog.Time);
B = case1BinaryLog(ia,:);


%% Contactenate the three tables

CaseTable = [B; case1MembraneLog; case1FluxLog ];


