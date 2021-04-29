 %% "Read the files"
dataBinary = importdata("/Volumes/OpenFOAM/multiformFoam/Case1/case1BinaryLog");
dataFlux = importdata("/Volumes/OpenFOAM/multiformFoam/Case1/case1FluxConcLog");
 
%% Find the unique rows with unique values in the first variable "Time"

[~,b]=unique(dataBinary(:,1),'stable');
output=dataBinary(b,:);


