% This script construcs models based on gene expression data in
% condition sepcific manner

%% Initiate the COBRA toolbox and load the relvant modle
initCobraToolbox;
load Reduced_Recon2.mat
model = Reduced_Recon2;
mg = model.genes;

%% Import the RPKM data from the excel file

[NUMERIC,TXT,RAW]=xlsread('Z:\AQAI - Ashfaq Ali\GSM\PAlmitate_control\PalmDataReconFormat.xlsx');
t = TXT(2:size(TXT,1),2); % 1st column from data file ID names
n1 = NUMERIC(:,4:15); % cytokine_24h_mean
n2=n1;

%% Initialize variables.
filename = 'Z:\AQAI - Ashfaq Ali\GSM\Palmitate_control\Differentially_expressed_ENTREZ.csv';
delimiter = ',';
startRow = 2;

%% Format string for each line of text:
formatSpec = '%*q%*q%*q%f%q%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);
%% Allocate imported array to column variable names
Fold_change = dataArray{:, 1};
entrezgene = dataArray{:, 2};


%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% Define highly expressed, moderately expressed and lowly expressed genes i
% 1, 0, -1 
Threshold_low = median(n1);
%Threshold_high = 5;
Threshold_high =  quantile(n1, 0.75);

%Find the index of highly exressed genes in the data 


for i=1:length(n1)
    for j = 1:size(n1,2)
        if n1(i,j) >= Threshold_high(j)
            n2(i,j) = 1;
        elseif n1(i,j)>=Threshold_low(j) && n1(i,j)<Threshold_high(j)
            n2(i,j) = 0;
        else
            n2(i,j) = -1;
        end
    end
end

%% Creat extra variables based high, low and moderate expression for each
% condition using minimum 3 sample approach

med_cont = round(median(n2(:,(1:2:12)),2));
med_palm = round(median(n2(:,(2:2:12)),2));


%% creat extra condition specific variable based on statictical
% differences/fold changes provided by the collatorators
x = find(ismember(t,entrezgene)); 
med_palm_up = find(ismember(t,entrezgene(Fold_change>0)));
med_palm_down = find(ismember(t,entrezgene(Fold_change<0)));

med_palm_stat = med_palm;
med_cont_stat = med_cont;

med_palm_stat(med_palm_up) = 2;
med_palm_stat(med_palm_down) = 0;
n2 = [n2,med_cont,med_palm,med_palm_stat ];

%% Collect the output models and reactions using iMAT
%Perform Flux variability analyses on each of the model and collect minmum
%and maximum fluxes for each of the models

collect_max_fluxes ={size(n2,2)};
collect_min_fluxes = {size(n2,2)};
collect_reactions ={size(n2,2)};
Collect_models_palm = {size(n2,2)};
solver='Shlomi';

for i=1:size(n2,2)
    expressionData.Locus=(t);
    expressionData.Data = n2(:,i);
    [Collect_models_palm{i},collect_reactions{i}] = createTissueSpecificModel(model,expressionData , solver);
    [collect_min_fluxes{i},collect_max_fluxes{i}] = fluxVariability(Collect_models_palm{i});
     
end
%Save the data in excel format for further statitical analyses in R

save('Z:\AQAI - Ashfaq Ali\GSM\Palmitate_control\Palmitate_modelling091116.mat_2', 'collect_max_fluxes', 'collect_min_fluxes', 'collect_reactions','Collect_models_palm', 't', 'n2','entrezgene', 'Fold_change');  