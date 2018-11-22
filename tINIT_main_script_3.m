% content = Cheng_readTXT('PalmDataReconFormat.csv',',');
% Genes = content(2:end,6);
% EXP = cellfun(@str2num,content(2:end,7:16));
% EXP0 = min(EXP,[],2); % For reference model
% EXP1 = min(EXP(:,1:2:10),[],2); % For control model
% EXP2 = min(EXP(:,2:2:10),[],2); % For stress model
% EXPlvl0 = cell(length(EXP0),1);
% EXPlvl1 = EXPlvl0;
% EXPlvl2 = EXPlvl0;
% EXPlvl0(EXP0<1) = {'None'};
% EXPlvl1(EXP1<1) = {'None'};
% EXPlvl2(EXP2<1) = {'None'};
% EXPlvl0(EXP0<10&EXP0>1) = {'Low'};
% EXPlvl1(EXP1<10&EXP1>1) = {'Low'};
% EXPlvl2(EXP2<10&EXP2>1) = {'Low'};
% EXPlvl0(EXP0<50&EXP0>10) = {'Medium'};
% EXPlvl1(EXP1<50&EXP1>10) = {'Medium'};
% EXPlvl2(EXP2<50&EXP2>10) = {'Medium'};
% EXPlvl0(EXP0>50) = {'High'};
% EXPlvl1(EXP1>50) = {'High'};
% EXPlvl2(EXP2>50) = {'High'};
% 
% 
% fid = fopen('RNAinputFortINIT.txt','w');
% fprintf(fid,'%s,%s,%s,%s,%s,%s\n','Gene','Tissue','Cell type','Level','Expression type','Reliability');
% for i = 1:length(1:length(EXPlvl0))
%     fprintf(fid,'%s,%s,%s,%s,%s,%s\n',Genes{i},'RefModel','Pancreas',EXPlvl0{i},'Staining','Supportive');
% end
% for i = 1:length(1:length(EXPlvl1))
%     fprintf(fid,'%s,%s,%s,%s,%s,%s\n',Genes{i},'RefModel','Pancreas',EXPlvl1{i},'Staining','Supportive');
% end
% for i = 1:length(1:length(EXPlvl2))
%     fprintf(fid,'%s,%s,%s,%s,%s,%s\n',Genes{i},'RefModel','Pancreas',EXPlvl2{i},'Staining','Supportive');
% end
% fclose(fid);

hpaData = parseHPA('RNAinputFortINIT_4.txt');
%tasksLC = parseTaskList_cheng('common_tasks.xlsx');

%exhpaData = parseHPA('normal_tissue.csv');
%hpaData = parseHPA('tINIT_main_script_3.csv');
tasksLC = parseTaskList('common_tasks.xlsx');

% model = importExcelModel('iCancer_Core.xlsx',false);
model = load('iCancer_Core.mat');
model = model.model
% Parameters setting
params.MSK_DPAR_OPTIMIZER_MAX_TIME = 1000; % added for Mosek 7
params.MSK_DPAR_MIO_TOL_REL_GAP    = 0.02;  % added for Mosek 7
tic()
for i=1:numel(hpaData.tissues)
    if ~exist([ strrep(hpaData.tissues{i},'/','-') ' - ' strrep(hpaData.celltypes{i},'/','-') '.mat'],'file')
        [model deletedRxnsInINIT essentialRxnForTasks addedRxnsForTasks]=getINITModel(model, hpaData.tissues{i}, hpaData.celltypes{i}, hpaData, [], [],[], false,false,tasksLC,params,[]);
        save([ strrep(hpaData.tissues{i},'/','-') ' - ' strrep(hpaData.celltypes{i},'/','-') '.mat'],'model','deletedRxnsInINIT','essentialRxnForTasks','addedRxnsForTasks');
    end
end
toc ()




