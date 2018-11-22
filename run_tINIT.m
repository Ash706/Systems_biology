tasksLC=parseTaskList('common_tasks.xlsx');
params.MSK_DPAR_OPTIMIZER_MAX_TIME = 30000;
params.MSK_DPAR_MIO_TOL_REL_GAP    = 0.02;
hpaData = parseHPA_cheng('ColonCancerGenes.csv');
for i=1:numel(hpaData.tissues)
        [model deletedRxnsInINIT essentialRxnForTasks addedRxnsForTasks]=getINITModel(model, hpaData.tissues{i}, hpaData.celltypes{i}, hpaData, [], [],[], true,true,tasksLC,params,[]);
        save(['/your/path/' strrep(hpaData.tissues{i},'/','-') ' - ' strrep(hpaData.celltypes{i},'/','-') '.mat'],'model','deletedRxnsInINIT','essentialRxnForTasks','addedRxnsForTasks');
end