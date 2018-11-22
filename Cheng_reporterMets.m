function repMets=Cheng_reporterMets(model,DEresult,printResults,outputFile)
% repMets = Cheng_reporterMets(model,DEresult,printResults,outputFile)
genes = DEresult.genes;
genePValues = DEresult.genePValues;
geneFoldChanges = DEresult.geneFoldChanges;


if nargin<3
    printResults=false;
end
if nargin<4
	outputFile=[];
end


%Check some stuff
if numel(genes)~=numel(genePValues)
    dispEM('The number of genes and the number of Z-scores must be the same'); 
end
if ~isfield(model,'rxnGeneMat')
    dispEM('The model structure must have a rxnGeneMat field'); 
end

%Remove the genes which are not in the model
geneFoldChanges(~ismember(genes,model.genes))=[];
genePValues(~ismember(genes,model.genes))=[];
genes(~ismember(genes,model.genes))=[];

%Convert p-values to Z-scores
geneZScores=norminv(genePValues)*-1;

%This is to prevent errors if the p-values are exactly 0 or 1
geneZScores(geneZScores==-inf)=-15;
geneZScores(geneZScores==inf)=15;

%For each metabolite, calculate the aggregate Z-score and keep track of
%the number of neighbouring genes
metZScores=nan(numel(model.mets),1);
metNGenes=nan(numel(model.mets),1);
meanZ=nan(numel(model.mets),1);
stdZ=nan(numel(model.mets),1);
for i=1:numel(model.mets)
    %Get the involved rxns
    I=model.S(i,:);
    
    %Get the involved genes
    [crap J]=find(model.rxnGeneMat(I~=0,:));
    
    %Find the genes in the gene list
    K=find(ismember(genes,model.genes(J)));
    
    %Calculate the aggregated Z-score for the metabolite
    if any(K)
        metZScores(i)=sum(geneZScores(K))/sqrt(numel(K));
        meanZ(i)=mean(geneZScores(K));
        stdZ(i)=std(geneZScores(K));
        metNGenes(i)=numel(K);
    end
end

%Remove the metabolites which have no Z-scores
mets=model.mets(~isnan(metZScores));
metNames=model.metNames(~isnan(metZScores));
metZScores(isnan(metZScores))=[];
metNGenes(isnan(metNGenes))=[];
meanZ(isnan(meanZ))=[];
stdZ(isnan(stdZ))=[];

%Then correct for background by calculating the mean Z-score for random
%sets of the same size as the ones that were found for the metabolites
sizes=unique(metNGenes);
nGenes=numel(genes);

for i=1:numel(sizes)
    %Sample 100000 sets for each size. Sample with replacement, which may or
    %may not be the best choice.
    I=ceil(rand(100000,sizes(i))*nGenes);
    J=geneZScores(I);
    K=sum(J,2)/sqrt(sizes(i));

    %Correct all elements of this size
    mK=mean(K);
    stdK=std(K);
    metZScores(metNGenes==sizes(i))=(metZScores(metNGenes==sizes(i))-mK)/stdK;
end

%Calculate the P-values
metPValues=1-normcdf(metZScores);

%Sort the results
[metZScores I]=sort(metZScores,'descend');
mets=mets(I);
metNames=metNames(I);
metPValues=metPValues(I);
metNGenes=metNGenes(I);
meanZ=meanZ(I);
stdZ=stdZ(I);

%Prepare the output
repMets.mets=mets;
repMets.metNames=metNames;
repMets.metZScores=metZScores;
repMets.metPValues=metPValues;
repMets.metNGenes=metNGenes;
repMets.meanZ=meanZ;
repMets.stdZ=stdZ;

%Call this function recursively if geneFoldChanges has been specified
repMets(1).test='all';
if any(geneFoldChanges)
    repMets=[repMets;reporterMets(model,genes(geneFoldChanges>=0),genePValues(geneFoldChanges>=0))];
    repMets=[repMets;reporterMets(model,genes(geneFoldChanges<0),genePValues(geneFoldChanges<0))];
    repMets(2).test='only up';
    repMets(3).test='only down';
end

%This is for printing results to the screen. For printing a full list of
%all scores, specify a output file
if printResults==true
    for i=1:numel(repMets)
        fprintf(['TOP 20 REPORTER METABOLITES\nTEST TYPE: ' repMets(i).test '\nID\tNAME\tP-VALUE\n']);
        for j=1:min(20,numel(repMets(i).mets))
            fprintf([repMets(i).mets{j} '\t' repMets(i).metNames{j} '\t' num2str(repMets(i).metPValues(j)) '\n']);
        end
        fprintf('\n');
    end
end

%Print results to file
if any(outputFile)
    fid=fopen(outputFile,'w');
    for i=1:numel(repMets)
        fprintf(fid,['REPORTER METABOLITES USING TEST TYPE: ' repMets(i).test '\n']);
        fprintf(fid,'ID\tNAME\tZ-SCORE\tP-VALUE\tNUMBER OF NEIGHBOURS\tAVERAGE Z-SCORE\tSTD Z-SCORE\n');
        for j=1:numel(repMets(i).mets)
            fprintf(fid,[repMets(i).mets{j} '\t' repMets(i).metNames{j} '\t' num2str(repMets(i).metZScores(j)) '\t' num2str(repMets(i).metPValues(j)) '\t' num2str(repMets(i).metNGenes(j)) '\t' num2str(repMets(i).meanZ(j)) '\t' num2str(repMets(i).stdZ(j)) '\n']);
        end
    end
    fclose(fid);
end
end
