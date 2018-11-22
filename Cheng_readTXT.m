function CellArray = Cheng_readTXT(fileName,delimiter)
% CellArray = Cheng_readTXT(fileName,delimiter)
if nargin < 2|| isempty(delimiter)
    delimiter = '\t';
end
fid = fopen(fileName,'r');
fileContent = textscan(fid,'%s','delimiter','\n');
firstLine = textscan(fileContent{1}{1},'%s','delimiter',delimiter);
CellArray = cell(length(fileContent{1}),length(firstLine{1}));
for i = 1:length(fileContent{1})
     iLine = textscan(fileContent{1}{i},'%s','delimiter',delimiter);
     if length(iLine{1}) == length(CellArray(i,:))
        CellArray(i,:) = iLine{1};
     else
        CellArray(i,1:length(iLine{1})) = iLine{1};
     end
end




end