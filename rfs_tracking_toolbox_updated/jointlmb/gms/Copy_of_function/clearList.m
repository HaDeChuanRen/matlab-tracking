function [ currList ] = clearList( histList )
if ~isempty(histList)
    histList(:,3) = histList(:,3)+1;
    for id=size(histList,1):-1:1
        if histList(id,3)>3
            histList(id,:)=[];
        end
    end
end
currList = histList;
end