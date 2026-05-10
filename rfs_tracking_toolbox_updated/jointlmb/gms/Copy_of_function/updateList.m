function [ currList ] = updateList( X,Y,histList )
for id = 1:size(histList,1)
    if(abs(histList(id,1)-X)<=20)&&(abs(histList(id,2)-Y)<=2)
        histList(id,1)=X;
        histList(id,2)=Y;
        histList(id,4) = histList(id,4)+1;
        histList(id,3)=0;
        currList = sortrows(histList,4);
        return
    end
end
persistent idx;
if isempty(idx)
    idx = 1;
end
currList = [histList;[X,Y,0,1,idx]];
idx = idx+1;
end

