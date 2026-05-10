function [] = plotList( histList )
for id = 1:size(histList,1)
    if histList(id,3)==0&&histList(id,4)==1
        plot(histList(id,2),histList(id,1),'wo','linewidth',2)
    elseif histList(id,3)==0&&histList(id,4)>1
        plot(histList(id,2),histList(id,1),'go','linewidth',2)
        text(histList(id,2)+1,histList(id,1),num2str(histList(id,end)),'Color','green')
    elseif histList(id,3)>=1
        plot(histList(id,2),histList(id,1),'blackx','linewidth',2)
    end
end
end