function showfigures(label,Td,preseam,costmap,SG,rr,step,stepname)
%SHOWFIGURES Summary of this function goes here
%   Detailed explanation goes here
global nbvar

figure(1); clf;
imagesc(Td); axis xy equal tight;
title('define terminal');
if size(preseam,1) > 0
    hold on;
    scatter(preseam(:,2),preseam(:,1));
end
print('-dpng',[num2str(rr),stepname,'_terminal_',num2str(step)]);

figure(2); clf;
imagesc(label); axis xy equal tight;
title('cut label');
if size(preseam,1) > 0
    hold on;
    scatter(preseam(:,2),preseam(:,1));
end
print('-dpng',[num2str(rr),stepname,'_label_',num2str(step)]);

figure(3); clf;
imagesc(costmap); axis xy equal tight;
title('cost map');
print('-dpng',[num2str(rr),stepname,'_costmap_',num2str(step)]);

figure(4);
for k = 1:nbvar
    clf;
    imagesc(SG(:,:,k)); axis xy equal tight;
    title(['SG',num2str(k)]);
    print('-dpng',[num2str(rr),stepname,'_SGvar',num2str(k),'_',num2str(step)]);
end

end

