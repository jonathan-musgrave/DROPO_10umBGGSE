f = figure(7);clf;


C = {[0.1647    0.5098    0.1686],[ 0.7930    0.1641    0.1641]};
plot(t*1E12,PFF(ind,:)./max(PFF(ind,:)),'Color',C{2},'linewidth',LW);
% hold on
% plot(t*1E12,PSHG(ind,:),'-','Color',CSHG,'linewidth',LW);
% hold off
xlabel('time (ps)','FontName',FontName,'FontSize',FS,'FontWeight','normal')
ylabel('Power (arb.)','FontName',FontName,'FontSize',FS,'FontWeight','normal')
set(gca,'FontName',FontName,'FontSize',FS,'FontWeight','normal','linewidth',1)
[~,inx] = max(PFF(ind,:));
t0 = t(inx)*1e12;
set(gca,'xlim',[-trange trange]+t0)
%set(gca,'ylim',[0 0.6])
set(gca,'Ycolor','k')
set(h,'Fontsize',FS);
ylim([0,1.1])
f.Position = [f.Position(1)   f.Position(2)   125    89];
grid on;


f = figure(8);clf;
plot(t*1E12,PSHG(ind,:)./max(PFF(ind,:)),'Color',C{1},'linewidth',LW);
% hold on
% plot(t*1E12,PSHG(ind,:),'-','Color',CSHG,'linewidth',LW);
% hold off
xlabel('time (ps)','FontName',FontName,'FontSize',FS,'FontWeight','normal')
ylabel('Power (arb.)','FontName',FontName,'FontSize',FS,'FontWeight','normal')
set(gca,'FontName',FontName,'FontSize',FS,'FontWeight','normal','linewidth',1)
[~,inx] = max(PSHG(ind,:));
t0 = t(inx)*1e12;
set(gca,'xlim',[-trange trange]+t0)
%set(gca,'ylim',[0 0.6])
set(gca,'Ycolor','k')
set(h,'Fontsize',FS);
ylim([0,1.1])
f.Position = [f.Position(1)   f.Position(2)   125    89];
grid on;
