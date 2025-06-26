figure(2);clf;
yyaxis left
plot(t*1E12,PFF(ind,:),'Color',CFF,'linewidth',LW);
hold on
plot(t*1E12,PSHG(ind,:),'-','Color',CSHG,'linewidth',LW);
hold off
xlabel('time (ps)','FontName',FontName,'FontSize',FS,'FontWeight','bold')
xlabel('time (ps)','FontName',FontName,'FontSize',FS,'FontWeight','normal')

ylabel('peak power (W)','FontName',FontName,'FontSize',FS,'FontWeight','bold')
ylabel('peak power (W)','FontName',FontName,'FontSize',FS,'FontWeight','normal')
set(gca,'FontName',FontName,'FontSize',FS,'FontWeight','bold','linewidth',LW)
set(gca,'FontName',FontName,'FontSize',FS,'FontWeight','normal','linewidth',1)
[~,inx] = max(PSHG(ind,:));
t0 = t(inx)*1e12;
set(gca,'xlim',[-trange trange]+t0)
%set(gca,'ylim',[0 0.6])
set(gca,'Ycolor','k')
%set(gca, 'YScale', 'log')
try fwhm(t*1E12,PSHG(ind,:))
FWHM_SHG = fwhm(t*1E12,PSHG(ind,:));
FWHM_FF  = fwhm(t*1e12,PFF(ind,:));
disp(strcat('FF FWHM = ',num2str(FWHM_FF),'ps'))
disp(strcat('SHG FWHM = ',num2str(FWHM_SHG),'ps'))
catch
   warning('Soliton Not formed skipping FWHM calculation') 
end
yyaxis right
plot(t*1E12,inst_w,'b--','linewidth',LW);
ylabel('chirp (THz)','FontName',FontName,'FontSize',FS,'FontWeight','bold')
set(gca,'FontName',FontName,'FontSize',FS,'FontWeight','bold','linewidth',LW)
set(gca,'FontName',FontName,'FontSize',FS,'FontWeight','normal','linewidth',1)
set(gca,'Ycolor','b')
%set(gca,'ylim',[-0.8 0.8])
h = legend('FF','SHG','chirp','location','northeast');legend boxoff 
set(h,'Fontsize',FS);
