frange = 6.*4;
detune_plot = -6.12;
detune_plot =-6;
FS = 8;
fontname = 'arial';
figure(5);clf;
try detune_ar;
    
   dontskip = 1;
catch
    warning('Skipping Figure 5 because there is no detune_ar variable meaning there was no pump sweep')
    dontskip  = 0;
end
if dontskip
    subplot(2,1,1);
    [~,detune_inx] = min(abs(detune_ar./alpha-detune_plot));
    SP = abs(fftshift(ifft(ifftshift(AAFF(detune_inx,:))))).^2/max(abs(fftshift(ifft(ifftshift(AAFF(detune_inx,:))))).^2);
    plot(w/2/pi/1e12,10*log10(SP),'Color',[0.7,0,0],'linewidth',LW)
    hold on
    SP1 = abs(fftshift(ifft(ifftshift(AASHG(detune_inx,:))))).^2/max(abs(fftshift(ifft(ifftshift(AASHG(detune_inx,:))))).^2);
    plot(w/2/pi/1e12,10*log10(SP1),'Color',[0.0,0.7,0],'linewidth',LW)
    hold off
    ylim([max(SP)]+[-140,5])
    set(gca,'xlim',[-frange frange])
    YL = get(gca,'ylim')
    yticks(YL(1):20:YL(2)+20)
    yticklabels('')
    xticks([-4,4])
    xticklabels({'-4','4'})
    xlabel('frequency (THz)')
    ylabel('20 db/div')
    set(gca,'fontname','arial')
    
%     ylabel('power (dBm)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
%     xlabel('frequency (THz)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
%     set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
%     h = legend('signal','pump','location','south');legend boxoff
%     set(h,'Fontsize',FS);
    
    frange = -min(t)./1e-12;
    subplot(2,1,2);
    [~,detune_inx] = min(abs(detune_ar./alpha-detune_plot));
    SP = abs((AAFF(detune_inx,:))).^2; 
    SP = SP./max(SP(:));
    plot(t/1e-12,(SP),'r','linewidth',LW)
    hold on
    SP1 = abs((AASHG(detune_inx,:))).^2; 
    SP1 = SP1./max(SP1(:));
    plot(t/1e-12,(SP1),'g','linewidth',LW)
    hold off
    set(gca,'xlim',[-frange frange])
    ylabel('power (Linear)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
    xlabel('time (ps)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
    set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
    h = legend('signal','pump','location','south');legend boxoff
    set(h,'Fontsize',FS);
    xlim([-2,2])
end