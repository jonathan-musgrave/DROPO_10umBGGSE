frange = 6;


figure(3);clf;
SP = abs(fftshift(ifft(ifftshift(AAFF(ind,:))))).^2;%/max(abs(fftshift(ifft(ifftshift(AAFF(ind,2:end))))).^2);
SP1 = abs(fftshift(ifft(ifftshift(AASHG(ind,:))))).^2;%/max(abs(fftshift(ifft(ifftshift(AASHG(ind,[1:Nw/2,Nw/2+2:end]))))).^2);
C = max([SP(1:Nw/2),SP1(1:Nw/2)]);
SP = SP./C;
SP1 = SP1./C;
plot(w/2/pi/1e12,10*log10(SP),'Color',CFF,'linewidth',LW)
hold on
%SP1 = SP1./max(SP1(1:Nw/2));
plot(w/2/pi/1e12,10*log10(SP1),'Color',CSHG,'linewidth',LW)
hold off
set(gca,'xlim',[-frange frange])
ylabel('power (dBm)','FontName',FontName,'FontSize',FS,'FontWeight','bold')
xlabel('frequency (THz)','FontName',FontName,'FontSize',FS,'FontWeight','bold')
set(gca,'FontName',FontName,'FontSize',FS,'FontWeight','bold','linewidth',LW)
ylabel('power (dBm)','FontName',FontName,'FontSize',FS,'FontWeight','normal')
xlabel('frequency (THz)','FontName',FontName,'FontSize',FS,'FontWeight','normal')
set(gca,'FontName',FontName,'FontSize',FS,'FontWeight','normal','linewidth',1)

% h = legend('FF','SHG','location','south');legend boxoff
set(h,'Fontsize',FS);
