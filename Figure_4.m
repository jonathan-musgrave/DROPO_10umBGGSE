figure(4);clf;
ax = 2;
sIIFF = abs(fftshift(ifft(ifftshift(AAFF,ax),[],ax),ax)).^2;
sIIFF = sIIFF./max(sIIFF(:));
IIFF = abs(AAFF).^2;
IIFF = IIFF./max(IIFF(:));

sIISHG = abs(fftshift(ifft(ifftshift(AASHG,ax),[],ax),ax)).^2;
sIISHG = sIISHG./max(sIISHG(:));
IISHG = abs(AASHG).^2;
IISHG = IISHG./max(IISHG(:));

if ~exist('detune_ar','var')
   yy = 1:Nrt;
   ylab = 'Round Trips';
else
   yy = detune_ar./alpha;
   ylab = 'Detuning (\alpha)';
end

SPi = {IIFF,10*log10(sIIFF),IISHG,10*log10(sIISHG)};
x   = {t./1e-12,w,t./1e-12,w};
y   = {yy,yy,yy,yy};
xl = {'Fast Time (ps)', '\omega (THz)','Fast Time (ps)', '\omega (THz)'};
yl  = {ylab,ylab,ylab,ylab};
tit = {'Signal time (1571 nm)','Frequency (1571 nm)','Pump time (785.5 nm)', 'Frequency (785.5 nm)'};
Clim = [-60,0];
for i = 1:4
    % figure(40+i)
    subplot(2,2,i)
    imagesc(y{i},x{i},(SPi{i})'); colorbar;
    ylabel(xl{i},'FontName',FontName,'FontSize',FS,'FontWeight','normal');
    xlabel(yl{i},'FontName',FontName,'FontSize',FS,'FontWeight','normal');
    title(tit{i},'FontName',FontName,'FontSize',FS,'FontWeight','normal')
    xlim([-7,7]);
    if y{i}(1)>0
    else
        set(gca,'YDIR','Normal')
    end
        set(gca,'XDIR','Reverse')
    if i == 2 || i ==4
        set(gca,'CLIM',Clim);
    else
        set(gca,'CLIM',[0,1]);
    end
end
