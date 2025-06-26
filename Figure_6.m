
if ~exist('detune_ar','var')
   xx = 1:Nrt;
   xlab = 'Round Trips';
else
   xx = detune_ar./alpha;
   xlab = 'Detuning (\alpha)';
end
ax = 2;
sPFF = abs(fftshift(fft(fftshift(AAFF,ax),[],ax),ax)./Nw).^2;
sPSHG = abs(fftshift(fft(fftshift(AASHG,ax),[],ax),ax)./Nw).^2;
 
FF0 = zeros(size(sPFF));
SHG0 = zeros(size(sPSHG));
FF_comb = zeros(size(sPFF));
SHG_comb = zeros(size(sPSHG));
if ~exist('nsave','var')
    IND = 1:Nrt;
else
    IND = 1:floor(Nrt./nsave);
end
for ind = IND
    FF0(ind)= (sPFF(ind,Nw/2+1));
    SHG0(ind)= (sPSHG(ind,Nw/2+1));
    
    FF_comb(ind)= sum(sPFF(ind,:))-(sPFF(ind,Nw/2+1));
    SHG_comb(ind)= sum(sPSHG(ind,:))-(sPSHG(ind,Nw/2+1));
end

figure(6);clf;
subplot(3,1,1)
plot(xx,FF0,'r','linewidth',LW)
hold on
plot(xx,SHG0,'g','linewidth',LW)
hold off
xlabel(xlab,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
ylabel('Fundamental power (W)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
if xx(1)>0
    set(gca,'XDIR','Reverse')
end
subplot(3,1,2)
plot(xx,FF_comb,'Color',[0.7,0.1,0.1],'linewidth',LW)
hold on
plot(xx,SHG_comb,'Color',[0.1,0.7,0.1],'linewidth',LW)
hold off
xlabel(xlab,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
ylabel('Comb Power (W)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
if xx(1)>0
    set(gca,'XDIR','Reverse')
end

subplot(3,1,3)
plot(xx,FF_comb+FF0,'Color',[0.5,0.1,0.1],'linewidth',LW)
hold on
plot(xx,SHG_comb+SHG0,'Color',[0.1,0.5,0.1],'linewidth',LW)
hold off
xlabel(xlab,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
ylabel('Total Power (W)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)

if xx(1)>0
    set(gca,'XDIR','Reverse')
end
f = figure(6);
f.Position = [100 0   560   760*3/2];

