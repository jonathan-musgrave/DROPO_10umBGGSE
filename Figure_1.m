
if ~exist('detune_ar','var')
   xx = 1:Nrt;
   xlab = 'Round Trips';
else
   xx = detune_ar./alpha;
   xlab = 'Detuning (\alpha)';
end

figure(1);clf;
subplot(2,1,1);
plot(xx,AA1,'Color',CFF,'linewidth',LW)
hold on
plot(xx,AA2,'Color',CSHG,'linewidth',LW)
hold off
xlabel(xlab,'FontName',FontName,'FontSize',FS,'FontWeight','bold')
ylabel('Peak Power (W)','FontName',FontName,'FontSize',FS,'FontWeight','bold')
set(gca,'FontName',FontName,'FontSize',FS,'FontWeight','bold','linewidth',LW)


if xx(1)>0
    set(gca,'XDIR','Reverse')
end

yy1 = sum(PFF,2).*(1-RR1).*dt./rT;
    yy1 = yy1./max(yy1);
yy2 = sum(PSHG,2).*(1-RR1).*dt./rT;
    yy2 = yy2./max(yy2)*0.7;
subplot(2,1,2);

figure(10);clf;
plot(xx,yy1,'Color',CFF,'linewidth',LW)
hold on
plot(xx,yy2,'Color',CSHG,'linewidth',LW)
hold off
xlabel(xlab,'FontName',FontName,'FontSize',FS,'FontWeight','bold')
ylabel('Transmitted Power (W)','FontName',FontName,'FontSize',FS,'FontWeight','bold')
set(gca,'FontName',FontName,'FontSize',FS,'FontWeight','bold','linewidth',LW)

xlabel(xlab,'FontName',FontName,'FontSize',FS,'FontWeight','normal')
ylabel('Transmitted Power (W)','FontName',FontName,'FontSize',FS,'FontWeight','normal')
set(gca,'FontName',FontName,'FontSize',FS,'FontWeight','normal','linewidth',1)
trange = 2;
ind = Nrt;
inst_w = -gradient(unwrap(angle(AAFF(ind,:))),dt)/2/pi/1E12;


if xx(1)>0
    set(gca,'XDIR','Reverse')
end


