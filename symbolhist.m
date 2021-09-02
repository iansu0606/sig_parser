function symbolhist(inputdata,setting,name,Databitgroup,factor)

symbol_rate=setting(1,1);
spb=setting(1,2);
sample_rate = symbol_rate*spb;
PAM_order=setting(1,3);

amp=factor(1,1);
mean=factor(1,2);
inputdata=inputdata*amp+mean;
figure;
hold on
xx=-(PAM_order*amp+mean)*1.5:(PAM_order*amp+mean)*0.0025:(PAM_order*amp+mean)*1.5;
for i=1:PAM_order
a(i,:)=hist((cell2mat(Databitgroup{i,1})*amp+mean),xx);
end

x=1:length(inputdata);
plot(x/sample_rate*1e12*spb,inputdata.','ob')
cc=[1 0 0;1 0.6 0;1 0.95 0;0 1 0;1 0 1];
for i=1:PAM_order
    plot((a(i,:)/max(max(a))*max(x)*0.1+max(x)*1.03)/sample_rate*1e12*spb,xx,'Linewidth',1,'color',cc(i,:))
end
axis([0 (max(x)*1.16)/sample_rate*1e12*spb -1.25*PAM_order*amp+mean 1.25*PAM_order*amp+mean])

xlabel ('time (ps)','FontSize',16,'FontWeight','bold');
if mean==0
ylabel ('amplitude (a.u.)','FontSize',16,'FontWeight','bold');
else
ylabel ('amplitude (volt)','FontSize',16,'FontWeight','bold');
end
set(gca,'FontSize',16);
set(gca,'FontWeight','bold');
title(name,'FontSize',20,'FontWeight','bold');