function eyeplot(inputeye,setting,name,factor)

symbol_rate=setting(1,1);
spb=setting(1,2);
sample_rate = symbol_rate*spb;
range=setting(1,3);
upsamplenumber=setting(1,4);
eyenumber=setting(1,5);


amp=factor(1,1);
mean=factor(1,2);
inputeye=inputeye*amp+mean;
inputeyerepeat=inputeye(:,1);
inputeyerepeat=[inputeyerepeat(2:end);inputeye(end,end)];
inputeye=[inputeye inputeyerepeat];

figure;
x=0:spb*eyenumber*upsamplenumber;
grid on
plot(x/upsamplenumber/sample_rate*1e12,inputeye.','b')
axis([0 spb*eyenumber/sample_rate*1e12 -1.25*range*amp+mean 1.25*range*amp+mean])
xlabel ('time (ps)','FontSize',16,'FontWeight','bold');
if mean==0
ylabel ('amplitude (a.u.)','FontSize',16,'FontWeight','bold');
else
ylabel ('amplitude (volt)','FontSize',16,'FontWeight','bold');
end
set(gca,'FontSize',16);
set(gca,'FontWeight','bold');
title(name,'FontSize',20,'FontWeight','bold');