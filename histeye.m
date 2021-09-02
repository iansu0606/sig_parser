function histeye(inputeye,setting,name,interpfactor,factor,viewx,viewy)

symbol_rate=setting(1,1);
spb=setting(1,2);
sample_rate = symbol_rate*spb;
range=setting(1,3);
upsamplenumber=setting(1,4);
eyenumber=setting(1,5);
interpx=interpfactor(1,1);
interpy=interpfactor(1,2);
histscale=interpfactor(1,3);

amp=factor(1,1);
mean=factor(1,2);
inputeye=inputeye*amp+mean;
x=0:spb*eyenumber*upsamplenumber-1;
xx=-2*range*amp+mean:(range*amp+mean)/histscale:2*range*amp+mean;
b=hist(inputeye.',xx);
brepeat=b(:,1);
brepeat=[brepeat(2:end);b(end,end)];
b=[b brepeat];
y3=-1.25*range*amp+mean:(range*amp+mean)/histscale/interpy:1.25*range*amp+mean;
x3=0:1/interpx:spb*eyenumber*upsamplenumber;
xp=0:spb*eyenumber*upsamplenumber;
[x2,y2]=meshgrid(xp,xx);
[x4,y4]=meshgrid(x3,y3);
d=interp2(x2,y2,b,x4,y4,'cubic');

figure;
surf(x3/upsamplenumber/sample_rate*1e12,y3,d,'EdgeColor','none','LineStyle','none')
colormap(hot)
axis([0 (spb*eyenumber)/sample_rate*1e12 -1.25*range*amp+mean 1.25*range*amp+mean -10 max(max(b))])
xlabel ('time (ps)','FontSize',16,'FontWeight','bold');
if mean==0
ylabel ('amplitude (a.u.)','FontSize',16,'FontWeight','bold');
else
ylabel ('amplitude (volt)','FontSize',16,'FontWeight','bold');
end
set(gca,'FontSize',16);
set(gca,'FontWeight','bold');
title(name,'FontSize',20,'FontWeight','bold');
view(viewx,viewy)