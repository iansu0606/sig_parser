function [eye2]=eyeanalyze(inputdata,Pam4levelmean,name,setting)
%% setting

sample_rate=setting(1,1);
spb=setting(1,2);
PAM_order=setting(1,3);
crosstime.precent=setting(1,4);
risefall.precent=setting(1,5);
BERthrehold=setting(1,6);
SNR.precent=setting(1,7);
eyeheight.prob=[];

for i=1:PAM_order-1
    eye2.amp(i)=Pam4levelmean(i+1)-Pam4levelmean(i);
    eye2.amp2(i)=(Pam4levelmean(i+1)-Pam4levelmean(i))/2;
    level(i)=(i-PAM_order/2)*2+1;
    refamp(i)=Pam4levelmean(i)+eye2.amp(i)-1;
end

%% Input

a=fix(length(inputdata)/2/spb);
inputdata=inputdata(1:2*spb*a);
inputeye=inputdata(1+spb/2:end-3*spb/2);
inputeyedouble=reshape(inputeye,2*spb,[]);
inputeyesingle=reshape(inputeye,spb,[]);
inputeyesinglemid=reshape(inputdata,spb,[]);

viewrange=PAM_order;
xsingle=1:spb;
xsingleintp=1:0.1:spb;
xdouble=1:spb*2;
xdoubleintp=1:0.1:spb*2;
y=-viewrange:viewrange/100:viewrange;
yintp=-viewrange:viewrange/800:viewrange;

% if PAM_order==2
%     divgroup=1;
% else
    divgroup=PAM_order-2;
% end

doublehist=hist(inputeyedouble.',y);
[xdouble2,ydouble2]=meshgrid(xdouble,y);
[xdouble4,ydouble4]=meshgrid(xdoubleintp,yintp);
doubleeye=interp2(xdouble2,ydouble2,doublehist,xdouble4,ydouble4,'cubic');

% figure;
% surf(xdoubleintp/sample_rate*1e12,yintp,doubleeye,'EdgeColor','none','LineStyle','none')
% colormap(hot)
% axis([0 spb*2/sample_rate*1e12 -viewrange viewrange])
% xlabel ('time (ps)','FontSize',16,'FontWeight','bold');
% ylabel ('normalized amplitude (v)','FontSize',16,'FontWeight','bold');
% set(gca,'FontSize',16);
% set(gca,'FontWeight','bold');
% view(0,90)

singlehist=hist(inputeyesingle.',y);
[xsingle2,ysingle2]=meshgrid(xsingle,y);
[xsingle4,ysingle4]=meshgrid(xsingleintp,yintp);
singaleye=interp2(xsingle2,ysingle2,singlehist,xsingle4,ysingle4,'cubic');

% figure;
% surf(xsingleintp/sample_rate*1e12,yintp,singaleye,'EdgeColor','none','LineStyle','none')
% colormap(hot)
% axis([0 spb/sample_rate*1e12 -viewrange viewrange])
% xlabel ('time (ps)','FontSize',16,'FontWeight','bold');
% ylabel ('normalized amplitude (v)','FontSize',16,'FontWeight','bold');
% set(gca,'FontSize',16);
% set(gca,'FontWeight','bold');
% view(0,90)


singlemid=hist(inputeyesinglemid.',y);
singleeyemid=interp2(xsingle2,ysingle2,singlemid,xsingle4,ysingle4,'cubic');

% figure;
% surf(xsingleintp/sample_rate*1e12,yintp,singleeyemid,'EdgeColor','none','LineStyle','none')
% colormap(hot)
% axis([0 spb/sample_rate*1e12 -viewrange viewrange])
% xlabel ('time (ps)','FontSize',16,'FontWeight','bold');
% ylabel ('normalized amplitude (v)','FontSize',16,'FontWeight','bold');
% set(gca,'FontSize',16);
% set(gca,'FontWeight','bold');
% view(0,90)

%% SNR 10%

SNR.histrange=[floor(ceil(length(xsingle))*(0.5-SNR.precent/100)) ceil(ceil(length(xsingle))*(0.5+SNR.precent/100))];
SNR.hist=inputeyesinglemid(SNR.histrange(1):SNR.histrange(2),:);
for j=1:SNR.histrange(2)-SNR.histrange(1)+1
    [SNR.levelmean SNR.levelstd SNR.Databitgroup]=RxMeanStd(SNR.hist(j,:),PAM_order);
    for i=1:PAM_order-1
        SNR.Qsweep(j,i)=(SNR.levelmean(1,i+1)-SNR.levelmean(1,i))/(SNR.levelstd(1,i+1)+SNR.levelstd(1,i));
    end
end

SNR.Q=10*log10(sum((SNR.Qsweep.^2).'));
eye2.SNR=min(SNR.Q);

%% Cross time

crosstime.range=[refamp-ones(1,PAM_order-1).*crosstime.precent/100;refamp+ones(1,PAM_order-1)*crosstime.precent/100];

for i=1:PAM_order-1
    for j=1:2
    [error crosstime.trackpoint(j,i)]=min(abs(y-ones(1,length(y))*crosstime.range(j,i)));
    end
end
    crosstime.trackhist=[];
for j=1:PAM_order-1
    for i=crosstime.trackpoint(1,1):crosstime.trackpoint(2,1)
        crosstime.trackhist(j,:)=sum(doublehist(crosstime.trackpoint(1,j):crosstime.trackpoint(2,j),:));
    end
end
crosstime.fronthist=crosstime.trackhist(:,1:spb);
crosstime.backhist=crosstime.trackhist(:,1+spb:end);
x=[];
for i=1:PAM_order-1
    x=[x;1:spb];
end
crosstime.front=sum((x.*crosstime.fronthist).')./sum(crosstime.fronthist.');
crosstime.back=sum((x.*crosstime.backhist).')./sum(crosstime.backhist.');
eye2.crosstime=(spb*ones(1,PAM_order-1)+crosstime.back-crosstime.front)/spb;

%% Cross amp

crossamp.pos=[mean(crosstime.front) mean(crosstime.back)+32];

for j=1:2
    [error crossamp.trackpoint(j)]=min(abs(xdoubleintp-ones(1,length(xdoubleintp))*crossamp.pos(1,j)));
end
crossamp.hist=[doubleeye(:,crossamp.trackpoint(1)) doubleeye(:,crossamp.trackpoint(2))];

crossamp.trackpoint2=[];
for j=1:divgroup
    [error crossamp.trackpoint2(j)]=min(abs(yintp-ones(1,length(yintp))*Pam4levelmean(1,1+j)));
end
crossamp.trackpoint2=[0 crossamp.trackpoint2 length(yintp)];
crossamp.group=cell(PAM_order-1,1);
refx=cell(PAM_order-1,1);
yintpb=yintp.';

for i=1:divgroup+1
    for j=1:2
        crossamp.group(i,j)=[{crossamp.hist(1+crossamp.trackpoint2(i):crossamp.trackpoint2(i+1),j)}];
    end
    refx(i,1)=[{yintpb(1+crossamp.trackpoint2(i):crossamp.trackpoint2(i+1),1)}];
end

for i=1:PAM_order-1
    for j=1:2
        crossamp.level(i,j)=sum(crossamp.group{i,j}.*refx{i,1})/sum(crossamp.group{i,j});
        crossamp.amp(i,j)=crossamp.level(i,j)-Pam4levelmean(1,i);
    end
    eye2.crossamp(1,i)=mean(crossamp.amp(i,:))/eye2.amp(1,i);
end

%% Rise Fall

risefall.threhold=[eye2.amp*risefall.precent/100+Pam4levelmean(1:PAM_order-1);eye2.amp*(1-risefall.precent/100)+Pam4levelmean(1:PAM_order-1)];

for i=1:PAM_order-1
    for j=1:2
    [error risefall.trackpoint(j,i)]=min(abs(yintp-ones(1,length(yintp))*risefall.threhold(j,i)));
    end
end

for i=1:PAM_order-1
        risefall.risehist(i,:)=singaleye(risefall.trackpoint(1,i),:);
        risefall.fallhist(i,:)=singaleye(risefall.trackpoint(2,i),:);
end

risefall.frontrisehist=risefall.risehist(:,1:ceil(length(xsingleintp)/2));
risefall.backrisehist=risefall.risehist(:,ceil(length(xsingleintp)/2):end);
risefall.frontfallhist=risefall.fallhist(:,1:ceil(length(xsingleintp)/2));
risefall.backfallhist=risefall.fallhist(:,ceil(length(xsingleintp)/2):end);


for i=1:PAM_order-1
    risefall.frontrisetime(1,i)=sum((risefall.frontrisehist(i,:).*xsingleintp(1:ceil(length(xsingleintp)/2))))/sum(risefall.frontrisehist(i,:));
    risefall.backrisetime(1,i)=sum((risefall.backrisehist(i,:).*xsingleintp(ceil(length(xsingleintp)/2):end)))/sum(risefall.backrisehist(i,:));
    risefall.frontfalltime(1,i)=sum((risefall.frontfallhist(i,:).*xsingleintp(1:ceil(length(xsingleintp)/2))))/sum(risefall.frontfallhist(i,:));
    risefall.backfalltime(1,i)=sum((risefall.backfallhist(i,:).*xsingleintp(ceil(length(xsingleintp)/2):end)))/sum(risefall.backfallhist(i,:));
end

eye2.rise=(risefall.backrisetime-risefall.frontrisetime)/spb;
eye2.fall=(risefall.backfalltime-risefall.frontfalltime)/spb;

%% Eye height

eyeheight.pos=(sum(crosstime.front.*sum(crosstime.fronthist.'))+sum(crosstime.back.*sum(crosstime.backhist.')))/sum(sum(crosstime.fronthist.')+sum(crosstime.backhist.'));

[error eyeheight.trackpoint]=min(abs(xsingleintp-ones(1,length(xsingleintp))*eyeheight.pos));
for i=1:PAM_order-1
    [error eyeheight.seppoint(1,i)]=min(abs(yintp-ones(1,length(yintp))*refamp(1,i)));
end

eyeheight.hist=singleeyemid(:,eyeheight.trackpoint);

eyeheight.grouplow=cell(PAM_order-1,1);
eyeheight.grouphigh=cell(PAM_order-1,1);
eyeheight.lowprob=cell(PAM_order-1,1);
eyeheight.highprob=cell(PAM_order-1,1);

for i=1:divgroup+1
    eyeheight.grouplow(i,1)=[{eyeheight.hist(1+crossamp.trackpoint2(i):eyeheight.seppoint(i),1)}];
    eyeheight.grouphigh(i,1)=[{eyeheight.hist(1+eyeheight.seppoint(i):crossamp.trackpoint2(i+1),1)}];
    eyeheight.grouplowsum(i,1)=sum(abs(cell2mat(eyeheight.grouplow(i,1))));
    eyeheight.grouphighsum(i,1)=sum(abs(cell2mat(eyeheight.grouphigh(i,1))));
    for k=length(eyeheight.grouplow{i,1}):-1:1
        eyeheight.lowprob{i,1}=[{sum(abs(eyeheight.grouplow{i,1}(k:end,1)))/eyeheight.grouplowsum(i,1)} eyeheight.lowprob{i,1}];
    end
    for k=1:length(eyeheight.grouphigh{i,1})
        eyeheight.highprob{i,1}=[eyeheight.highprob{i,1} {sum(abs(eyeheight.grouphigh{i,1}(1:k,1)))/eyeheight.grouphighsum(i,1)}];
    end
end

for i=1:divgroup+1
    eyeheight.lowpos(1,i)=max(find(cell2mat(eyeheight.lowprob{i,1})>BERthrehold));
    eyeheight.highpos(1,i)=min(find(cell2mat(eyeheight.highprob{i,1})>BERthrehold));
    eyeheight.prob=[eyeheight.prob cell2mat(eyeheight.lowprob{i,1}) cell2mat(eyeheight.highprob{i,1})];
end

eyeheight.BERprob=ones(1,length(yintp))*BERthrehold;
figure;
semilogy(yintp,eyeheight.prob,yintp,eyeheight.BERprob,'-.m','Linewidth',3);
xlabel ('normalized amplitude (v)','FontSize',16,'FontWeight','bold');
ylabel ('Probability','FontSize',16,'FontWeight','bold');
axis([-PAM_order PAM_order 1e-8 1])
set(gca,'FontSize',16);
set(gca,'FontWeight','bold');
eyehight.name=['Bathtub of eyehight ',name];
title(eyehight.name,'FontSize',20,'FontWeight','bold');

eyeheight.lowposre=eyeheight.lowpos+crossamp.trackpoint2(1:1:PAM_order-1);
eyeheight.highposre=eyeheight.highpos+eyeheight.seppoint;

eye2.eyeheight=((eyeheight.highposre-eyeheight.lowposre)*PAM_order/800)./eye2.amp;
%     eye2.height(i)=Pam4levelmean(i+1)-3*Pam4levelstd(i+1)-Pam4levelmean(i)-3*Pam4levelstd(i);

%% Eye width

eyewidth.pos=eyeheight.seppoint;
for i=1:divgroup+1
    eyewidth.hist(i,:)=singleeyemid(eyewidth.pos(1,i),:);
    eyewidth.histleft(i,:)=eyewidth.hist(i,1:ceil(length(eyewidth.hist)/2));
    eyewidth.histright(i,:)=eyewidth.hist(i,ceil(length(eyewidth.hist)/2):end);
end

eyewidth.histleftsum=sum(abs(eyewidth.histleft.'));
eyewidth.histrightsum=sum(abs(eyewidth.histright.'));
    
eyewidth.leftprob=[];
eyewidth.rightprob=[];

for k=length(eyewidth.histleft):-1:1
    eyewidth.leftprob=[(sum(abs(eyewidth.histleft(:,k:end).'))./eyewidth.histleftsum).' eyewidth.leftprob];
end
for k=1:length(eyewidth.histright)
    eyewidth.rightprob=[eyewidth.rightprob (sum(abs(eyewidth.histright(:,1:k).'))./eyewidth.histrightsum).'];
end
    
eyewidth.prob=[];
for i=1:divgroup+1
    eyewidth.leftpos(1,i)=max(find(eyewidth.leftprob(i,:)>BERthrehold));
    eyewidth.rightpos(1,i)=min(find(eyewidth.rightprob(i,:)>BERthrehold));
    eyewidth.prob(i,:)=[eyewidth.leftprob(i,:) eyewidth.rightprob(i,2:end)];
end

eyewidth.BERprob=ones(1,length(xsingleintp))*BERthrehold;

figure;
semilogy(xsingleintp/sample_rate*1e12,eyewidth.prob,xsingleintp/sample_rate*1e12,eyewidth.BERprob,'-.m','Linewidth',3);
xlabel ('time (ps)','FontSize',16,'FontWeight','bold');
ylabel ('Probability','FontSize',16,'FontWeight','bold');
axis([0 max(xsingleintp)/sample_rate*1e12 1e-8 1])
set(gca,'FontSize',16);
set(gca,'FontWeight','bold');
eyewidth.name=['Bathtub of eyewidth ',name];
title(eyewidth.name,'FontSize',20,'FontWeight','bold');

eyewidth.leftposre=eyewidth.leftpos;
eyewidth.rightposre=eyewidth.rightpos+ceil(length(eyewidth.hist)/2);
eye2.eyewidth=(eyewidth.rightposre-eyewidth.leftposre)./length(xsingleintp);
% clear trackpoint;
% for i=1:PAM_order-1
%     [error trackpoint(i)]=min(abs(y5-ones(1,length(y5))*eye2.crossamp(i)));
% end
% trackhist=[];
% trackhist=b(trackpoint,:);
% xup=1:1:spb/2;
% xdown=spb/2+1:1:spb;
% yup=trackhist(:,1:length(xup));
% ydown=trackhist(:,1+length(xup):end);
% upvalue=cell(PAM_order-1,1);
% downvalue=cell(PAM_order-1,1);
% 
% for k=1:PAM_order-1
%     upvaluetemp=[];
%     downvaluetemp=[];
%     for i=1:length(xup)
%         upvaluetemp=[upvaluetemp ones(1,yup(k,i))*xup(i)];
%         downvaluetemp=[downvaluetemp ones(1,ydown(k,i))*xdown(i)];
%     end
%     upvalue{k,1}=[{upvaluetemp}];
%     downvalue{k,1}=[{downvaluetemp}];
% end
% for x=1:PAM_order-1
%     upmean(1,x)=mean(cell2mat(upvalue{x,1}));
%     upstd(1,x)=std(cell2mat(upvalue{x,1}));
%     downmean(1,x)=mean(cell2mat(downvalue{x,1}));
%     downstd(1,x)=std(cell2mat(downvalue{x,1}));
% end
% eye2.width=(-upmean-3*upstd+downmean-3*downstd)/sample_rate*1e12;