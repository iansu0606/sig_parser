clear all
clc
close all


%% Path

% for n=12:2:20
% 
%         name=['C:\Users\Jim\Desktop\EQ\2km\-',num2str(n),'\Rx.txt'];
%         Rxpower=load(name);
%         name2=['C:\Users\Jim\Desktop\EQ\2km\-',num2str(n),'\Tx.txt'];
%         Txpower=load(name2);



% for power=0:10:80
%     if power==0
%         name=['C:\Users\Jim\Desktop\EQ\data\bitrate50G10km\RxPower100.txt'];
%         Rxpower=load(name);
%         name2=['C:\Users\Jim\Desktop\EQ\data\bitrate50G10km\TxPower100.txt'];
%         Txpower=load(name2);
%     elseif power==5
%         name=['C:\Users\Jim\Desktop\EQ\data\bitrate50G10km\RxPower105.txt'];
%         Rxpower=load(name);
%         name2=['C:\Users\Jim\Desktop\EQ\data\bitrate50G10km\TxPower105.txt'];
%         Txpower=load(name2);
%     else
%         name=['C:\Users\Jim\Desktop\EQ\data\bitrate50G10km\RxPower1',num2str(power),'.txt'];
%         Rxpower=load(name);
%         name2=['C:\Users\Jim\Desktop\EQ\data\bitrate50G10km\TxPower1',num2str(power),'.txt'];
%         Txpower=load(name2);
%     end
    
% Rxpower=load('C:\Users\Jim\Desktop\EQ\PAM450Gbps\RxPower100.txt');
% Txpower=load('C:\Users\Jim\Desktop\EQ\PAM450Gbps\TxPower100.txt');


%% Setting
PAM_order = 4;
PRBS_length = 11;
spb = 32;                                                                       % sample per symbol
symbol_rate = 50e9;                                                             % bit_rate/2 = baud_rate
sample_rate = symbol_rate*spb;                                                  % bit_rate*sample_count_per_bit

reflength=100;

tapsweep=0;                                                                     % EQ tap
FFEtap=40;
DFEtap=FFEtap-1;

crosstimeprecent=10;
risefallprecent=10;
BERthrehold=1e-4;
SNRprecent=10;

eqenable=1;
filter_type=0;
filterbandwidth=symbol_rate/2;

if PAM_order==2
    Txpower=load ('C:\Users\Jim\Desktop\EQ\2km\-12\Tx.txt');
    Rxpower=load ('C:\Users\Jim\Desktop\EQ\2km\-12\Rx.txt');
    Pam4Level1 = 0;
    Pam4Level2 = 0.09;
    Tx.Pam4level=[Pam4Level1 Pam4Level2];
    viewrange=0.02;
elseif PAM_order==4
%     Rxpower=load('C:\Users\Jim\Desktop\EQ\data\bitrate50G10km\RxPower100.txt');
%     Txpower=load('C:\Users\Jim\Desktop\EQ\data\bitrate50G10km\TxPower100.txt');
%     Txpower=load('C:\Users\Jim\Desktop\EQ\data\50G_2km\Tx_10.txt');
%     Rxpower=load('C:\Users\Jim\Desktop\EQ\data\50G_2km\Rx_10.txt');

    Txpower=load('D:\A-Ted\20160824\56GEML\tx_5a.txt');           %Tx signal
    Rxpower=load('D:\A-Ted\20160824\56GEML\rx_5a.txt');           %Rx signal
    Pam4Level1 = 0.03;
    Pam4Level2 = 0.06;
    Pam4Level3 = 0.09;
    Pam4Level4 = 0.12;
    viewrange=0.02;
    Tx.Pam4level=[Pam4Level1 Pam4Level2 Pam4Level3 Pam4Level4];
end

Rx.Data = Rxpower(:,2);
Tx.Ref = Txpower(:,2);

analyzesetting=[sample_rate,spb,PAM_order,crosstimeprecent,risefallprecent,BERthrehold,SNRprecent];

%% filter

    beta=filterbandwidth/symbol_rate*2;
    if filter_type ~= 0
        matching_filter = rcosfir(0.3,10,spb,1,'sqrt');
        data_scope_up_temp2 = conv(Rx.Data,matching_filter);
        Rx.Datafilter = data_scope_up_temp2((length(matching_filter)-1)/2+1:end-(length(matching_filter)-1)/2);
    else
        Rx.Datafilter = Rx.Data;
    end

% fvtool(matching_filter, 'Analysis', 'impulse')
%% Correlation

tic;

FindDelay = Tx.Ref(1:1000);
for start = 1:500;
    y = corrcoef(Rx.Datafilter(start:start-1+1000),FindDelay);
    SyncChoice(start,1) = y(2,1);
end
[PeakMax ShiftLength] = max(SyncChoice);

Tx.Refcorrelate = Tx.Ref(1:(2^PRBS_length-1-floor(ShiftLength/spb))*spb,1);
Rx.Datacorrelate = Rx.Datafilter(ShiftLength:(2^PRBS_length-1-floor(ShiftLength/spb))*spb+ShiftLength-1,1);

%% Signal normalize

Tx.Refcorrelateshift = Tx.Refcorrelate-mean(Tx.Refcorrelate);
Tx.normalized = Tx.Refcorrelateshift/mean(abs(Tx.Refcorrelateshift))*log2(PAM_order);

Rx.Datacorrelateshift = Rx.Datacorrelate-mean(Rx.Datacorrelate);
Rx.normalized = Rx.Datacorrelateshift/mean(abs(Rx.Datacorrelateshift))*log2(PAM_order);

% signalfft(Tx.normalized,sample_rate)
% signalfft(Rx.normalized,sample_rate)

%% Eye Diagram of Training and Data

Tx.refeye = reshape(Tx.Refcorrelateshift(1:2*spb*floor(length(Tx.Refcorrelate)/2/spb)),2*spb,[]);
Rx.inputeye = reshape(Rx.normalized(1:2*spb*floor(length(Tx.Refcorrelate)/2/spb)),2*spb,[]);

eyeplot(Rx.inputeye.',sample_rate,spb,'Eyediagram of Data before EQ',PAM_order)
histeye(Rx.inputeye,sample_rate,spb,'Eyediagram of Data before EQ',PAM_order)

%% ADC of Training

Tx.Refreshape = reshape(Tx.Refcorrelate,spb,[]);
Tx.Refbit = Tx.Refreshape(spb/2,:);                                         % symbol select
[Tx.Symbol]=BittoLevel(Tx.Refbit,Tx.Pam4level,PAM_order);
Tx.Train = pammod(Tx.Symbol,PAM_order,0);                          % PAM-4 Gemeration

%% Find min BER pos

Qsweep=zeros(spb,PAM_order-1);
BER.estsweep=zeros(spb,PAM_order-1);
BER.esttotalsweep=zeros(1,spb);
Rx.Datareshape = reshape(Rx.normalized,spb,[]);

for eyeopenpos=1:spb

    Rx.Databit = Rx.Datareshape(eyeopenpos,:);
    [Rx.Pam4levelmean Rx.Pam4levelstd Rx.Databitgroup]=RxMeanStd(Rx.Databit,PAM_order);

    for i=1:PAM_order-1
        Qsweep(eyeopenpos,i)=(Rx.Pam4levelmean(1,i+1)-Rx.Pam4levelmean(1,i))/(Rx.Pam4levelstd(1,i+1)+Rx.Pam4levelstd(1,i));
        BER.estsweep(eyeopenpos,i)=0.5*erfc(Qsweep(eyeopenpos,i)/sqrt(2));
    end
    
    BER.esttotalsweep(1,eyeopenpos)=sum(BER.estsweep(eyeopenpos,:));
end

[eyeopenmin eyeopenpos]=min(BER.esttotalsweep);

%% ADC of Data and level shift

Rx.Databit = Rx.Datareshape(eyeopenpos,:);
[Rx.Pam4levelmean Rx.Pam4levelstd Rx.Databitgroup]=RxMeanStd(Rx.Databit,PAM_order);
Rx.Databitnormalize=Rx.Databit;

symbolhist(Rx.Databitnormalize,sample_rate,PAM_order,spb,'Symbol of Data before EQ',Rx.Databitgroup)
% eyeanaorigin=eyeanalyze(Rx.normalized,Rx.Pam4levelmean,'before EQ',analyzesetting);
% eye(1,:)=eyecomanalyze(Rx.normalized,sample_rate,spb,PAM_order);

%% BER Estimate w/o EQ

Q=zeros(1,PAM_order-1);
BER.est=zeros(1,PAM_order-1);

for i=1:PAM_order-1
    Q(1,i)=(Rx.Pam4levelmean(1,i+1)-Rx.Pam4levelmean(1,i))/(Rx.Pam4levelstd(1,i+1)+Rx.Pam4levelstd(1,i));
    BER.est(1,i)=0.5*erfc(Q(1,i)/sqrt(2));
end
BER.esttotal=sum(BER.est);

%% BER Count w/o EQ

    [Rx.Symbol]=BittoLevel(Rx.Databit,Rx.Pam4levelmean,PAM_order);
    BER.count=sum(sum(dec2bin(Rx.Symbol) ~= dec2bin(Tx.Symbol)))/length(dec2bin(Rx.Symbol))/log2(PAM_order);
    SNR.woEQ=10*log10(sum(Q.^2));
    SNR.woEQD=10*log10((Q.^2));
    
    
    if eqenable==1;
%% EQ
fftap=5;
dfetap=5;
Ref=Tx.normalized(1:spb*reflength).';
% Ref=[];
if PAM_order==2
    lmssize=0.001;
    lmssize2=lmssize/1e2;
elseif PAM_order==4
    lmssize=0.0005;
    lmssize2=lmssize/1e5;
end

[Rx.DataFFE data_rx_eq errCount]=ffeblock(Rx.normalized,Ref,PAM_order,lmssize,fftap,0);
% Rx.DataFFE=Rx.DataFFE/mean(abs(Rx.DataFFE))*log2(PAM_order);
Rx.DataDFE=dfeblock(Rx.normalized,Ref,PAM_order,lmssize,lmssize2,fftap,dfetap);
% Rx.DataDFE=Rx.DataDFE/mean(abs(Rx.DataDFE))*log2(PAM_order);

figure
hold on
plot(Tx.normalized,'o')
plot(Rx.DataFFE,'or')
plot(Rx.normalized,'og')
% plot(Rx.DataDFE,'oblack')

%% recorrelate

for eqtype=1:2
    if eqtype==1
        Rx.DataEQ=Rx.DataFFE;
    elseif eqtype==2
        Rx.DataEQ=Rx.DataDFE;
    end

FindDelay = Tx.Ref(1:1000);
for start = 1:500;
    y = corrcoef(Rx.DataEQ(start:start-1+1000),FindDelay);
    SyncChoice(start,1) = y(2,1);
end
[PeakMax ShiftLength2] = max(SyncChoice);

Rx.DataEQcorrelate=Rx.DataEQ(ShiftLength2:(2^PRBS_length-2-floor(ShiftLength2/spb))*spb+ShiftLength2-1,1);
Rx.DataEQcorrelate2=Rx.DataEQcorrelate(1+reflength*spb:end);
Tx.EQ=Tx.normalized(1:(2^PRBS_length-2-floor(ShiftLength2/spb))*spb,1);

Tx.TrainEQ=Tx.Train(1:(2^PRBS_length-2-floor(ShiftLength2/spb)));
Tx.SymbolEQ=Tx.Symbol(1:(2^PRBS_length-2-floor(ShiftLength2/spb)));
    
%% Eye Diagram of Training and Data

Qsweep=zeros(spb,PAM_order-1);
BER.estsweep=zeros(spb,PAM_order-1);
BER.esttotalsweep=zeros(1,spb);
Rx.DataEQreshape = reshape(Rx.DataEQcorrelate,spb,[]);    

for eyeopenpos=1:spb

    Rx.DatabitEQ = Rx.DataEQreshape(eyeopenpos,:);
    [Rx.Pam4levelmeanEQ Rx.Pam4levelstdEQ Rx.DatabitEQgroup]=RxMeanStd(Rx.DatabitEQ,PAM_order);

    for i=1:PAM_order-1
        Qsweep(eyeopenpos,i)=(Rx.Pam4levelmeanEQ(1,i+1)-Rx.Pam4levelmeanEQ(1,i))/(Rx.Pam4levelstdEQ(1,i+1)+Rx.Pam4levelstdEQ(1,i));
        BER.estsweep(eyeopenpos,i)=0.5*erfc(Qsweep(eyeopenpos,i)/sqrt(2));
    end
    
    BER.esttotalsweep(1,eyeopenpos)=sum(BER.estsweep(eyeopenpos,:));
end

[eyeopenmin eyeopenposEQ]=min(BER.esttotalsweep);

%% Eye Diagram

    Rx.outputeye = reshape(Rx.DataEQcorrelate(length(Ref)+1:2*spb*floor(length(Tx.EQ)/2/spb)),2*spb,[]);
    Rx.DataEQreshape = reshape(Rx.DataEQcorrelate,spb,[]);    
    Rx.DatabitEQ=Rx.DataEQreshape(eyeopenposEQ,:);
    Rx.DatabitEQ=Rx.DatabitEQ(1,reflength+1:end);
    [Rx.Pam4levelmeanEQ Rx.Pam4levelstdEQ Rx.DatabitEQgroup]=RxMeanStd(Rx.DatabitEQ,PAM_order);
    if eqtype==1
        eyeplot(Rx.outputeye.',sample_rate,spb,'Eyediagram of Data after FFE',PAM_order)
%         eye(2,:)=eyecomanalyze(Rx.DataEQcorrelate2,sample_rate,spb,PAM_order);
        histeye(Rx.outputeye,sample_rate,spb,'Eyediagram of Data after FFE',PAM_order)
        symbolhist(Rx.DatabitEQ,sample_rate,PAM_order,spb,'Symbol of Data after FFE',Rx.DatabitEQgroup)
%         eyeanaFFE=eyeanalyze(Rx.DataEQcorrelate2,Rx.Pam4levelmeanEQ,'after FFE',analyzesetting);
        Rx.Pam4levelmeanFFE=Rx.Pam4levelmeanEQ;
        Rx.Pam4levelstdFFE=Rx.Pam4levelstdEQ;
    elseif eqtype==2
        eyeplot(Rx.outputeye.',sample_rate,spb,'Eyediagram of Data after DFE',PAM_order)
%         eye(3,:)=eyecomanalyze(Rx.DataEQcorrelate,sample_rate,spb,PAM_order);
        histeye(Rx.outputeye,sample_rate,spb,'Eyediagram of Data after DFE',PAM_order)
        symbolhist(Rx.DatabitEQ,sample_rate,PAM_order,spb,'Symbol of Data after DFE',Rx.DatabitEQgroup)
%         eyeanaDFE=eyeanalyze(Rx.DataEQcorrelate,Rx.Pam4levelmeanEQ,'after DFE',analyzesetting);
        Rx.Pam4levelmeanDFE=Rx.Pam4levelmeanEQ;
        Rx.Pam4levelstdDFE=Rx.Pam4levelstdEQ;
    end

%% BER Estimate After EQ

QEQ=zeros(1,PAM_order-1);
BER.estEQ=zeros(1,PAM_order-1);

for i=1:PAM_order-1
    QEQ(1,i)=(Rx.Pam4levelmeanEQ(1,i+1)-Rx.Pam4levelmeanEQ(1,i))/(Rx.Pam4levelstdEQ(1,i+1)+Rx.Pam4levelstdEQ(1,i));
    BER.estEQ(1,i)=0.5*erfc(QEQ(1,i)/sqrt(2));
end

BER.esttotalEQ=sum(BER.estEQ);

%% BER Count After EQ

    [Rx.SymbolEQ]=BittoLevel(Rx.DatabitEQ,Rx.Pam4levelmeanEQ,PAM_order);
    BER.countEQ=sum(sum(dec2bin(Rx.SymbolEQ) ~= dec2bin(Tx.SymbolEQ(reflength+1:end))))/length(dec2bin(Rx.SymbolEQ))/log2(PAM_order);
    SNR.EQ=10*log10(sum(QEQ.^2));
    if eqtype==1
        SNR.QFFE=QEQ;
        SNR.FFED=10*log10((QEQ.^2));
        BER.FFE=BER.estEQ;
    elseif eqtype==2
        SNR.QDFE=QEQ;
        SNR.DFED=10*log10((QEQ.^2));
        BER.DFE=BER.estEQ;
    end
%     SNR.EQ = 10*log10(sum(Tx.TrainEQ(1+reflength:end).^2)/sum((Rx.DatabitEQ-Tx.TrainEQ(1+reflength:end)).^2));
    
if tapsweep==1;
    BERsweep.estEQ(ffen,dfen+1)=BER.esttotalEQ;
    BERsweep.countEQ(ffen,dfen+1)=BER.countEQ;
    BERsweep.SNREQ(ffen,dfen+1)=SNR.EQ;
end

    if eqtype==1
        SNR.FFE=SNR.EQ;
        BER.countFFE=BER.countEQ;
        BER.esttotalFFE=BER.esttotalEQ;
    elseif eqtype==2
        SNR.DFE=SNR.EQ;
        BER.countDFE=BER.countEQ;
        BER.esttotalDFE=BER.esttotalEQ;
    end
end

%% Output

% if tapsweep~=1
%     Result=[{'  '},{'SNR(dB)'},{'SNR +-10% (dB)'},{'BER count'},{'BER Q Factor'};
%         {'Without EQ'},SNR.woEQ,eyeanaorigin.SNR,BER.count,BER.esttotal;
%         {'With FFE'},SNR.FFE,eyeanaFFE.SNR,BER.countFFE,BER.esttotalFFE;
%         {'With DFE'},SNR.DFE,eyeanaDFE.SNR,BER.countDFE,BER.esttotalDFE;] 
% for i=PAM_order-1:-1:1
%     EyeAnalyze=[{'  '},{'eye amplitude'},{'eye height'},{'eye width'},{'  '};
%         {'Without EQ'},eyeanaorigin.amp2(1,i),eyeanaorigin.eyeheight(1,i),eyeanaorigin.eyewidth(1,i),{'  '};
%         {'With FFE'},eyeanaFFE.amp2(1,i),eyeanaFFE.eyeheight(1,i),eyeanaFFE.eyewidth(1,i),{'  '};
%         {'With DFE'},eyeanaDFE.amp2(1,i),eyeanaDFE.eyeheight(1,i),eyeanaDFE.eyewidth(1,i),{'  '};
%         {'  '},{'Eye Rise Time'},{'Eye Fall Time'},{'Eye Crossing Time'},{'Eye Crossing Amp'};
%         {'Without EQ'},eyeanaorigin.rise(1,i),eyeanaorigin.fall(1,i),eyeanaorigin.crosstime(1,i),eyeanaorigin.crossamp(1,i);
%         {'With FFE'},eyeanaFFE.rise(1,i),eyeanaFFE.fall(1,i),eyeanaFFE.crosstime(1,i),eyeanaFFE.crossamp(1,i);
%         {'With DFE'},eyeanaDFE.rise(1,i),eyeanaDFE.fall(1,i),eyeanaDFE.crosstime(1,i),eyeanaDFE.crossamp(1,i);] 
% end
% %         {'  '},{'Vertical eye open'},{'Horizontal eye open'},{'FECVertical eye open'},{'FEC Horizontal eye open'};
% %         {'Without EQ'},eye(1,1).vopen,eye(1,1).hopen,eye(1,1).vopenfec,eye(1,1).hopenfec;
% %         {'With FFE'},eye(2,1).vopen,eye(2,1).hopen,eye(2,1).vopenfec,eye(2,1).hopenfec;
% %         {'With DFE'},eye(3,1).vopen,eye(3,1).hopen,eye(3,1).vopenfec,eye(3,1).hopenfec;] 
% %     
% elseif tapsweep==1;
%     for i=1:FFEtap-1
%         for j=i+1:FFEtap
%             BERsweep.estEQ(i,j)=1;
%             BERsweep.countEQ(i,j)=1;
%             BERsweep.SNREQ(i,j)=1;
%         end
%     end
%     BERsweep.estEQMax = min(min(BERsweep.estEQ));
%     BERsweep.countEQMax = min(min(BERsweep.countEQ));
%     BERsweep.SNREQMax = max(max(BERsweep.SNREQ));
%     
%     Result=[{'  '},{'SNR(dB)'},{'BER count'},{'BER Q Factor'};
%         {'Without EQ'},SNR.woEQ,BER.count,BER.esttotal;
%         {'With LMS EQ'},BERsweep.SNREQMax,BERsweep.countEQMax,BERsweep.estEQMax ;] 
% end
    end
% %     OptFFEBER((fftap-5)/5+1,lmstime)=BER.esttotalFFE;
% %     OptFFEBER((n-12)/2+1,lmstime)=BER.esttotalFFE;
%     BeforeEQ((n-12)/2+1,1)=BER.esttotal;
%     OptFFEBER((n-12)/2+1,1)=BER.esttotalFFE;
%     OptDFEBER((n-12)/2+1,1)=BER.esttotalDFE;
%     BeforeEQ((power-0)/10+1,1)=BER.esttotal;
%     OptFFEBER((power-0)/10+1,1)=BER.esttotalFFE;
%     OptDFEBER((power-0)/10+1,1)=BER.esttotalDFE;
% % end
% end
% 
% save('OOK_BeforeEQ.txt','BeforeEQ','-ascii');
% save('OOK_FFE_5_powersweep.txt','OptFFEBER','-ascii')
% save('OOK_DFE_5_powersweep.txt','OptDFEBER','-ascii')

% save('PAM4_BeforeEQ.txt','BeforeEQ','-ascii');
% save('PAM4_FFE_5_powersweep.txt','OptFFEBER','-ascii')
% save('PAM4_DFE_5_powersweep.txt','OptDFEBER','-ascii')

% % 
% figure;
% grid on
% hold on
% % plot(dfesize,log10(OptDFEBER),'-o');
% plot(log10(OptFFEBER),'-o');
% plot(log10(OptDFEBER),'-o');
% % plot(log10(OptDFEBER),'-ro');
% xlabel ('lms step','FontSize',16,'FontWeight','bold');
% ylabel ('BER (dB)','FontSize',16,'FontWeight','bold');
% set(gca,'FontSize',16);
% set(gca,'FontWeight','bold');
% title('Optimized BER after EQ with different input power','FontSize',20,'FontWeight','bold');
% % % 
% % OptDFEBERoutput=[OptFFEBER OptDFEBER];

% save('PAM4_DFE_10_powersweep.txt','OptDFEBER','-ascii')