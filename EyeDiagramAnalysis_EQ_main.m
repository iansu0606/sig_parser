clear all
clc
close all

%% Switching
pulsetest=0;
ploten=1;
filter_type=0; % raise cosine
eqenable=1;
fullwave=0;
output=1;
volterra_filter=0;                                                         % volterra filter switch
%% Signal setting
PAM_order = 4;
PRBS_length = 11;
if pulsetest==1
    PAM_order = 2;
    PRBS_length = 5;
end
spb = 32;                                                                       % sample per symbol
symbol_rate = 53.125e9;                                                             % bit_rate/2 = baud_rate
sample_rate = symbol_rate*spb;                                                  % bit_rate*sample_count_per_bit
correlatelength=500;
upsamplenumber=1;
eyenumber=2;
calculatelength=2^PRBS_length;
basicsetting=[symbol_rate,spb,PAM_order,upsamplenumber,eyenumber,correlatelength,PRBS_length];

%% EQ setting
fftap=3;
dfetap=1;
if PAM_order==2
    lmssize=0.08;
    lmssize2=0.05;
elseif PAM_order==4
    lmssize=0.001;
    lmssize2=lmssize/10;    
end

%% Volterra Setting
reflength=100;
weighting_nunber=100;
tap_number=7;                                                             %volterra tap number
tap_order=3;                                                               % order of volterra filter
tap=tap_number;                                                            % memory = tap number-1 1st
tap_2=tap;                                                                 % memory = tap number-1 2nd
tap_3=tap;                                                                 % memory = tap number-1 3rd
third_term_type=1;                                                         % 1是3階全部 2是3階只有自己  1 is type 2, 2 is type1
repeat_time=1;

%% Analysis setting
crosstimeprecent=10;
risefallprecent=10;
BERthrehold=2e-4;
SNRprecent=10;
analyzesetting=[sample_rate,spb,PAM_order,crosstimeprecent,risefallprecent,BERthrehold,SNRprecent];

%% Load file

if PAM_order==2
    Txpower=load('I:\HW_ted_105_2\HW\tx_dfe_1.txt');
    Rxpower=load('I:\HW_ted_105_2\HW\rx_dfe_1.txt');
    Pam4Level1 = 0.03;
    Pam4Level2 = 0.06;
    Tx.Pam4level=[Pam4Level1 Pam4Level2];
elseif PAM_order==4
%     Txpower=load('D:\1_Ted\Delta_simul\20201127\DR4_4Sig_tx01_real.sig');           %Tx signal_c1
%     Rxpower=load('D:\1_Ted\Delta_simul\20201127\DR4_4Sig_rx3_real.sig');           %Rx signal%sig
    Txpower=load('D:\1_Ted\Delta_simul\20210524_DML\CWDMTx_1_real.sig');           %Tx signal_c1
    Rxpower=load('D:\1_Ted\Delta_simul\20210524_DML\CWDMSigPlt2_real.sig');           %Rx signal%sig
%% VCSEL
     %Txpower=load('D:\1_Ted\Delta_simul\20210308\lstmpTx_1_real.sig');           %Tx signal_c1
    % Rxpower=load('D:\1_Ted\Delta_simul\20210308\lstmpRx_1_real.sig');           %Rx signal%sig
    %Rxpower= Rxpower(:,2);
    Pam4Level1 = -0.6;
    Pam4Level2 = -0.2;
    Pam4Level3 = 0.2;
    Pam4Level4 = 0.6;
    Tx.Pam4level=[Pam4Level1 Pam4Level2 Pam4Level3 Pam4Level4];
end
Rx.Data = Rxpower(:,2);
Tx.Ref = Txpower(:,2);
Rx.power=-4;
ER=6.5;
OMA=Rx.power-10*log10(0.5*(ER+1)/(ER-1));

%% Upsample

x=1:length(Tx.Ref);
x2=1:1/upsamplenumber:(length(Tx.Ref)+1-1/upsamplenumber);
Tx.Ref=interp1(x,Tx.Ref,x2,'cubic').';tap_number
Rx.Data=interp1(x,Rx.Data,x2,'cubic').';

%% raise cosine filter

    if filter_type ~= 0
        matching_filter = rcosfir(0.1,5,spb,1,'sqrt');
        data_scope_up_temp2 = conv(Rx.Data,matching_filter);
        Rx.Datafilter = data_scope_up_temp2((length(matching_filter)-1)/2+1:end-(length(matching_filter)-1)/2);
    else
        Rx.Datafilter = Rx.Data;
    end
    
%% Correlation
tic;

[PeakMax ShiftLength]=headcorr(Rx.Data,Tx.Ref,correlatelength);
if length(Rx.Data)>=calculatelength*spb*upsamplenumber
    Tx.Refcorrelate=Tx.Ref(1:(calculatelength-ceil(ShiftLength/spb/upsamplenumber))*spb*upsamplenumber,1);
    Rx.Datacorrelate=Rx.Data(ShiftLength:length(Tx.Refcorrelate)+ShiftLength-1);
else
    Rx.Datacorrelate=Rx.Data(ShiftLength:end-ceil((ShiftLength-1)/spb/upsamplenumber)*spb*upsamplenumber+(ShiftLength-1));
    Tx.Refcorrelate=Tx.Ref(1:length(Rx.Datacorrelate),1);
    fprintf('1\n');
end

%% Signal normalize

[Tx.normalized,Tx.norfactor]=datanormalize(Tx.Refcorrelate,PAM_order,pulsetest);
[Rx.normalized,Rx.norfactor]=datanormalize(Rx.Datacorrelate,PAM_order,pulsetest);

%% Find min BER pos

Tx.Refreshape = reshape(Tx.Refcorrelate,spb*upsamplenumber,[]);
Rx.Datareshape = reshape(Rx.normalized,spb*upsamplenumber,[]);
[eyeopenmin eyeopenpos]=BERsweep(Rx.Datareshape,basicsetting);

%% Choose Decision Points

Tx.Refbit = Tx.Refreshape(floor(spb*upsamplenumber/2)+1,:);                                         % symbol select
eyeopenpos=16;
Rx.Databit = Rx.Datareshape(eyeopenpos,:);

%% BER Estimate w/o EQ

    [Rx.Pam4levelmean Rx.Pam4levelstd Rx.Databitgroup]=RxMeanStd(Rx.Databit,PAM_order);
    [BER.esttotal,BER.est,SNR.Q]=QBER(Rx.Pam4levelmean,Rx.Pam4levelstd,basicsetting);
    SNR.woEQ=10*log10(sum(SNR.Q.^2));
    SNR.woEQD=10*log10((SNR.Q.^2));
%% TDECQ calculation
Qt=3.414;
BER_q=2.64E-07;
a=(erfcinv(BER_q*8/3)*(sqrt(2)/Qt)).^2;
TDECQ_q=sqrt(a./(a-1));
TDECQ_OFC=10*log10(TDECQ_q)
% %
% sigmal_simulate=0.04;
% sigma_ideal=0.0586;
% sigma_added=sqrt(sigma_ideal^2-sigmal_simulate^2);
% TDECQ_teacher=10*log10(sigma_ideal/sigma_added);
%% BER Count w/o EQ

    [Tx.Symbol]=BittoLevel(Tx.Refbit,Tx.Pam4level,PAM_order);
    [Rx.Symbol]=BittoLevel(Rx.Databit,Rx.Pam4levelmean,PAM_order);
    BER.count=sum(sum(dec2bin(Rx.Symbol,PAM_order) ~= dec2bin(Tx.Symbol,PAM_order)))/length(dec2bin(Rx.Symbol,PAM_order))/log2(PAM_order);

%% Eye Diagram of Data
 
if ploten==1
    Tx.refeye=Tx.normalized(floor(spb*upsamplenumber/2)+1:end);
    Rx.inputeye=Rx.normalized(1+floor(spb*upsamplenumber/2):end-ceil(spb*upsamplenumber/2));
    if mod(length(Tx.refeye),eyenumber*spb*upsamplenumber)~=0
        Tx.refeye=Tx.refeye(1:end-mod(length(Tx.refeye),eyenumber*spb*upsamplenumber));
    end
    if mod(length(Rx.inputeye),eyenumber*spb*upsamplenumber)~=0
        Rx.inputeye=Rx.inputeye(1:end-mod(length(Rx.inputeye),eyenumber*spb*upsamplenumber));
    end
    Tx.refeye=reshape(Tx.refeye,eyenumber*spb*upsamplenumber,[]);
    Rx.inputeye=reshape(Rx.inputeye,eyenumber*spb*upsamplenumber,[]);
    eyeplot(Tx.refeye.',basicsetting,'Transceiver Eye',Tx.norfactor)
    eyeplot(Rx.inputeye.',basicsetting,'Eyediagram of Data before EQ',Rx.norfactor)
    histeye(Rx.inputeye,basicsetting,'Eyediagram of Data before EQ',[10,8,100],Rx.norfactor,0,90)
    symbolhist(Rx.Databit,basicsetting,'Symbol of Data before EQ',Rx.Databitgroup,Rx.norfactor)
end
% pause
%% EQ

if eqenable==1
    if fullwave==1
        FFERef=[];
        [Rx.DataFFE data_rx_eq errCount]=ffefullwave(Rx.Databit,Rx.normalized,FFERef,PAM_order,lmssize,fftap,0,basicsetting);
        [BER.FFE,SNR.FFE,Rx.Pam4levelmeanFFE,Rx.Pam4levelstdFFE]=EQanalysis(Rx.DataFFE,basicsetting,Tx.Ref,0,Tx.Symbol,ploten,'FFE');
        
        DFERef=[];
        Rx.DataDFE=dfeblockfullwave(Rx.Databit,Rx.normalized,DFERef,PAM_order,lmssize2,fftap,dfetap,basicsetting);
        [BER.DFE,SNR.DFE,Rx.Pam4levelmeanDFE,Rx.Pam4levelstdDFE]=EQanalysis(Rx.DataDFE,basicsetting,Tx.Ref,0,Tx.Symbol,ploten,'DFE');
    else
        Ref=Tx.Ref(1+floor(spb*upsamplenumber/2):correlatelength+floor(spb*upsamplenumber/2));
        FFERef=[];
        [Rx.DataFFE data_rx_eq errCount]=ffeblock(Rx.Databit,FFERef,PAM_order,lmssize,fftap,0);
        Rx.DataFFEre=resample(Rx.DataFFE,sample_rate*upsamplenumber,symbol_rate).';
        [BER.FFE,SNR.FFE,Rx.Pam4levelmeanFFE,Rx.Pam4levelstdFFE]=EQanalysis(Rx.DataFFEre,basicsetting,Ref,0,Tx.Symbol,ploten,'FFE');

        DFERef=[];
        Rx.DataDFE=dfeblock(Rx.Databit,DFERef,PAM_order,lmssize2,fftap,dfetap);
        Rx.DataDFEre=resample(Rx.DataDFE,sample_rate*upsamplenumber,symbol_rate).';
        [BER.DFE,SNR.DFE,Rx.Pam4levelmeanDFE,Rx.Pam4levelstdDFE]=EQanalysis(Rx.DataDFEre,basicsetting,Ref,0,Tx.Symbol,ploten,'DFE');
    end
end

%% Volterra

if volterra_filter==1
    Ref=Tx.Ref(1+floor(spb*upsamplenumber/2):correlatelength+floor(spb*upsamplenumber/2));
    SymL=spb;
    T_SymC=reflength*upsamplenumber;
    D_SymC=length(Rx.normalized)/spb-T_SymC;
    txTrainSig=Tx.normalized(1:spb*upsamplenumber*reflength).';
    RX_data_c=Rx.normalized;
    RX_data=Rx.normalized;
    [RX_in_wiener,w,R_det]=RXwiener_v2star(RX_data_c,RX_data,SymL,txTrainSig,tap_order,tap,tap_2,tap_3,third_term_type,T_SymC,D_SymC,repeat_time,weighting_nunber);
    Rx.vot=RX_in_wiener(1+floor(spb*upsamplenumber/2):end-ceil(spb*upsamplenumber/2));
    [BER.Vol,SNR.Vol,Rx.Pam4levelmeanVol,Rx.Pam4levelstdVol]=EQanalysis(Rx.vot,basicsetting,Ref,0,Tx.Symbol(1:end-1),ploten,'Vol');
end

%% Output
if output==1
    Result=[{'  '},{'S/N(dB)'},{'BER count'},{'BER Q Factor'};
        {'Without EQ'},SNR.Q,BER.count,BER.esttotal;
        {'With FFE'},SNR.FFE.Q,BER.FFE.countEQ,BER.FFE.esttotalEQ;
        {'With DFE'},SNR.DFE.Q,BER.DFE.countEQ,BER.DFE.esttotalEQ;]
%         {'With Volt'},SNR.Vol.EQ,BER.Vol.countEQ,BER.Vol.esttotalEQ; 
end
%% plot spec - Ted
% signalfft(Rx.normalized,sample_rate)
% %
% A=[0.2 0.8 -0.1 0.1 0];
% B=1;
% w = -pi:pi/10000:pi;
% H = freqz(A,B,w);
% subplot(211)
% plot(w,abs(H))
% axis([0 pi 0 1.5]); grid;
% ylabel('Magnitude Response')
% subplot(212)
% plot(w,angle(H))
% axis([-pi pi -pi pi]); grid;
% ylabel('Phase Response (rad)')
% xlabel('hat(\omega)')
%
% Y = filter(B,A,Rx.normalized);
% figure;plot(Y)
%% esimate eye linearity
V_averg_noEQ=1/4*sum(Rx.Pam4levelmean);
ES_noEQ_1=(Rx.Pam4levelmean(2)-V_averg_noEQ)/(Rx.Pam4levelmean(1)-V_averg_noEQ);
ES_noEQ_2=(Rx.Pam4levelmean(3)-V_averg_noEQ)/(Rx.Pam4levelmean(4)-V_averg_noEQ);
S_noEQ_min=1/2*min([Rx.Pam4levelmean(4)-Rx.Pam4levelmean(3),Rx.Pam4levelmean(3)-...
    Rx.Pam4levelmean(2),Rx.Pam4levelmean(2)-Rx.Pam4levelmean(1)]);
Rlm_noEQ=6*S_noEQ_min/(Rx.Pam4levelmean(4)-Rx.Pam4levelmean(1));
eye_noEQ_height_1=Rx.Pam4levelmean(2)-Rx.Pam4levelmean(1);
eye_noEQ_height_2=Rx.Pam4levelmean(3)-Rx.Pam4levelmean(2);
eye_noEQ_height_3=Rx.Pam4levelmean(4)-Rx.Pam4levelmean(3);
eye_noEQ_height=[eye_noEQ_height_1,eye_noEQ_height_2,eye_noEQ_height_3]
%
V_averg_FFE=1/4*sum(Rx.Pam4levelmeanFFE);
ES_FFE_1=(Rx.Pam4levelmeanFFE(2)-V_averg_FFE)/(Rx.Pam4levelmeanFFE(1)-V_averg_FFE);
ES_FFE_2=(Rx.Pam4levelmeanFFE(3)-V_averg_FFE)/(Rx.Pam4levelmeanFFE(4)-V_averg_FFE);
S_FFE_min=1/2*min([Rx.Pam4levelmeanFFE(4)-Rx.Pam4levelmeanFFE(3),Rx.Pam4levelmeanFFE(3)-...
    Rx.Pam4levelmeanFFE(2),Rx.Pam4levelmeanFFE(2)-Rx.Pam4levelmeanFFE(1)]);
Rlm_FFE=6*S_FFE_min/(Rx.Pam4levelmeanFFE(4)-Rx.Pam4levelmeanFFE(1));
eye_FFE_height_1=Rx.Pam4levelmeanFFE(2)-Rx.Pam4levelmeanFFE(1);
eye_FFE_height_2=Rx.Pam4levelmeanFFE(3)-Rx.Pam4levelmeanFFE(2);
eye_FFE_height_3=Rx.Pam4levelmeanFFE(4)-Rx.Pam4levelmeanFFE(3);
eye_FFE_height=[eye_FFE_height_1,eye_FFE_height_2,eye_FFE_height_3]
%
V_averg_DFE=1/4*sum(Rx.Pam4levelmeanDFE);
ES_DFE_1=(Rx.Pam4levelmeanDFE(2)-V_averg_DFE)/(Rx.Pam4levelmeanDFE(1)-V_averg_DFE);
ES_DFE_2=(Rx.Pam4levelmeanDFE(3)-V_averg_DFE)/(Rx.Pam4levelmeanDFE(4)-V_averg_DFE);
S_DFE_min=1/2*min([Rx.Pam4levelmeanDFE(4)-Rx.Pam4levelmeanDFE(3),Rx.Pam4levelmeanDFE(3)-...
    Rx.Pam4levelmeanDFE(2),Rx.Pam4levelmeanDFE(2)-Rx.Pam4levelmeanDFE(1)]);
Rlm_DFE=6*S_DFE_min/(Rx.Pam4levelmeanDFE(4)-Rx.Pam4levelmeanDFE(1));
eye_DFE_height_1=Rx.Pam4levelmeanDFE(2)-Rx.Pam4levelmeanDFE(1);
eye_DFE_height_2=Rx.Pam4levelmeanDFE(3)-Rx.Pam4levelmeanDFE(2);
eye_DFE_height_3=Rx.Pam4levelmeanDFE(4)-Rx.Pam4levelmeanDFE(3);
eye_DFE_height=[eye_DFE_height_1,eye_DFE_height_2,eye_DFE_height_3]
%
Rlm_total=[Rlm_noEQ,Rlm_FFE,Rlm_DFE]
%% TDECQ after EQ
Qt=3.414;
BER_q=BER.DFE.esttotalEQ;
a=(erfcinv(BER_q*8/3)*(sqrt(2)/Qt)).^2;
TDECQ_q=sqrt(a./(a-1));
TDECQ_OFC=10*log10(TDECQ_q)
 Result=[{'  '},{'S/N(dB)'},{'BER count'},{'BER Q Factor'};
        {'Without EQ'},SNR.Q,BER.count,BER.esttotal;
        {'With FFE'},SNR.FFE.Q,BER.FFE.countEQ,BER.FFE.esttotalEQ;
        {'With DFE'},SNR.DFE.Q,BER.DFE.countEQ,BER.DFE.esttotalEQ;]
 %%
 Rx.Pam4levelmeanDFE
 Rx.Pam4levelstdDFE
 %% cascade BW fmax=0.24/T, T is rise time from 20 to 80%
 a=0.24;
 f1=40e9;
 f2=40e9;
 f3=40e9;
 T1=a/f1;
 T2=a/f2;
 T3=0;%a/f3;
 T=sqrt(T1^2+T2^2+T3^2);
 f_cascade=a/T;
 