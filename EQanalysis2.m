function [BER,SNR,Pam4levelmeanEQ,Pam4levelstdEQ]=EQanalysis(inputdata,setting,ref,reflength,TxSymbol,ploten,name)

symbol_rate=setting(1,1);
spb=setting(1,2);
PAM_order=setting(1,3);
upsamplenumber=setting(1,4);
eyenumber=setting(1,5);
correlatelength=setting(1,6);


%% recorrelate

[PeakMax ShiftLength2]=headcorr(inputdata,ref,correlatelength);

Rx.DataEQcorrelate=inputdata(ShiftLength2:end-ceil((ShiftLength2-1)/spb/upsamplenumber)*spb*upsamplenumber+ShiftLength2-1,1);

%% Find min BER pos

Rx.DataEQreshape = reshape(Rx.DataEQcorrelate,spb*upsamplenumber,[]);    
[eyeopenmin eyeopenposEQ]=BERsweep(Rx.DataEQreshape,setting);

if eyeopenposEQ>floor(spb*upsamplenumber/2)+1
    countstart=1;
else
    countstart=0;
end

%% Choose Decision Points

    Tx.SymbolEQ=TxSymbol(reflength+1:end-ceil((ShiftLength2-1)/spb/upsamplenumber)*spb*upsamplenumber);
    Rx.DatabitEQ=Rx.DataEQreshape(eyeopenposEQ,:);
    Rx.DatabitEQ=Rx.DatabitEQ(1,reflength+1:end-1);

%% BER Estimate After EQ
    [Pam4levelmeanEQ Pam4levelstdEQ Rx.DatabitEQgroup]=RxMeanStd(Rx.DatabitEQ,PAM_order);
    [BER.esttotalEQ,BER.estEQ,SNR.Q]=QBER(Pam4levelmeanEQ,Pam4levelstdEQ,setting);
    SNR.EQ=10*log10(sum(SNR.Q.^2));
    SNR.D=10*log10((SNR.Q.^2));
        
%% BER Count After EQ
    
    [Rx.SymbolEQ]=BittoLevel(Rx.DatabitEQ,Pam4levelmeanEQ,PAM_order);
    BER.countEQ=sum(sum(dec2bin(Rx.SymbolEQ(1:end)) ~= dec2bin(Tx.SymbolEQ(1+countstart:countstart+length(Rx.SymbolEQ)))))/length(dec2bin(Rx.SymbolEQ))/log2(PAM_order);
    
%% Eye Diagram
if ploten==1
if mod(length(Rx.DataEQcorrelate),eyenumber*spb*upsamplenumber)~=0
    Rx.outputeye=Rx.DataEQcorrelate(1:end-mod(length(Rx.DataEQcorrelate),eyenumber*spb*upsamplenumber));
else
    Rx.outputeye=Rx.DataEQcorrelate;
end
    Rx.outputeye=reshape(Rx.outputeye,eyenumber*spb*upsamplenumber,[]);
    name1=['Eyediagram of Data after ',name];
    name2=['Symbol of Data after ',name];
    eyeplot(Rx.outputeye.',setting,name1,[1 0])
    histeye(Rx.outputeye,setting,name1,10,8,[1 0],0,90)
    symbolhist(Rx.DatabitEQ,setting,name2,Rx.DatabitEQgroup,[1 0])
end