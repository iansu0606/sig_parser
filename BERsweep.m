function [eyeopenmin, eyeopenpos]=BERsweep(inputdata,setting)

spb=setting(1,2);
PAM_order=setting(1,3);
upsamplenumber=setting(1,4);

Qsweep=zeros(spb*upsamplenumber,PAM_order-1);
BER.estsweep=zeros(spb*upsamplenumber,PAM_order-1);
BER.esttotalsweep=zeros(1,spb*upsamplenumber);

for eyeopenpos=1:spb*upsamplenumber

    Rx.Databit = inputdata(eyeopenpos,:);
    [Rx.Pam4levelmean Rx.Pam4levelstd Rx.Databitgroup]=RxMeanStd(Rx.Databit,PAM_order);

    for i=1:PAM_order-1
        Qsweep(eyeopenpos,i)=(Rx.Pam4levelmean(1,i+1)-Rx.Pam4levelmean(1,i))/(Rx.Pam4levelstd(1,i+1)+Rx.Pam4levelstd(1,i));
        BER.estsweep(eyeopenpos,i)=0.5*erfc(Qsweep(eyeopenpos,i)/sqrt(2));
    end
    
    BER.esttotalsweep(1,eyeopenpos)=sum(BER.estsweep(eyeopenpos,:));
end

[eyeopenmin, eyeopenpos]=min(BER.esttotalsweep);

end