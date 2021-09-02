function [esttotal,est,Q]=QBER(Pam4levelmean,Pam4levelstd,setting)
%BER=0.5*erfc(Q/sqrt(2))
PAM_order=setting(1,3);

Q=zeros(1,PAM_order-1);
BER.est=zeros(1,PAM_order-1);

for i=1:PAM_order-1
    Q(1,i)=(Pam4levelmean(1,i+1)-Pam4levelmean(1,i))/(Pam4levelstd(1,i+1)+Pam4levelstd(1,i));
    est(1,i)=0.5*erfc(Q(1,i)/sqrt(2));
end

esttotal=sum(est);

