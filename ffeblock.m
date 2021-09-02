function [DataFFE,data_rx_eq,errCount]=ffeblock(inputdata,Ref,PAM_order,lmssize,fftap,dfetap)

SigConst = pammod(0:PAM_order-1,PAM_order);                                 % signal constellation
eqlms = dfe(fftap+1,dfetap,lms(lmssize));
eqlms.SigConst = SigConst;                                                  % Set signal constellation.
eqlms.Weights = [1 zeros(1,fftap)];
eqlms.ResetBeforeFiltering = 0;  

for i=1:5
[DataFFE data_rx_eq errCount]= equalize(eqlms,inputdata,Ref);
end
% [DataFFE data_rx_eq errCount]= equalize(eqlms,inputdata);



x1=0:1:fftap;
x2=0:0.1:fftap;
fitcurve=interp1(x1,eqlms.Weights,x2,'cubic');


figure;
hold on
grid on
stem(x1,eqlms.Weights,'--o','Linewidth',3,'Markersize',6)
plot(x2,fitcurve,':r','Linewidth',4)
axis([0 fftap -0.2 1.2])
xlabel ('Number of taps','FontSize',14,'FontWeight','bold');
ylabel ('Weight','FontSize',14,'FontWeight','bold');
set(gca,'FontSize',14);
set(gca,'FontWeight','bold');
title('Tap Weights of FFE','FontSize',20,'FontWeight','bold');