function [DataDFE]=dfeblock(inputdata,Ref,PAM_order,lmssize2,fftap,dfetap)

SigConst = pammod(0:PAM_order-1,PAM_order);                                 % signal constellation
eqlms = dfe(fftap+1,dfetap,lms(lmssize2));
eqlms.SigConst = SigConst;                                                  % Set signal constellation.
eqlms.Weights = [1 zeros(1,fftap) zeros(1,dfetap)];
eqlms.ResetBeforeFiltering = 0;

for q=1:5
    [DataDFE data_rx_eq errCount]= equalize(eqlms,inputdata,Ref);
end
tapweight=eqlms.Weights;
tapweight(1,dfetap+1:end)=eqlms.Weights(1,1:fftap+1);
for i=1:dfetap
    tapweight(1,i)=eqlms.Weights(1,end-i+1);
end

x1=0-dfetap:1:fftap;
x2=0-dfetap:0.1:fftap;
fitcurve=interp1(x1,tapweight,x2,'cubic');

figure;
hold on
grid on
% stem(eqlms.Weights,'--o','Linewidth',2,'Markersize',6)
stem(x1,tapweight,'--o','Linewidth',3,'Markersize',6)
plot(x2,fitcurve,':r','Linewidth',4)
axis([-dfetap fftap -0.2 1.2])
xlabel ('Number of taps','FontSize',14,'FontWeight','bold');
ylabel ('Weight','FontSize',14,'FontWeight','bold');
set(gca,'FontSize',14);
set(gca,'FontWeight','bold');
title('Tap Weights of DFE','FontSize',20,'FontWeight','bold');