% FFE response
clear all
close all
%
w = -pi:pi/100:pi;
H = freqz([-0.2 0.1 0.8 -0.2 -0.1],1,w);
subplot(211)
plot(w,abs(H))
axis([0 pi 0 1.5]); grid;
ylabel('Magnitude Response')
subplot(212)
plot(w,angle(H))
axis([-pi pi -pi pi]); grid;
ylabel('Phase Response (rad)')
xlabel('hat(\omega)')