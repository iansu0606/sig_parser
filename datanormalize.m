function [output,factor]=datanormalize(inputdata,PAM_order,pulsetest)
if pulsetest~=1
    factor_mean=mean(inputdata);
    datashift=inputdata-factor_mean;
    factor_amp=mean(abs(datashift))/log2(PAM_order);
    output=datashift/factor_amp;
    factor=[factor_amp factor_mean];
else
    factor_amp=(max(inputdata)-min(inputdata))/2;
    factor_mean=(max(inputdata)+min(inputdata))/2;
    datashift=inputdata-factor_mean;
    output=datashift/factor_amp;
    factor=[factor_amp factor_mean];
end