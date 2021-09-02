function [PeakMax, ShiftLength]=headcorr(inputdata,ref,correlatelength)

FindDelay = ref(1:correlatelength);
for start = 1:length(inputdata)-correlatelength
    y = corrcoef(inputdata(start:start-1+correlatelength),FindDelay);
    SyncChoice(start,1) = y(2,1);
end
[PeakMax, ShiftLength] = max(SyncChoice);
end



%%  Test Bench     
%{
inputdata = Rx.Data;
ref =Tx.Ref;
FindDelay = ref(1:correlatelength);
for start = 1:length(inputdata)-correlatelength
    y = corrcoef(inputdata(start:start-1+correlatelength),FindDelay);
    SyncChoice(start,1) = y(2,1);
end
[PeakMax, ShiftLength] = max(SyncChoice);
%}