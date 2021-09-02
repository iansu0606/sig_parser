function [Pam4levelmean Pam4levelstd group]=RxMeanStd(Databit,PAMorder)

for i=1:PAMorder-1
    level(i)=(i-PAMorder/2)*2;
end

group=cell(PAMorder,1);
temp=zeros(1,length(Databit));

for bitn=1:length(Databit)
    for i=1:PAMorder-1
        if Databit(bitn)>level(i)
            temp(bitn)=i;
        end
    end
    group{temp(bitn)+1,1}=[group{temp(bitn)+1,1} {Databit(bitn)}];
end

for x=1:PAMorder
    Pam4levelmean(1,x)=mean(cell2mat(group{x,1}));
    Pam4levelstd(1,x)=std(cell2mat(group{x,1}));
end


% function [Pam4levelmean Pam4levelstd group]=RxMeanStd(Databit,Train,PAMorder)
% 
% group=cell(PAMorder,1);
% 
% for bitn=1:length(Databit)
%     group{Train(bitn)+1,1}=[group{Train(bitn)+1,1} {Databit(bitn)}];
% end
% 
% for x=1:PAMorder
%     Pam4levelmean(1,x)=mean(cell2mat(group{x,1}));
%     Pam4levelstd(1,x)=std(cell2mat(group{x,1}));
% end