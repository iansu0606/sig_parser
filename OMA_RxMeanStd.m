function [OMA_Pam4level, Pam4level, group, OMA_index]=RxMeanStd(Databit,PAMorder,OMA_decision)

for i=1:PAMorder-1
    level(i)=(i-PAMorder/2)*2;
end

OMA_index = cell(PAMorder,1);
group = cell(PAMorder,1);
temp = zeros(1,length(Databit));

for bitn=1:length(Databit)
    for i=1:PAMorder-1
        if Databit(bitn)>level(i)
            temp(bitn)=i;
        end
    end
    group{temp(bitn)+1,1}=[group{temp(bitn)+1,1} {Databit(bitn)}];
    OMA_index{temp(bitn)+1,1} = [OMA_index{temp(bitn)+1,1} {OMA_decision(bitn)}];
end

for x=1:PAMorder
    OMA_Pam4levelmean(1,x) = mean(cell2mat(OMA_index{x,1})); 
    OMA_Pam4levelstd(1,x) = std(cell2mat(OMA_index{x,1})); 
    Pam4levelmean(1,x)=mean(cell2mat(group{x,1}));
    Pam4levelstd(1,x)=std(cell2mat(group{x,1}));
end
Pam4level = [Pam4levelmean; Pam4levelstd];
OMA_Pam4level = [OMA_Pam4levelmean; OMA_Pam4levelstd];
end