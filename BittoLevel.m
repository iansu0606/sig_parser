function [symbol]=BittoLevel(Refbit,Pam4level,PAMorder)

level=zeros(1,PAMorder);

for i=1:PAMorder-1
    level(i)=(i-1);
    levelmean(1,i)=(Pam4level(1,i)+Pam4level(1,i+1))/2;
end
    level(PAMorder)=PAMorder-1;
    
 for i=1:length(Refbit);
    for j=1:PAMorder-1
        if Refbit(1,i) >= levelmean(1,j)
            symbol(1,i) = level(j+1);
        elseif Refbit(1,i) < levelmean(1,1)
            symbol(1,i) = level(1);
        end
    end
 end