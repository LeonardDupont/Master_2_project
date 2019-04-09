function [ISIs] = simple_ISIs(binaryspikes)

    
    wherearethey = find(binaryspikes == 1); 
    S = length(wherearethey);
    ISIs = zeros(S-1,1);
    for k = 1:S-1
        ISIs(k) = 30/(wherearethey(k+1)-wherearethey(k)); 
    end


end