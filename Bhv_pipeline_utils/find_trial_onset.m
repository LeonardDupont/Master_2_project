function [idx] = find_trial_onset(time,fs)

    difft = zeros(length(time),1);
    for k = 1:length(time)
        if k==1
            difft(k) = -time(k);
        else
            difft(k) = time(k) - time(k-1); 
        end
    end

    idx = find(abs(difft) > 2/fs); 
end