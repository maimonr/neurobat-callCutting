function merged = merge_wins(wins,fs,thresh)
if size(wins,1)>1
    diffs = (wins(2:length(wins),1)- wins(1:length(wins)-1,2));
    if any(diffs/fs < thresh)
        w = find(diffs/fs < thresh,1);
        wins(w,:) = [wins(w,1),wins(w+1,2)];
        wins_del = 1:length(wins) ~= w+1;
        merged = wins(wins_del,:);
        merged = merge_wins(merged,fs,thresh);
    else
        merged = wins;
    end
else
    merged = wins;
end

end

    