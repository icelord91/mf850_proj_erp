
function signal_gen()

    cids = unique(T.compid); % company id's
    % leverage example
    y = sortrows(T(T.compid == 1, :), 'Date');
    
    
    
end