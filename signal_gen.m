
% function signals = signal_gen(T)
    
    cids = unique(T.compid); % company id's
    shape = [height(T), 1]; % signal shape
    
%     ver = str2double(regexp(version, '^\d+.\d+', 'match')); % check matlab version
    
    days = unique(T.Date);
    H = length(days);
    W = width(T);
    template = array2table(nan([H, W]));
    template.Properties.VariableNames = T.Properties.VariableNames;
    template.Date = days;
    template.Industry = repmat({''}, [H,1]);
    template.Sector = repmat({''}, [H,1]);
    
    %% momentum
    mom1m = nan(shape);
    func = @(ts) erplib.shift(log(1 + ts.RETMONTH) - log(1 + ts.retmonth_spx), 1);
    tmpl = template;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date); % aligning missing dates
        tmpl(pos, :) = ts;
        res = func(tmpl);
        mom1m(loc) = res(pos, :);
    end;
    mom1m = mom1m / 0.5776;
    
    mom3m = nan(shape);
    func = @(ts) erplib.shift(erplib.ma(log(1 + ts.RETMONTH) - log(1 + ts.retmonth_spx), 3), 1);
    tmpl = template;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date);
        tmpl(pos, :) = ts;
        res = func(tmpl);
        mom3m(loc) = res(pos, :);
    end;
    mom3m = mom3m / 0.3538;
    
    mom6m = nan(shape);
    func = @(ts) erplib.shift(erplib.ma(log(1 + ts.RETMONTH) - log(1 + ts.retmonth_spx), 6), 1);
    tmpl = template; 
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date);
        tmpl(pos, :) = ts;
        res = func(tmpl);
        mom6m(loc) = res(pos, :);
    end;
    mom6m = mom6m / 0.2711;
    
    mom12m = nan(shape);
    func = @(ts) erplib.shift(erplib.ma(log(1 + ts.RETMONTH) - log(1 + ts.retmonth_spx), 12), 1);
    tmpl = template;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date);
        tmpl(pos, :) = ts;
        res = func(tmpl);
        mom12m(loc) = res(pos, :);
    end;
    mom12m = mom12m / 0.2169;
    
    chmom = nan(shape);
    func = @(ts) erplib.ma(ts.mom1m, 3) -  erplib.ma(ts.mom1m, 6);
    days = unique(T.Date);
    H = length(days);
    W = 2;
    tmpl = array2table(nan([H, W]));
    tmpl.Properties.VariableNames = {'Date', 'mom1m'};
    tmpl.Date = days;
    M = table(T.Date, mom1m, 'VariableNames', {'Date', 'mom1m'});
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = M(loc, :);
        pos = ismember(tmpl.Date, ts.Date);
        tmpl(pos, :) = ts;
        res = func(tmpl);
        chmom(loc) = res(pos, :);
    end;
    chmom = chmom / 0.4228;
    
    indmom = nan(shape);
    days = unique(T.Date);
    sectors = unique(T.Sector);
    H = length(days);
    W = 1 + length(sectors);
    S = array2table(nan([H, W]));
    S.Properties.VariableNames = {'Date', sectors{:}};
    S.Date = days;
    M = table(T.Date, mom1m, 'VariableNames', {'Date', 'mom1m'}); % industry benchmark
    for j = 1:length(sectors) % cross-sectional calculation
        block = M(strcmp(T.Sector, sectors(j)), :);
        S(:, sectors(j)) = ...
            array2table(arrayfun(@(day) nanmean(block(block.Date == day, :).mom1m), days));
%         S(:, sectors(j)).Variables = ...
%             arrayfun(@(day) nanmean(block(block.Date == day, 'mom1m').Variables), days);
    end;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(S.Date, ts.Date);
        res = S(:, ts.Sector(1));
        indmom(loc) = table2array(res(pos, :));
%         res = S(:, ts.Sector(1)).Variables;
%         indmom(loc) = res(pos, :);
    end;
	indmom = indmom / 0.2642;
    
    
    %% trend
    beta = nan(shape);
    func = @(ts) erplib.shift( ...
        (erplib.ma(ts.RETMONTH .* ts.retmonth_spx, 6) ...
            ./ erplib.std(ts.RETMONTH, 6) ...
            ./ erplib.std(ts.retmonth_spx, 6)), 1); % biased
    tmpl = template; 
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date);
        tmpl(pos, :) = ts;
        res = func(tmpl);
        beta(loc) = res(pos, :);
    end;
    beta = real(beta);
    beta(~isfinite(beta), :) = nan;
    beta = beta / 2.3904;
    
    
    %% volatility
    gvol = nan(shape);
    func = @(ts) log(erplib.ma(ts.mom1m .^ 2, 3));
    days = unique(T.Date);
    H = length(days);
    W = 2;
    tmpl = array2table(nan([H, W]));
    tmpl.Properties.VariableNames = {'Date', 'mom1m'};
    tmpl.Date = days;
    M = table(T.Date, mom1m, 'VariableNames', {'Date', 'mom1m'});
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = M(loc, :);
        pos = ismember(tmpl.Date, ts.Date);
        tmpl(pos, :) = ts;
        res = func(tmpl);
        gvol(loc) = res(pos, :);
    end;
    gvol = gvol / 5.8917;
    
    idiovol = nan(shape);
    func = @(ts) log(erplib.ma((ts.mom1m - ts.indmom) .^ 2, 3));
    days = unique(T.Date);
    H = length(days);
    W = 3;
    tmpl = array2table(nan([H, W]));
    tmpl.Properties.VariableNames = {'Date', 'mom1m', 'indmom'};
    tmpl.Date = days;
    Mi = table(T.Date, mom1m, indmom, 'VariableNames', {'Date', 'mom1m', 'indmom'});
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = Mi(loc, :);
        pos = ismember(tmpl.Date, ts.Date);
        tmpl(pos, :) = ts;
        res = func(tmpl);
        idiovol(loc) = res(pos, :);
    end;
    idiovol = idiovol / 5.0320;
    
    vol3m = nan(shape);
    func = @(ts) log(erplib.shift(erplib.ma((ts.RETMONTH - ts.retmonth_spx) .^ 2, 3), 1));
    tmpl = template;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date);
        tmpl(pos, :) = ts;
        res = func(tmpl);
        vol3m(loc) = res(pos, :);
    end;
    vol3m = vol3m / 5.7736;
    
    vol6m = nan(shape);
    func = @(ts) log(erplib.shift(erplib.ma((ts.RETMONTH - ts.retmonth_spx) .^ 2, 6), 1));
    tmpl = template;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date);
        tmpl(pos, :) = ts;
        res = func(tmpl);
        vol6m(loc) = res(pos, :);
    end;
    vol6m = vol6m / 4.9238;
    
    indvol = nan(shape);
    days = unique(T.Date);
    sectors = unique(T.Sector);
    H = length(days);
    W = 1 + length(sectors);
    S = array2table(nan([H, W]));
    S.Properties.VariableNames = {'Date', sectors{:}};
    S.Date = days;
    M = table(T.Date, idiovol, 'VariableNames', {'Date', 'idiovol'});
    for j = 1:length(sectors) % cross-sectional
        block = M(strcmp(T.Sector, sectors(j)), :);
        S(:, sectors(j)) = ...
            array2table(arrayfun(@(day) nanmean(block(block.Date == day, :).idiovol), days));
    end;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(S.Date, ts.Date);
        res = S(:, ts.Sector(1));
        indvol(loc) = table2array(res(pos, :));
    end;
    indvol = indvol / 0.6224;
    
    
    %% size
    sz = log(T.marketcap);
    sz = sz / 8.4234;
    
    sznl = log(T.marketcap) .^ 3;
    sznl = sznl / 10844.7730;
    
    indsz = nan(shape);
    days = unique(T.Date);
    sectors = unique(T.Sector);
    H = length(days);
    W = 1 + length(sectors);
    S = array2table(nan([H, W]));
    S.Properties.VariableNames = {'Date', sectors{:}};
    S.Date = days;
    M = table(T.Date, sz, 'VariableNames', {'Date', 'sz'});
    for j = 1:length(sectors) % cross-sectional
        block = M(strcmp(T.Sector, sectors(j)), :);
        S(:, sectors(j)) = ...
            array2table(arrayfun(@(day) nanmean(block(block.Date == day, :).sz), days));
    end;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(S.Date, ts.Date);
        res = S(:, ts.Sector(1));
        indsz(loc) = table2array(res(pos, :));
    end;
    indsz = indsz / 0.1205;
    
    
    %% trading activity
    tvol1m = nan(shape);
    func = @(ts) erplib.shift(log(ts.Adj_Volume), 1);
    tmpl = template;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date);
        tmpl(pos, :) = ts;
        res = func(tmpl);
        tvol1m(loc) = res(pos, :);
    end;
    tvol1m(~isfinite(tvol1m)) = nan;
    tvol1m = tvol1m / 8.9045;
    
    tvol3m = nan(shape);
    func = @(ts) erplib.shift(log(erplib.ma(ts.Adj_Volume, 3)), 1);
    tmpl = template;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date);
        tmpl(pos, :) = ts;
        res = func(tmpl);
        tvol3m(loc) = res(pos, :);
    end;
    tvol3m(~isfinite(tvol3m)) = nan;
    tvol3m = tvol3m / 8.2009;
    
    turn1m = nan(shape);
    func = @(ts) log(erplib.shift(ts.Adj_Volume, 1) ./ ts.shareswa);
    tmpl = template;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date);
        tmpl(pos, :) = ts;
        res = func(tmpl);
        turn1m(loc) = res(pos, :);
    end;
    turn1m(~isfinite(turn1m)) = nan;
    turn1m = turn1m / 5.2891;
    
    turn6m = nan(shape);
    func = @(ts) log(erplib.ma(erplib.shift(ts.Adj_Volume, 1), 6) ...
                ./ erplib.ma(ts.shareswa, 6));
    tmpl = template;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date);
        tmpl(pos, :) = ts;
        res = func(tmpl);
        turn6m(loc) = res(pos, :);
    end;
    turn6m(~isfinite(turn6m)) = nan;
    turn6m = turn6m / 4.0797;
    
    sdturn = nan(shape);
    func = @(ts) erplib.std(ts.turn1m, 12); % TODO: chi-sqr -> norm
    days = unique(T.Date);
    H = length(days);
    W = 2;
    tmpl = array2table(nan([H, W]));
    tmpl.Properties.VariableNames = {'Date', 'turn1m'};
    tmpl.Date = days;
    N = table(T.Date, turn1m, 'VariableNames', {'Date', 'turn1m'});
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = N(loc, :);
        pos = ismember(tmpl.Date, ts.Date);
        tmpl(pos, :) = ts;
        res = func(tmpl);
        sdturn(loc) = res(pos, :);
    end;
    sdturn = sdturn / 0.3314;
    
    
    %% earning
    etp = nan(shape);
    func = @(ts) erplib.ma(ts.pe .^ -1, 12);
    tmpl = template;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date); % aligning missing dates
        tmpl(pos, :) = ts;
        res = func(tmpl);
        etp(loc) = res(pos, :);
    end;
    etp = etp / 1.1591;
    
    stp = nan(shape);
    func = @(ts) erplib.ma(ts.ps .^ -1, 12);
    tmpl = template;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date);
        tmpl(pos, :) = ts;
        res = func(tmpl);
        stp(loc) = res(pos, :);
    end;
    stp = stp / 5.7266;
    
    indetp = nan(shape);
    days = unique(T.Date);
    sectors = unique(T.Sector);
    H = length(days);
    W = 1 + length(sectors);
    S = array2table(nan([H, W]));
    S.Properties.VariableNames = {'Date', sectors{:}};
    S.Date = days;
    M = table(T.Date, etp, 'VariableNames', {'Date', 'etp'});
    for j = 1:length(sectors) % cross-sectional
        block = M(strcmp(T.Sector, sectors(j)), :);
        S(:, sectors(j)) = ...
            array2table(arrayfun(@(day) nanmean(block(block.Date == day, :).etp), days));
    end;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(S.Date, ts.Date);
        res = S(:, ts.Sector(1));
        indetp(loc) = table2array(res(pos, :));
    end;
    indetp = indetp / 0.3587;
    
    
    %% growth
    agr = nan(shape);
    func = @(ts) (erplib.ma(ts.assets, 12) -  erplib.ma(ts.assets, 24)) ...
                ./ erplib.ma(ts.assets, 24);
%     func = @(ts) log(1 + ...
%         (erplib.ma(ts.assets, 12) -  erplib.ma(ts.assets, 24)) ./ erplib.ma(ts.assets, 24));
    tmpl = template;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date); % aligning missing dates
        tmpl(pos, :) = ts;
        res = func(tmpl);
        agr(loc) = res(pos, :);
    end;
	arg = agr / 10.9546;
    
    egr12m = nan(shape);
    func = @(ts) (erplib.ma(ts.epsusd, 6) -  erplib.ma(ts.epsusd, 12)) ...
                ./ erplib.ma(ts.epsusd, 12);
    tmpl = template;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date); % aligning missing dates
        tmpl(pos, :) = ts;
        res = func(tmpl);
        egr12m(loc) = res(pos, :);
    end;
    egr12m(egr12m < -10) = -10;
    egr12m(egr12m > 10) = 10;
    egr12m = egr12m / 5.4837;
    
    egr24m = nan(shape);
    func = @(ts) (erplib.ma(ts.epsusd, 12) -  erplib.ma(ts.epsusd, 24)) ...
                ./ erplib.ma(ts.epsusd, 24);
    tmpl = template;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date); % aligning missing dates
        tmpl(pos, :) = ts;
        res = func(tmpl);
        egr24m(loc) = res(pos, :);
    end;
    egr24m(egr24m < -10) = -10;
    egr24m(egr24m > 10) = 10;
    egr24m = egr24m / 4.7879;
    
    rdtsz = nan(shape);
    func = @(ts) erplib.ma(ts.rnd ./ ts.marketcap, 12); % TODO: better shape after log
    tmpl = template;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date); % aligning missing dates
        tmpl(pos, :) = ts;
        res = func(tmpl);
        rdtsz(loc) = res(pos, :);
    end;
    rdtsz = rdtsz / 0.07632;
    
    indagr = nan(shape);
    days = unique(T.Date);
    sectors = unique(T.Sector);
    H = length(days);
    W = 1 + length(sectors);
    S = array2table(nan([H, W]));
    S.Properties.VariableNames = {'Date', sectors{:}};
    S.Date = days;
    M = table(T.Date, agr, 'VariableNames', {'Date', 'agr'});
    for j = 1:length(sectors) % cross-sectional
        block = M(strcmp(T.Sector, sectors(j)), :);
        S(:, sectors(j)) = ...
            array2table(arrayfun(@(day) nanmean(block(block.Date == day, :).agr), days));
    end;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(S.Date, ts.Date);
        res = S(:, ts.Sector(1));
        indagr(loc) = table2array(res(pos, :));
    end;
    indagr = indagr / 0.0941;
    
    
    %% earning variability
    envr = nan(shape);
    func = @(ts) erplib.std(ts.epsusd, 12) ./ erplib.ma(ts.epsusd, 12);
    tmpl = template;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date); % aligning missing dates
        tmpl(pos, :) = ts;
        res = func(tmpl);
        envr(loc) = res(pos, :);
    end;
    envr = real(envr);
    envr(envr < -30) = -30;
    envr(envr > 30) = 30;
    envr = envr / 12.1284;
    
    
    %% dividend
    dy = nan(shape);
    func = @(ts) erplib.ma(ts.divyield, 12);
    tmpl = template;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(tmpl.Date, ts.Date); % aligning missing dates
        tmpl(pos, :) = ts;
        res = func(tmpl);
        dy(loc) = res(pos, :);
    end;
    dy = dy / 0.0630;
    
    %% leverage
    dte = T.de;
    dte(abs(dte) > 100) = 100;
    dte = dte / 18.739;
    
    inddte = nan(shape);
    days = unique(T.Date);
    sectors = unique(T.Sector);
    H = length(days);
    W = 1 + length(sectors);
    S = array2table(nan([H, W]));
    S.Properties.VariableNames = {'Date', sectors{:}};
    S.Date = days;
    M = table(T.Date, dte, 'VariableNames', {'Date', 'dte'});
    for j = 1:length(sectors) % cross-sectional
        block = M(strcmp(T.Sector, sectors(j)), :);
        S(:, sectors(j)) = ...
            array2table(arrayfun(@(day) nanmean(block(block.Date == day, :).dte), days));
    end;
    for j = 1:length(cids)
        loc = logical(T.compid == cids(j));
        ts = T(loc, :);
        pos = ismember(S.Date, ts.Date);
        res = S(:, ts.Sector(1));
        inddte(loc) = table2array(res(pos, :));
    end;
    inddte = inddte / 0.1424;
    
    % sentimental
    
    
    %% aggregation
    signals = table( ...
        mom1m, mom3m, mom6m, mom12m, chmom, indmom, ...
        beta, ...
        gvol, idiovol, vol3m, vol6m, indvol, ...
        sz, sznl, indsz, ...
        tvol1m, tvol3m, turn1m, turn6m, sdturn, ...
        etp, stp, indetp, ...
        agr, egr12m, egr24m, rdtsz, indagr, ...
        envr, ...
        dy, ...
        dte, inddte);
    
    
    % ex-return
    exrtn = log(1 + T.RETMONTH) - log(1 + T.retmonth_spx);
    
    save('data.mat', 'T', 'signals', 'exrtn');
    
    %% range check
    names = signals.Properties.VariableNames;
    fprintf('------------------------------------------------------------\n');
    for j = 1:length(names)
        ts = table2array(signals(:, names{j}));
        q = quantile(ts, [0.025 0.975]);
        fprintf('%12s: [%7.3f, %7.3f] [%8.4f] (%4.2f%%)\n', ...
            names{j}, q(1), q(2), q(2) - q(1), 100 * sum(isfinite(ts)) / length(ts));
    end;
    
    
% end

