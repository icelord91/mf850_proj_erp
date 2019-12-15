
classdef erplib
    
    % TECHNICAL ANALYSIS BASIC FORMULA
    methods (Static)

        % Shift series
        function ts_ = shift(ts, n)
            
            ts_ = nan(size(ts));

            % y_(max(1, n+1):min(numel(ts), numel(ts)+n)) ...
            if n == 0
                ts_ = ts;
            elseif n > 0
                ts_(n+1:end) = ts(1:end-n);
            else
                ts_(1:numel(ts)+n) = ts(-n+1:end);
            end;

        end
        
        % Arithmetic average
        function mu = ma(ts, n)

            % base mean calculation
            mu = filter(ones(1,n)/n, 1, ts);

            % head correction
            mu(1) = ts(1);
            for i = 2:min(n-1, numel(ts))
                mu(i) = mu(i-1) * (i-1) / i + ts(i) / i;
            end;

        end

        % Standard deviation
        function sigma = std(ts, n)

            % base mean calculation (derived form mean calculation)
            sigma = sqrt(filter(ones(1,n)/n, 1, ts.^2) - filter(ones(1,n)/n, 1, ts).^2);

            % head correction
            sigma(1) = 0;
            for i = 2:min(n-1, numel(ts))
                sigma(i) = std(ts(1:i), 1);
            end;

        end

        % Exponential moving average
        function ema_ = ema(ts, n)

            len = length(ts);
            n = min(n, len);
            alpha = 2 / (n + 1);
            beta = 1 - alpha;
            % geometric upper triangle transformation matrix
            tr = alpha * diag(beta .^ -[1:len]) * triu(ones(len)) * diag(beta .^ [1:len]);
            % first term (number) correction
            tr(:,len) = tr(:,len) + beta .^ [len:-1:1]';
            ema_ = flip(tr * flip(ts));

        end
        
        % Trim
        function x_ = trim(x, varargin)
            
            l = -1;
            r = 2;
            if nargin > 2
                l = range(1);
                r = range(2);
            elseif nargin > 1
                r = range(1);
            end;
            
            x_ = nan(size(x));
            mu = nanmean(x);
            sig = nanstd(x);
            loc = (abs(x - mu) > l * sig) & (abs(x - mu) < r * sig);
            x_(loc, :) = x(loc, :);
            
        end
        
        % Normalization
        function x_ = normalize(x)
            
            mu = nanmean(x);
            sig = nanstd(x);
            x_ = (x - mu) ./ sig;

        end
        
    end
    
end