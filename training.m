
[X, loc] = rmmissing(signals);
X = X.Variables;
y = exrtn(~loc, :);

loc = isnan(y);
X = X(~loc, :);
y = y(~loc, :);

loc = abs(log(1+y) - mean(log(1+y))) < 6 * std(log(1+y));
loc = loc & abs(log(1+y) - mean(log(1+y))) > .25 * std(log(1+y));
X = X(loc, :);
y = y(loc, :);

net = feedforwardnet([10 4 2]);
net = train(net, X', y');

% view(net)