
%% data loading
load('data');

X = table2array(signals);
loc = sum(~isfinite(X), 2) > 0;
X = X(~loc, :);
% [X, loc] = rmmissing(table2array(signals));
% [X, loc] = rmmissing(signals);
% X = X.Variables;
y = exrtn(~loc, :);

loc = isnan(y);
X = X(~loc, :);
y = y(~loc, :);

loc = abs(log(1+y) - mean(log(1+y))) < 6 * std(log(1+y));
loc = loc & abs(log(1+y) - mean(log(1+y))) > .25 * std(log(1+y));
X = X(loc, :);
y = y(loc, :);

%% neural network for regression
hiddenLayerSize = [32 16 4];
net = feedforwardnet(hiddenLayerSize);

net.trainFcn = 'trainbfg'; % 'trainlm' (default) | 'trainbfg' | 'trainbr' | 'traingdm' | 'traingdx'
net.performFcn = 'mse'; % 'mae' | 'sse' | 'sae'

net.performParam.normalization = 'none'; % 'none' (default) | 'standard' | 'percent'
net.performParam.regularization = 0.01;
net.trainParam.max_fail = 16;
net.trainParam.min_grad = 1e-7;

net.divideParam.trainRatio = 0.8;
net.divideParam.valRatio = 0.1;
net.divideParam.testRatio = 0.1;

net = train(net, X', y');

R = sqrt(1 - mean((y - net(X')') .^ 2) / var(y))

%% neural network for classification
hiddenLayerSize = [30 18 6];
net = patternnet(hiddenLayerSize);

net.trainFcn = 'trainlm'; % 'trainlm' (default) | 'trainbfg' | 'trainbr' | 'traingdm' | 'traingdx'
net.performFcn = 'mse'; % 'mae' | 'sse' | 'sae'

net.performParam.normalization = 'none'; % 'none' (default) | 'percent'
net.performParam.regularization = 0.01;
net.trainParam.max_fail = 16;

net.divideParam.trainRatio = 0.8;
net.divideParam.valRatio = 0.1;
net.divideParam.testRatio = 0.1;

y = (sign(y) + 1) / 2;
net = train(net, X', y');

% Acc = 

%% performance
% view(net)

yhat = net(X');
err = gsubtract(yhat, y');
performance = perform(net, y', yhat);

% R = sqrt(1 - mean((y - yhat') .^ 2) / var(y))