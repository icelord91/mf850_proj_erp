%fulldata = readtable("fulldata.csv");

%returns = fulldata.RETMONTH_end;

%colnames = strtrim(char(fulldata.Properties.VariableNames));

%inputs = removevars(fulldata,{'Date', 'compid', 'Industry', 'RETMONTH_end'});
%inputs = fulldata(:,[2 3 5 14 19 20 65 86 88]);

%direction = sign(returns);

layers = [sequenceInputLayer(9),
          tanhLayer,
          fullyConnectedLayer(6),
          tanhLayer,
          fullyConnectedLayer(4),
          tanhLayer,
          fullyConnectedLayer(2),
          fullyConnectedLayer(1),
          regressionLayer];
      
options = trainingOptions('sgdm', ...
                          'plots', 'training-progress', ...
                          'MaxEpochs',300, ...
                          'InitialLearnRate',0.0000001, ...
                          'MiniBatchSize',256, ...
                          'Shuffle','every-epoch', ...
                          );

input_array = normalize(table2array(inputs).', 2);

net = trainNetwork(input_array, returns.', layers, options);

guess = predict(net, input_array);




sent=inputs(:,{'sentiment_bullish_end', 'sentiment_neutral_end', 'sentiment_bearish_end'});
input_sent = table2array(sent).';

sent_layers=[sequenceInputLayer(3),
             tanhLayer,
             fullyConnectedLayer(5),
             tanhLayer,
             fullyConnectedLayer(3),
             tanhLayer,
             fullyConnectedLayer(2),
             tanhLayer,
             fullyConnectedLayer(1),
             regressionLayer];

sent_net = trainNetwork(input_sent, returns.', sent_layers, options);