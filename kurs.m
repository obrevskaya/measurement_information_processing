clear all
close all
clc

##absciss = -pi:0.01:pi;
##plot(absciss, sin(absciss))
data = load('kyrs1_2.txt');
##plot(data)

origData = data;

q = 0.05;


data = num2cell(data, 1);
observNum = length(data{1});
samplesNum = length(data);
outliers = {};

for sample = 1:samplesNum
  data{sample} = [data{sample}, (1:observNum)'];
  outliers{end+1} = [];

  while 1
    sampleVals = data{sample};
    M = mean(sampleVals (:, 1));
    S = std(sampleVals (:, 1));

    [val, idx] = min(sampleVals (:, 1));
    P = 1 - erf(((M - val) / S) / sqrt(2));

    assert(P<1)
    assert(P>0)
    allDone = 1;

    if P < q
      absIdx = sampleVals(idx,2);
      allDone = 0;
      sampleVals(idx, :) = [];
      outliers{sample} = [outliers{sample}, absIdx];

    endif

    [val, idx] = max(sampleVals (:,1));
    P = 1 - erf(((val - M) / S) / sqrt(2));

    assert(P<1)
    assert(P>0)

    if P < q
      absIdx = sampleVals(idx,2);
      allDone = 0;
      sampleVals(idx, :) = [];
      outliers{sample} = [outliers{sample}, absIdx];

    endif

    if allDone
      break
    endif

    data{sample} = sampleVals;
  end
end

figure
hold on

for sample = 1:samplesNum
  plot (origData(:, sample))
  curoutliers = outliers{sample};
  plot(curoutliers, origData(curoutliers, sample), 'ro')
end















