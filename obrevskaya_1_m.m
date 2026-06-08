clear all
close all
clc
pkg load statistics

##absciss = -pi:0.01:pi;
##plot(absciss, sin(absciss))
data = load('kyrs1_11.txt');
plot(data)

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

for sample = 1:samplesNum
  curData = data{sample}(:,1);
  M = mean(curData);
  S = std(curData);
  n = length(curData);
  skewness = sum((curData - M).^3) / (n * S^3);
  absSkew = abs(skewness);
  DA = 6*(n - 1)/((n + 1) *(n + 3));
  critSkew = 3 * sqrt(DA);

  kurtosis = sum((curData - M).^4) / (n * S^4) - 3;
  absKurt = abs(kurtosis);
  DE = 24*n*(n - 2)*(n-3)/((n+1)^2*(n + 3) *(n + 5));
  critKurt = 5 * sqrt(DE);

  fprintf('Sample #%d:\n', sample);
  fprintf('  Skewness = %.5f, its std = %.5f)\n', skewness, DA);
  fprintf('  |A| = %.5f, 3 sqrt(Da) = %.5f)\n', absSkew, critSkew);
  if absSkew > critSkew
      fprintf('  |A| > 3 sqrt(Da), hypothesis rejected\n');
  else
      fprintf('  |A| <= 3 sqrt(Da), hypothesis accepted\n');
  end

  fprintf('  Kurtosis = %.5f, its std = %.5f\n', kurtosis, DE);
  fprintf('  |E| = %.5f, 5 sqrt(De) = %.5f)\n', absKurt, critKurt);
  if absKurt > critKurt
      fprintf('  |E| > 5 sqrt(De), hypothesis rejected\n');
  else
      fprintf('  |E| <= 5 sqrt(De), hypothesis accepted\n');
  end
  fprintf('\n');
end

## punct 3.2
filteredData = data;
data = [];
for sample =1:samplesNum
  data = [data;filteredData{sample}(:,1)];
endfor
n = length(data)
k = round(log(n)*1.44+1)##kol-vo intervalov
minVal = min (data);
delta = (max(data)-minVal)/k;
Nm = zeros(1,k);
Pm = zeros(1,k);

M = mean(data);
S = std(data);
for i = 1:k
  leftBoarder = minVal + (i-1) * delta;
  rightBoarder = minVal + i * delta;
  if i == 1
    leftBoarder = -inf;
  endif
  if i == k
    rightBoarder = inf;
  endif
  Nm(i) = nnz(data <= rightBoarder & data > leftBoarder);
  Pm(i) = normcdf((rightBoarder - M)/S) - normcdf((leftBoarder - M)/S);
endfor

assert(sum(Nm)==n)
assert(abs(sum(Pm)-1) < eps)

figure
hold on
bar(Nm)
Em = Pm * n;## ?theoriticheskie chastoty
plot(Em)
xlabel('k,номер интервала')
ylabel('p, плотность вероятности')
legend('р расч','p теор')

Chi2 = sum(((Nm - Em).^2)./Em)
fprintf('квантиль Пирсона: %f\n', chi2inv(1-0.05, k-3));
if Chi2 > chi2inv(1-0.05, k-3)
  fprintf('Chi2 > chi2inv, hypothesis rejected\n');
else
  fprintf('Chi2 <= chi2inv, hypothesis accepted\n');
end

####08/10
data = filteredData{1}(:,1);
data = sort(data);
n = length(data);
Fn = (1:n)/n;
Fn = Fn';## transponirovali vector
F = normcdf((data - mean(data))/std(data));

figure
hold on
stairs(data, Fn)
plot(data, F)

[maxDelta, maxDeltaIdx] = max(abs(F - Fn));
[maxDelta2, maxDeltaIdx2] = max(abs(F(2:end) -Fn(1:(end-1))));
if maxDelta2 > maxDelta
  maxDelta = maxDelta2;
  maxDeltaIdx = maxDeltaIdx2;
endif
clear maxDelta2, maxDeltaIdx2
maxDelta
lambda = maxDelta * sqrt(n)

plot(repmat(data(maxDeltaIdx),1,2), [0,1])
xlabel('случайная величина')
ylabel('F, функция распределения СВ')
legend('F эмпир','F теор')

fprintf('квантиль Колмогорова: %f\n', kolminv(1-0.1));
if lambda > kolminv(1-0.1)
  fprintf('lambda > lambda_(1-q), hypothesis rejected\n');
else
  fprintf('lambda <= lambda_(1-q), hypothesis accepted\n');
end

####15/10 punct 4

data = [];
for sample =1:samplesNum
  data = [data;filteredData{sample}(:,1)];
endfor
n = length(data)
n = round(n/2);
r = zeros(1,n);
for i = 1:n
  r(i) = Rx(data, i-1);
endfor

figure
plot(r)
xlabel('интервал времени')
ylabel('R, корреляционная функция')

##22/10
##punct 5 1 abzac
minLen = length(filteredData{1});
for sample  = 2:samplesNum
  minLen = min(minLen, length(filteredData{sample}));
endfor
minLen

data = zeros(minLen, samplesNum);
for sample = 1:samplesNum
  data(:,sample) = filteredData{sample}(1:minLen,1);
endfor
D = var(data);
g = max(D)/sum(D)
g_kohran = 0.2434
fprintf('квантиль Кохрана\n');
if g > g_kohran
  fprintf('g > g_(1-q), hypothesis rejected\n');
else
  fprintf('g <= g_(1-q), hypothesis accepted\n');
end
##https://rtlam.blogspot.com/p/tabellen.html
##0.2434 g bolshe, process ne statsionaren

M = zeros(1,samplesNum);
D = zeros(1,samplesNum);
n = 0;
for sample = 1:samplesNum
  curData = filteredData{sample}(:,1);
  M(sample) = mean(curData);
  D(sample) = var(curData);
  n = n + length(curData);
end
clear curData

## punct 6
f = (var(M)/mean(D)) * n
fprintf('квантиль Фишера: %f\n', finv(1-0.05, samplesNum-1, samplesNum * (n - 1)));

if f > finv(1-0.05, samplesNum-1, samplesNum * (n - 1))
  fprintf('F > F(1-q), hypothesis rejected\n');
else
  fprintf('F <= F(1-q), hypothesis accepted\n');
end


rx = abs(r/mean(D));
rq = 0.195;
figure
plot(rx)
hold on
plot(xlim,[rq, rq])
legend('r_x', 'r_{1-q} = 0.195')

















