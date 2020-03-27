function basicStatsStruct = calcBasicStatistics(data,basicStatsStruct)


if(size(data,1) > 0)
  
  dataSort = sort(data);
  n = length(data);
  n05      = max(round(n*0.05),1);
  n25      = max(round(n*0.25),1);
  n75      = max(round(n*0.75),1);
  n95      = max(round(n*0.95),1);

  basicStatsStruct.min    = min(data);
  basicStatsStruct.mean   = mean(data);
  basicStatsStruct.median = median(data);
  basicStatsStruct.max    = max(data);
  basicStatsStruct.n = length(data);
  
  basicStatsStruct.p05 = dataSort(n05,1);
  basicStatsStruct.p25 = dataSort(n25,1);
  basicStatsStruct.p75 = dataSort(n75,1);
  basicStatsStruct.p95 = dataSort(n95,1);  
end
