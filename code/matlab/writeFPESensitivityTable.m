function success = writeFPESensitivityTable(tableFolderName,...
                                        groupData,groupMetaData,...
                                        metricNameList,metricSubFields,...
                                        trialsToProcess,trialTypeNames,...
                                        phaseNames,modificationName)
                                
success = 0;




for indexTrialsToProcess=1:1:size(groupData,2)

  indexTrial = 0;
  for k=1:1:length(trialTypeNames)
    if(contains(trialsToProcess{indexTrialsToProcess},trialTypeNames{k}))
      indexTrial = k;
    end
  end

  for indexPhase=1:1:size(groupData,3)

    csvData = zeros(length(groupMetaData)*length(metricSubFields) , length(metricNameList));

    
    csvDataLabels = cell(length(groupMetaData)*length(metricSubFields),3);
    

    indexRow = 1;
    for indexGroup=1:1:length(groupMetaData)      
      for indexMetric=1:1:length(metricNameList)
          metricName = metricNameList{indexMetric};

          for indexField=1:1:length(metricSubFields)

              fieldName = metricSubFields{indexField};

              indexRow = indexGroup + indexField - 1;

              if( isfield( groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName),fieldName)==1)            
                if(isfield( groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName),'median') == 1)
                  if(isempty( groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).median) == 0)
                    minData    = groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).min;
                    meanData   = groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).mean;
                    maxData    = groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).max;
                    medianData = groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).median;

                    p05Data = groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).p05;
                    p25Data = groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).p25;
                    p75Data = groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).p75;
                    p95Data = groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).p95;

                    csvData(indexRow,indexMetric) = maxData;
                    csvDataLabels{indexRow,1} = groupData(indexGroup,indexTrialsToProcess,indexPhase).name;
                    csvDataLabels{indexRow,2} = fieldName;
                    csvDataLabels{indexRow,3} = 'max';
                    

                    

                  end
                end
              end
             
          end
      end

    end

    tableName = [tableFolderName,'tableFPESensitivity',...
      trialTypeNames{indexTrial},...
      phaseNames{indexPhase},'Group',...
      modificationName,'.csv'];

    fid =fopen(tableName,'w');

    fprintf(fid,',,');
    for i=1:1:size(metricNameList,2)
      fprintf(fid,',%s',metricNameList{1,i});      
    end
    fprintf(fid,',\n');
    

    emptyLine = ',,,';
    for i=1:1:size(csvData,1)
      emptyLine = [emptyLine,','];
    end
    emptyLine = [emptyLine,'\n'];



    for i=1:1:size(csvData,1)

      fprintf(fid,'%s', csvDataLabels{i,1});
      fprintf(fid,',%s',csvDataLabels{i,2});
      fprintf(fid,',%s',csvDataLabels{i,3});


      for j=1:1:size(csvData,2)
        fprintf(fid,',%1.6f',csvData(i,j));
      end
      fprintf(fid,',\n');   
    end



    fclose(fid);
  end
end
  
  success = 1;                                