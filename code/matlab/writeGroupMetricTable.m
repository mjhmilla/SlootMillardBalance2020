function success = writeGroupMetricTable(tableFolderName,...
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

      csvData = zeros(9 ,length(metricNameList)*length(groupMetaData)*length(metricSubFields));
      csvRowLabels = {'median','p25p75','min','mean','max','p05','p25','p75','p95'};
      csvColumnLabelsA = cell(1,length(metricNameList)*length(groupMetaData)*length(metricSubFields));
      csvColumnLabelsB = cell(1,length(metricNameList)*length(groupMetaData)*length(metricSubFields));
      csvColumnLabelsC = cell(1,length(metricNameList)*length(groupMetaData)*length(metricSubFields));

      idxColumn = 1;
      for indexMetric=1:1:length(metricNameList)
        metricName = metricNameList{indexMetric};

          idxRow = 1;
          for indexField=1:1:length(metricSubFields)
            for indexGroup=1:1:length(groupMetaData)

              fieldName = metricSubFields{indexField};

              csvColumnLabelsA{1,idxColumn} = metricName;
              csvColumnLabelsB{1,idxColumn} = fieldName;            
              csvColumnLabelsC{1,idxColumn} = groupMetaData(indexGroup).name;

              if( isfield( groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName),fieldName)==1)            
                if(isfield( groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName),'median') == 1)
                  if(isempty( groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).median) == 0)
                    minData    = groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).min;
                    meanData  = groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).mean;
                    maxData    = groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).max;
                    medianData = groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).median;

                    p05Data = groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).p05;
                    p25Data = groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).p25;
                    p75Data = groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).p75;
                    p95Data = groupData(indexGroup,indexTrialsToProcess,indexPhase).(metricName).(fieldName).p95;

                    csvData(1,idxColumn) = medianData;
                    csvData(2,idxColumn) = p75Data-p25Data;
                    csvData(3,idxColumn) = minData;
                    csvData(4,idxColumn) = meanData;
                    csvData(5,idxColumn) = maxData;
                    csvData(6,idxColumn) = p05Data;
                    csvData(7,idxColumn) = p25Data;
                    csvData(8,idxColumn) = p75Data;
                    csvData(9,idxColumn) = p95Data;              
                  end
                end
              end
              idxColumn=idxColumn+1;          
            end
          end
      end
      idxColumn = idxColumn-1;

      

      tableName = [tableFolderName,'table',...
        trialTypeNames{indexTrial},...
        phaseNames{indexPhase},'Group',...
        modificationName,'.csv'];

      fid =fopen(tableName,'w');

      fprintf(fid,',,');
      for i=1:1:size(csvData,1)
        fprintf(fid,',%s',csvRowLabels{1,i});      
      end
      fprintf(fid,',\n');
      
      
      emptyLine = ',,,';
      for i=1:1:size(csvData,1)
        emptyLine = [emptyLine,','];
      end
      emptyLine = [emptyLine,'\n'];
      
      lastLabel = csvColumnLabelsA{1,1};
      
      for j=1:1:size(csvData,2)
        if( j>1)
          if(strcmp(lastLabel,csvColumnLabelsA{1,j})==0)
            fprintf(fid,emptyLine);
          end
        end
        
        fprintf(fid,'%s',csvColumnLabelsA{1,j});
        fprintf(fid,',%s',csvColumnLabelsB{1,j});
        fprintf(fid,',%s',csvColumnLabelsC{1,j});
        
        lastLabel = csvColumnLabelsA{1,j};
        
        for i=1:1:size(csvData,1)
          fprintf(fid,',%1.6f',csvData(i,j));
        end
        fprintf(fid,',\n');   
      end
            


      fclose(fid);
    end
  end
  
  success = 1;
