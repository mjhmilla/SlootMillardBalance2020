function success = writeParticipantMetricTable(tableFolderName,...
                                        subjectData,...
                                        metricNameList,metricSubFields,...
                                        trialsToProcess,trialTypeNames,...
                                        phaseNames,modificationName)
success = 0;
for indexTrialsToProcess=1:1:size(subjectData,2)
    
  indexTrial = 0;
  for k=1:1:length(trialTypeNames)
    if(contains(trialsToProcess{indexTrialsToProcess},trialTypeNames{k}))
      indexTrial = k;
    end
  end

  for indexPhase=1:1:size(subjectData,3)

    tableName = [tableFolderName,'table',...
      trialTypeNames{indexTrial},...
      phaseNames{indexPhase},'Part',...
      modificationName,'.csv'];

    fid =fopen(tableName,'w');

    for indexMetric=1:1:length(metricNameList)
      metricName = metricNameList{indexMetric};

      for indexField=1:1:length(metricSubFields)
        numberOfSubjects = size(subjectData,1);

        fieldName = metricSubFields{indexField};

        csvData = zeros(4, numberOfSubjects);
        csvRowLabels = {'median','p25p75','p25','p75','n'};
        csvColLabels = cell(1,numberOfSubjects+2);
        csvColLabels{1,1} = metricName;
        csvColLabels{1,2} = fieldName;
        subjectIndex = 1;

        for indexSubject = 1:1:size(subjectData,1)
          csvColLabels{1,2+indexSubject}= ...
            subjectData(indexSubject,indexTrialsToProcess,indexPhase).('subjectId');


          if( isfield( subjectData(indexSubject,indexTrialsToProcess,indexPhase).(metricName),fieldName)==1)            
            if(isfield( subjectData(indexSubject,indexTrialsToProcess,indexPhase).(metricName).(fieldName),'median') == 1)
              if(isempty( subjectData(indexSubject,indexTrialsToProcess,indexPhase).(metricName).(fieldName).median) == 0)

                minData    = subjectData(indexSubject,indexTrialsToProcess,indexPhase).(metricName).(fieldName).min;
                meanData   = subjectData(indexSubject,indexTrialsToProcess,indexPhase).(metricName).(fieldName).mean;
                maxData    = subjectData(indexSubject,indexTrialsToProcess,indexPhase).(metricName).(fieldName).max;
                medianData = subjectData(indexSubject,indexTrialsToProcess,indexPhase).(metricName).(fieldName).median;
                p05Data = subjectData(indexSubject,indexTrialsToProcess,indexPhase).(metricName).(fieldName).p05;
                p25Data = subjectData(indexSubject,indexTrialsToProcess,indexPhase).(metricName).(fieldName).p25;
                p75Data = subjectData(indexSubject,indexTrialsToProcess,indexPhase).(metricName).(fieldName).p75;
                p95Data = subjectData(indexSubject,indexTrialsToProcess,indexPhase).(metricName).(fieldName).p95;
              
                nData = subjectData(indexSubject,indexTrialsToProcess,indexPhase).(metricName).(fieldName).n;
                
                csvData(1,indexSubject) = medianData;
                csvData(2,indexSubject) = p75Data-p25Data;
                csvData(3,indexSubject) = p25Data;
                csvData(4,indexSubject) = p75Data;
                csvData(5,indexSubject) = nData;
              end
            end
          end

        end
        
        emptyRow = ',';
        fprintf(fid,'%s',csvColLabels{1,1});
        for i=2:1:size(csvColLabels,2)
          fprintf(fid,',%s',csvColLabels{1,i});
          emptyRow = [emptyRow,','];
        end
        fprintf(fid,'\n');

        for i=1:1:size(csvData,1)
          fprintf(fid,',%s',csvRowLabels{1,i});
          fprintf(fid,',%1.6f',csvData(i,1));
          for j=2:1:size(csvData,2)
            fprintf(fid,',%1.6f',csvData(i,j));
          end
          fprintf(fid,'\n');   
        end
        fprintf(fid,'%s\n',emptyRow);        
      end




    end



    fclose(fid);
  end
end

  success = 1;
