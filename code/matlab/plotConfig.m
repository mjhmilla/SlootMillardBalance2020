%%
% Setup plot parameters
%%
%totalWidth = 177.13668/10; %Frontiers journal text width.
totalWidth = 29.7;

%Landscape
pageWidth  = 29.7;
pageHeight = 21.0;


numberOfVerticalPlotRows = 2;
numberOfHorizontalPlotColumns = 3;

%if(flag_usingOctave == 0)
set(groot, 'defaultAxesFontSize',8);
set(groot, 'defaultTextFontSize',8);
set(groot, 'defaultAxesLabelFontSizeMultiplier',1.2);
set(groot, 'defaultAxesTitleFontSizeMultiplier',1.2);
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTitleFontWeight','bold');  
set(groot, 'defaultFigurePaperUnits','centimeters');
set(groot, 'defaultFigurePaperSize',[pageWidth pageHeight]);
set(groot,'defaultFigurePaperType','A4');
%end
plotHorizMarginCm = 2;
plotVertMarginCm  = 2;

plotHeight= ((pageHeight-plotVertMarginCm)/numberOfVerticalPlotRows);

plotWidth = ((pageWidth-plotHorizMarginCm)/numberOfHorizontalPlotColumns);
plotWidth  = plotWidth/pageWidth;
plotHeight = plotHeight/pageHeight;

plotHorizMargin = plotHorizMarginCm/pageWidth;
plotVertMargin  = plotVertMarginCm/pageHeight;

topLeft = [0/pageWidth pageHeight/pageHeight];

subPlotPanel=zeros(numberOfVerticalPlotRows,numberOfHorizontalPlotColumns,4);
subPlotPanelIndex = zeros(numberOfVerticalPlotRows,numberOfHorizontalPlotColumns);

idx=1;
scaleVerticalMargin = 0.;
for(ai=1:1:numberOfVerticalPlotRows)
  if(ai > 1)
    scaleVerticalMargin = 1;
  end
  for(aj=1:1:numberOfHorizontalPlotColumns)
      subPlotPanelIndex(ai,aj) = idx;
      subPlotPanel(ai,aj,1) = topLeft(1) + plotHorizMargin...
                            + (aj-1)*(plotWidth);
      %-plotVertMargin*scaleVerticalMargin ...                             
      subPlotPanel(ai,aj,2) = topLeft(2) -plotHeight ...                            
                            + (ai-1)*(-plotHeight);
      subPlotPanel(ai,aj,3) = (plotWidth-plotHorizMargin);
      subPlotPanel(ai,aj,4) = (plotHeight-plotHorizMargin);
      idx=idx+1;
  end
end




