function rsquared = multiRegress2var(DV, IV1, IV2)

DV = tiedrank(DV,0);
IV1 = tiedrank(IV1,0);
IV2 = tiedrank(IV2,0);


tableData = table;
tableData.a = DV;
    
tableData.b = IV1;
tableData.c = IV2;

    
analysisString = 'a~ b + c';

linearModel = fitlm(tableData,analysisString);                                                      

rsquared = linearModel.Rsquared.Ordinary;