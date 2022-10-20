function rsquared = multiRegress3var(DV, IV1, IV2, IV3)

DV = tiedrank(DV,0);
IV1 = tiedrank(IV1,0);
IV2 = tiedrank(IV2,0);
IV3 = tiedrank(IV3,0);


tableData = table;
tableData.a = DV;
    
tableData.b = IV1;
tableData.c = IV2;
tableData.d = IV3;

    
analysisString = 'a~b+c+d';

linearModel = fitlm(tableData,analysisString);                                                      

rsquared = linearModel.Rsquared.Ordinary;