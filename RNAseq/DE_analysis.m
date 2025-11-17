%% Load data
% PC Filename
file = 'C:\Users\Tom\Documents\iCloudDrive\Documents\My Papers\p38\Data\RNAseq\gene_count2.xlsx';
cd 'C:\Users\Tom\Documents\iCloudDrive\Documents\My Papers\p38\Data\RNAseq';

% Mac Filename
%file = '/Users/thomaskuhlman/Documents/My Papers/p38/Data/RNAseq/gene_count2.xlsx';
%cd 'Users/thomaskuhlman/Documents/My Papers/p38/Data/RNAseq';


% Load data
opts = detectImportOptions(file);
opts = setvartype(opts, "gene_chr", 'string');
opts.VariableNamingRule = 'preserve';
geneCountTable = readtable(file, opts);

% Sample names and counts
counts = geneCountTable{:,2:32};
geneCountTable = sortrows(geneCountTable,'gene_id','ascend');
gene_names = geneCountTable{:,33};


% Estimate pseudo-reference with geometric mean row by row
pseudoRefSample = geomean(counts,2);
nz = pseudoRefSample > 0;
ratios = bsxfun(@rdivide, counts(nz,:), pseudoRefSample(nz));
sizeFactors = median(ratios,1);


% Transform to common scale
normCounts = bsxfun(@rdivide, counts, sizeFactors);


decutoff = 1;
sig_cutoff = 0.01;

%% C20mM Analysis

% Consider the mean
meanUntreated = mean([normCounts(:,1) normCounts(:,2)], 2);
meanTreated = mean([normCounts(:,3) normCounts(:,4)], 2);

% Consider the dispersion
dispUntreated = std([normCounts(:,1) normCounts(:,2)], 0, 2) ./meanUntreated;
dispTreated = std([normCounts(:,3) normCounts(:,4)], 0, 2) ./meanTreated;


% Plot on loglog scale
%{
figure
loglog(meanTreated, dispTreated, 'r.');
hold on;
loglog(meanUntreated, dispUntreated, 'b.');
xlabel('log2(Mean)');
ylabel('log2(Dispersion)');
legend('Treated','Untreated','Location','southwest')
%}


% Compute the mean and the log2FC
meanBase = (meanTreated + meanUntreated)/2;
foldChange = meanTreated ./ meanUntreated;
log2FC = log2(foldChange);


% Plot mean vs. fold change (MA plot)
%{
mairplot(normCounts(:,3), meanUntreated,'Type','MA','Plotonly',true);
set(get(gca,'Xlabel'),'String','mean of normalized counts')
set(get(gca,'Ylabel'),'String','log2(fold change)')
%}


% Create a table with statistics about each gene
geneTableC20 = table(gene_names, meanBase, meanTreated, meanUntreated, foldChange, log2FC);
geneTableC20.Properties.RowNames = geneCountTable.gene_id;


% DE analysis
diffTableLocalC20 = rnaseqde(geneCountTable,["negative-1","negative-2"],["C20mM-1","C20mM-2"],"VarianceLink","local",IDColumns="gene_id");


lowCountThreshold = 10;
lowCountGenes = all(counts < lowCountThreshold, 2);
% Plot P-Value enrichment
%{
figure
histogram(diffTableLocalC20.PValue,100)
xlabel('P-Value')
ylabel('Frequency')
title('P-value enrichment')

figure
histogram(diffTableLocalC20.PValue(~lowCountGenes),100)
xlabel('P-value')
ylabel('Frequency')
title('P-value enrichment without low count genes')
%}


% Create a table with significant genes
sig = diffTableLocalC20.AdjustedPValue < sig_cutoff;
diffTableLocalSigC20 = diffTableLocalC20(sig,:);
diffTableLocalSigC20 = sortrows(diffTableLocalSigC20,'AdjustedPValue');


% Find all DE genes
deC20 = abs(diffTableLocalSigC20.Log2FoldChange) >= decutoff;
deGenesC20 = sortrows(diffTableLocalSigC20(deC20,:),'Log2FoldChange','descend');
numSigGenesDEC20 = sum(deC20);

% Find up-regulated genes
upC20 = diffTableLocalSigC20.Log2FoldChange >= decutoff;
upGenesC20 = sortrows(diffTableLocalSigC20(upC20,:),'Log2FoldChange','descend');
numSigGenesUpC20 = sum(upC20);


% find down-regulated genes
downC20 = diffTableLocalSigC20.Log2FoldChange <= -decutoff;
downGenesC20 = sortrows(diffTableLocalSigC20(downC20,:),'Log2FoldChange','ascend');
numSigGenesDowC20 = sum(downC20);



%% C40mM Analysis

% Consider the mean
meanUntreated = mean(normCounts(:,1:2), 2);
meanTreated = mean(normCounts(:,5:6), 2);


% Consider the dispersion
dispUntreated = std(normCounts(:,1:2), 0, 2) ./meanUntreated;
dispTreated = std(normCounts(:,5:6), 0, 2) ./meanTreated;


% Plot on loglog scale
%{
figure
loglog(meanTreated, dispTreated, 'r.');
hold on;
loglog(meanUntreated, dispUntreated, 'b.');
xlabel('log2(Mean)');
ylabel('log2(Dispersion)');
legend('Treated','Untreated','Location','southwest')
%}


% Compute the mean and the log2FC
meanBase = (meanTreated + meanUntreated)/2;
foldChange = meanTreated ./ meanUntreated;
log2FC = log2(foldChange);


% Plot mean vs. fold change (MA plot)
%{
mairplot(normCounts(:,3), meanUntreated,'Type','MA','Plotonly',true);
set(get(gca,'Xlabel'),'String','mean of normalized counts')
set(get(gca,'Ylabel'),'String','log2(fold change)')
%}


% Create a table with statistics about each gene
geneTableC40 = table(gene_names, meanBase, meanTreated, meanUntreated, foldChange, log2FC);
geneTableC40.Properties.RowNames = geneCountTable.gene_id;


% DE analysis
diffTableLocalC40 = rnaseqde(geneCountTable,["negative-1","negative-2"],["C40mM-1","C40mM-2"],"VarianceLink","local",IDColumns="gene_id");


lowCountThreshold = 10;
lowCountGenes = all(counts < lowCountThreshold, 2);
% Plot P-Value enrichment
%{
figure
histogram(diffTableLocalC40.PValue,100)
xlabel('P-Value')
ylabel('Frequency')
title('P-value enrichment')

figure
histogram(diffTableLocalC40.PValue(~lowCountGenes),100)
xlabel('P-value')
ylabel('Frequency')
title('P-value enrichment without low count genes')
%}


% Create a table with significant genes
sig = diffTableLocalC40.AdjustedPValue < sig_cutoff;
diffTableLocalSigC40 = diffTableLocalC40(sig,:);
diffTableLocalSigC40 = sortrows(diffTableLocalSigC40,'AdjustedPValue');

% Find all DE genes
deC40 = abs(diffTableLocalSigC40.Log2FoldChange) >= decutoff;
deGenesC40 = sortrows(diffTableLocalSigC40(deC40,:),'Log2FoldChange','descend');
numSigGenesDEC40 = sum(deC40);

% Find up-regulated genes
upC40 = diffTableLocalSigC40.Log2FoldChange >= decutoff;
upGenesC40 = sortrows(diffTableLocalSigC40(upC40,:),'Log2FoldChange','descend');
numSigGenesUpC40 = sum(upC40);

% find down-regulated genes
downC40 = diffTableLocalSigC40.Log2FoldChange <= -decutoff;
downGenesC40 = sortrows(diffTableLocalSigC40(downC40,:),'Log2FoldChange','ascend');
numSigGenesDowC40 = sum(downC40);



%{
figure
scatter(log2(geneTableC40.meanBase), diffTableLocalC40.Log2FoldChange, 3, diffTableLocalC40.AdjustedPValue, 'o')
colormap(flipud(cool(256)))
colorbar;
ylabel('log2(Fold Change)')
xlabel('log2(Mean of normalized counts)')
title('Fold change by FDR')
%}


% Volcano plot for C40 vs negative controls
%{
sig_level = 0.05;
figure
patch([-8 -1 -1 -8],[0 0 -log10(sig_level) -log10(sig_level)], [0.7 0.7 0.7])
hold on
patch([1 8 8 1],[0 0 -log10(sig_level) -log10(sig_level)], [0.7 0.7 0.7])
patch([-1 1 1 -1],[0 0 80 80], [0.7 0.7 0.7])
scatter(diffTableLocalC40.Log2FoldChange, -log10(diffTableLocalC40.AdjustedPValue), 'filled')
plot([-1 -1], [0 80], '--r')
plot([1 1], [0 80], '--r')
plot([-8 8], [-log10(sig_level) -log10(sig_level)], '--r')
xlabel('Log2(Fold Change)')
ylabel('-log10(Padj)')
%}



%% P36SB Analysis

% Consider the mean
meanUntreated = mean(normCounts(:,1:2), 2);
meanTreated = mean(normCounts(:,9:10), 2);

% Consider the dispersion
dispUntreated = std(normCounts(:,1:2), 0, 2) ./meanUntreated;
dispTreated = std(normCounts(:,9:10), 0, 2) ./meanTreated;



% Plot on loglog scale
%{
figure
loglog(meanTreated, dispTreated, 'r.');
hold on;
loglog(meanUntreated, dispUntreated, 'b.');
xlabel('log2(Mean)');
ylabel('log2(Dispersion)');
legend('Treated','Untreated','Location','southwest')
%}


% Compute the mean and the log2FC
meanBase = (meanTreated + meanUntreated)/2;
foldChange = meanTreated ./ meanUntreated;
log2FC = log2(foldChange);


% Plot mean vs. fold change (MA plot)
%{
mairplot(normCounts(:,3), meanUntreated,'Type','MA','Plotonly',true);
set(get(gca,'Xlabel'),'String','mean of normalized counts')
set(get(gca,'Ylabel'),'String','log2(fold change)')
%}


% Create a table with statistics about each gene
geneTable36SB = table(gene_names, meanBase, meanTreated, meanUntreated, foldChange, log2FC);
geneTable36SB.Properties.RowNames = geneCountTable.gene_id;


% DE analysis
diffTableLocal36SB = rnaseqde(geneCountTable,["negative-1","negative-2"],["P36SB-1","P36SB-2"],"VarianceLink","local",IDColumns="gene_id");


lowCountThreshold = 10;
lowCountGenes = all(counts < lowCountThreshold, 2);
% Plot P-Value enrichment
%{
figure
histogram(diffTableLocalC120.PValue,100)
xlabel('P-Value')
ylabel('Frequency')
title('P-value enrichment')

figure
histogram(diffTableLocalC120.PValue(~lowCountGenes),100)
xlabel('P-value')
ylabel('Frequency')
title('P-value enrichment without low count genes')
%}


% Create a table with significant genes
sig = diffTableLocal36SB.AdjustedPValue < sig_cutoff;
diffTableLocalSig36SB = diffTableLocal36SB(sig,:);
diffTableLocalSig36SB = sortrows(diffTableLocalSig36SB,'AdjustedPValue');

% Find all DE genes
de36SB = abs(diffTableLocalSig36SB.Log2FoldChange) >= decutoff;
deGenes36SB = sortrows(diffTableLocalSig36SB(de36SB,:),'Log2FoldChange','descend');
numSigGenesDE36SB = sum(de36SB);

% Find up-regulated genes
up36SB = diffTableLocalSig36SB.Log2FoldChange >= decutoff ;
upGenes36SB = sortrows(diffTableLocalSig36SB(up36SB,:),'Log2FoldChange','descend');
numSigGenesUp36SB = sum(up36SB);


% find down-regulated genes
down36SB = diffTableLocalSig36SB.Log2FoldChange <= -decutoff;
downGenes36SB = sortrows(diffTableLocalSig36SB(down36SB,:),'Log2FoldChange','ascend');
numSigGenesDow36SB = sum(down36SB);



%% 0.09 Analysis

% Consider the mean
meanUntreated = mean(normCounts(:,1:2), 2);
meanTreated = mean(normCounts(:,11:12), 2);

% Consider the dispersion
dispUntreated = std(normCounts(:,1:2), 0, 2) ./meanUntreated;
dispTreated = std(normCounts(:,11:12), 0, 2) ./meanTreated;


% Plot on loglog scale
%{
figure
loglog(meanTreated, dispTreated, 'r.');
hold on;
loglog(meanUntreated, dispUntreated, 'b.');
xlabel('log2(Mean)');
ylabel('log2(Dispersion)');
legend('Treated','Untreated','Location','southwest')
%}


% Compute the mean and the log2FC
meanBase = (meanTreated + meanUntreated)/2;
foldChange = meanTreated ./ meanUntreated;
log2FC = log2(foldChange);


% Plot mean vs. fold change (MA plot)
%{
mairplot(normCounts(:,3), meanUntreated,'Type','MA','Plotonly',true);
set(get(gca,'Xlabel'),'String','mean of normalized counts')
set(get(gca,'Ylabel'),'String','log2(fold change)')
%}


% Create a table with statistics about each gene
geneTable09 = table(gene_names, meanBase, meanTreated, meanUntreated, foldChange, log2FC);
geneTable09.Properties.RowNames = geneCountTable.gene_id;


% DE analysis
diffTableLocal09 = rnaseqde(geneCountTable,["negative-1","negative-2"],["P9-1","P9-2"],"VarianceLink","local",IDColumns="gene_id");

lowCountThreshold = 10;
lowCountGenes = all(counts < lowCountThreshold, 2);



% Create a table with significant genes
sig = diffTableLocal09.AdjustedPValue < sig_cutoff;
diffTableLocalSig09 = diffTableLocal09(sig,:);
diffTableLocalSig09 = sortrows(diffTableLocalSig09,'AdjustedPValue');

% Find All 09 DE genes
de09 = abs(diffTableLocalSig09.Log2FoldChange) >= decutoff;
DEGenes09 = sortrows(diffTableLocalSig09(de09,:),'Log2FoldChange','descend');
numSigGenesDE09 = sum(de09);

% Find up-regulated genes
up09 = diffTableLocalSig09.Log2FoldChange >= decutoff;
upGenes09 = sortrows(diffTableLocalSig09(up09,:),'Log2FoldChange','descend');
numSigGenesUp09 = sum(up09);


% find down-regulated genes
down09 = diffTableLocalSig09.Log2FoldChange <= -decutoff;
downGenes09 = sortrows(diffTableLocalSig09(down09,:),'Log2FoldChange','ascend');
numSigGenesDow09 = sum(down09);




%% 0.18 Analysis

% Consider the mean
meanUntreated = mean(normCounts(:,1:2), 2);
meanTreated = mean(normCounts(:,14:15), 2);

% Consider the dispersion
dispUntreated = std(normCounts(:,1:2), 0, 2) ./meanUntreated;
dispTreated = std(normCounts(:,14:15), 0, 2) ./meanTreated;


% Plot on loglog scale
%{
figure
loglog(meanTreated, dispTreated, 'r.');
hold on;
loglog(meanUntreated, dispUntreated, 'b.');
xlabel('log2(Mean)');
ylabel('log2(Dispersion)');
legend('Treated','Untreated','Location','southwest')
%}


% Compute the mean and the log2FC
meanBase = (meanTreated + meanUntreated)/2;
foldChange = meanTreated ./ meanUntreated;
log2FC = log2(foldChange);


% Plot mean vs. fold change (MA plot)
%{
mairplot(normCounts(:,3), meanUntreated,'Type','MA','Plotonly',true);
set(get(gca,'Xlabel'),'String','mean of normalized counts')
set(get(gca,'Ylabel'),'String','log2(fold change)')
%}


% Create a table with statistics about each gene
geneTable18 = table(gene_names, meanBase, meanTreated, meanUntreated, foldChange, log2FC);
geneTable18.Properties.RowNames = geneCountTable.gene_id;


% DE analysis
diffTableLocal18 = rnaseqde(geneCountTable,["negative-1","negative-2"],["P18-1","P18-2"],"VarianceLink","local",IDColumns="gene_id");

lowCountThreshold = 10;
lowCountGenes = all(counts < lowCountThreshold, 2);



% Create a table with significant genes
sig = diffTableLocal18.AdjustedPValue < sig_cutoff;
diffTableLocalSig18 = diffTableLocal18(sig,:);
diffTableLocalSig18 = sortrows(diffTableLocalSig18,'AdjustedPValue');

% Find all DE genes
de18 = abs(diffTableLocalSig18.Log2FoldChange) >= decutoff;
deGenes18 = sortrows(diffTableLocalSig18(de18,:),'Log2FoldChange','descend');
numSigGenesDE18 = sum(de18);

% Find up-regulated genes
up18 = diffTableLocalSig18.Log2FoldChange >= decutoff;
upGenes18 = sortrows(diffTableLocalSig18(up18,:),'Log2FoldChange','descend');
numSigGenesUp18 = sum(up18);


% find down-regulated genes
down18 = diffTableLocalSig18.Log2FoldChange <= -decutoff;
downGenes18 = sortrows(diffTableLocalSig18(down18,:),'Log2FoldChange','ascend');
numSigGenesDow18 = sum(down18);




%% 0.27 Analysis

% Consider the mean
meanUntreated = mean(normCounts(:,1:2), 2);
meanTreated = mean(normCounts(:,16:17), 2);

% Consider the dispersion
dispUntreated = std(normCounts(:,1:2), 0, 2) ./meanUntreated;
dispTreated = std(normCounts(:,16:17), 0, 2) ./meanTreated;


% Plot on loglog scale
%{
figure
loglog(meanTreated, dispTreated, 'r.');
hold on;
loglog(meanUntreated, dispUntreated, 'b.');
xlabel('log2(Mean)');
ylabel('log2(Dispersion)');
legend('Treated','Untreated','Location','southwest')
%}


% Compute the mean and the log2FC
meanBase = (meanTreated + meanUntreated)/2;
foldChange = meanTreated ./ meanUntreated;
log2FC = log2(foldChange);


% Plot mean vs. fold change (MA plot)
%{
mairplot(normCounts(:,3), meanUntreated,'Type','MA','Plotonly',true);
set(get(gca,'Xlabel'),'String','mean of normalized counts')
set(get(gca,'Ylabel'),'String','log2(fold change)')
%}


% Create a table with statistics about each gene
geneTable27 = table(gene_names, meanBase, meanTreated, meanUntreated, foldChange, log2FC);
geneTable27.Properties.RowNames = geneCountTable.gene_id;


% DE analysis
diffTableLocal27 = rnaseqde(geneCountTable,["negative-1","negative-2"],["P27-1","P27-2"],"VarianceLink","local",IDColumns="gene_id");

lowCountThreshold = 10;
lowCountGenes = all(counts < lowCountThreshold, 2);



% Create a table with significant genes
sig = diffTableLocal27.AdjustedPValue < sig_cutoff;
diffTableLocalSig27 = diffTableLocal27(sig,:);
diffTableLocalSig27 = sortrows(diffTableLocalSig27,'AdjustedPValue');

% Find all DE genes
de27 = abs(diffTableLocalSig27.Log2FoldChange) >= decutoff;
deGenes27 = sortrows(diffTableLocalSig27(de27,:),'Log2FoldChange','descend');
numSigGenesDE27 = sum(de27);

% Find up-regulated genes
up27 = diffTableLocalSig27.Log2FoldChange >= decutoff;
upGenes27 = sortrows(diffTableLocalSig27(up27,:),'Log2FoldChange','descend');
numSigGenesUp27 = sum(up27);


% find down-regulated genes
down27 = diffTableLocalSig27.Log2FoldChange <= -decutoff;
downGenes27 = sortrows(diffTableLocalSig27(down27,:),'Log2FoldChange','ascend');
numSigGenesDow27 = sum(down27);




%% 0.36 Analysis

% Consider the mean
meanUntreated = mean(normCounts(:,1:2),2);
meanTreated = mean(normCounts(:,18:19), 2);

% Consider the dispersion
dispUntreated = std(normCounts(:,1:2), 0, 2) ./meanUntreated;
dispTreated = std(normCounts(:,18:19), 0, 2) ./meanTreated;


% Plot on loglog scale
%{
figure
loglog(meanTreated, dispTreated, 'r.');
hold on;
loglog(meanUntreated, dispUntreated, 'b.');
xlabel('log2(Mean)');
ylabel('log2(Dispersion)');
legend('Treated','Untreated','Location','southwest')
%}


% Compute the mean and the log2FC
meanBase = (meanTreated + meanUntreated)/2;
foldChange = meanTreated ./ meanUntreated;
log2FC = log2(foldChange);


% Plot mean vs. fold change (MA plot)
%{
mairplot(normCounts(:,3), meanUntreated,'Type','MA','Plotonly',true);
set(get(gca,'Xlabel'),'String','mean of normalized counts')
set(get(gca,'Ylabel'),'String','log2(fold change)')
%}


% Create a table with statistics about each gene
geneTable36 = table(gene_names, meanBase, meanTreated, meanUntreated, foldChange, log2FC);
geneTable36.Properties.RowNames = geneCountTable.gene_id;


% DE analysis
diffTableLocal36 = rnaseqde(geneCountTable,["negative-1","negative-2"],["P36-1","P36-2"],"VarianceLink","local",IDColumns="gene_id");


lowCountThreshold = 10;
lowCountGenes = all(counts < lowCountThreshold, 2);
% Plot P-Value enrichment
%{
figure
histogram(diffTableLocal36.PValue,100)
xlabel('P-Value')
ylabel('Frequency')
title('P-value enrichment')

figure
histogram(diffTableLocal36.PValue(~lowCountGenes),100)
xlabel('P-value')
ylabel('Frequency')
title('P-value enrichment without low count genes')
%}


% Create a table with significant genes
sig = diffTableLocal36.AdjustedPValue < sig_cutoff;
diffTableLocalSig36 = diffTableLocal36(sig,:);
diffTableLocalSig36 = sortrows(diffTableLocalSig36,'AdjustedPValue');

% Find all DE genes
de36 = abs(diffTableLocalSig36.Log2FoldChange) >= decutoff;
deGenes36 = sortrows(diffTableLocalSig36(de36,:),'Log2FoldChange','descend');
numSigGenesDE36 = sum(de36);

% Find up-regulated genes
up36 = diffTableLocalSig36.Log2FoldChange >= decutoff;
upGenes36 = sortrows(diffTableLocalSig36(up36,:),'Log2FoldChange','descend');
numSigGenesUp36 = sum(up36);

% find down-regulated genes
down36 = diffTableLocalSig36.Log2FoldChange <= -decutoff;
downGenes36 = sortrows(diffTableLocalSig36(down36,:),'Log2FoldChange','ascend');
numSigGenesDow36 = sum(down36);


%{
figure
scatter(log2(geneTable36.meanBase), diffTableLocal36.Log2FoldChange, 3, diffTableLocal36.AdjustedPValue, 'o')
colormap(flipud(cool(256)))
colorbar;
ylabel('log2(Fold Change)')
xlabel('log2(Mean of normalized counts)')
title('Fold change by FDR')
%}


% Volcano plot for 0.36 vs negative controls
%{
sig_level = 0.05;
figure
patch([-8 -1 -1 -8],[0 0 -log10(sig_level) -log10(sig_level)], [0.7 0.7 0.7])
hold on
patch([1 8 8 1],[0 0 -log10(sig_level) -log10(sig_level)], [0.7 0.7 0.7])
patch([-1 1 1 -1],[0 0 80 80], [0.7 0.7 0.7])
scatter(diffTableLocal36.Log2FoldChange, -log10(diffTableLocal36.AdjustedPValue), 'filled')
plot([-1 -1], [0 80], '--r')
plot([1 1], [0 80], '--r')
plot([-8 8], [-log10(sig_level) -log10(sig_level)], '--r')
xlabel('Log2(Fold Change)')
ylabel('-log10(Padj)')
%}



%% 0.45 Analysis

% Consider the mean
meanUntreated = mean(normCounts(:,1:2), 2);
meanTreated = mean(normCounts(:,20:21), 2);

% Consider the dispersion
dispUntreated = std(normCounts(:,1:2), 0, 2) ./meanUntreated;
dispTreated = std(normCounts(:,20:21), 0, 2) ./meanTreated;


% Plot on loglog scale
%{
figure
loglog(meanTreated, dispTreated, 'r.');
hold on;
loglog(meanUntreated, dispUntreated, 'b.');
xlabel('log2(Mean)');
ylabel('log2(Dispersion)');
legend('Treated','Untreated','Location','southwest')
%}


% Compute the mean and the log2FC
meanBase = (meanTreated + meanUntreated)/2;
foldChange = meanTreated ./ meanUntreated;
log2FC = log2(foldChange);


% Plot mean vs. fold change (MA plot)
%{
mairplot(normCounts(:,3), meanUntreated,'Type','MA','Plotonly',true);
set(get(gca,'Xlabel'),'String','mean of normalized counts')
set(get(gca,'Ylabel'),'String','log2(fold change)')
%}


% Create a table with statistics about each gene
geneTable45 = table(gene_names, meanBase, meanTreated, meanUntreated, foldChange, log2FC);
geneTable45.Properties.RowNames = geneCountTable.gene_id;


% DE analysis
diffTableLocal45 = rnaseqde(geneCountTable,["negative-1","negative-2"],["P45-1","P45-2"],"VarianceLink","local",IDColumns="gene_id");

lowCountThreshold = 10;
lowCountGenes = all(counts < lowCountThreshold, 2);



% Create a table with significant genes
sig = diffTableLocal45.AdjustedPValue < sig_cutoff;
diffTableLocalSig45 = diffTableLocal45(sig,:);
diffTableLocalSig45 = sortrows(diffTableLocalSig45,'AdjustedPValue');

% Find all DE genes
de45 = abs(diffTableLocalSig45.Log2FoldChange) >= decutoff;
deGenes45 = sortrows(diffTableLocalSig45(de45,:),'Log2FoldChange','descend');
numSigGenesDE45 = sum(de45);

% Find up-regulated genes
up45 = diffTableLocalSig45.Log2FoldChange >= decutoff;
upGenes45 = sortrows(diffTableLocalSig45(up45,:),'Log2FoldChange','descend');
numSigGenesUp45 = sum(up45);


% find down-regulated genes
down45 = diffTableLocalSig45.Log2FoldChange <= -decutoff;
downGenes45 = sortrows(diffTableLocalSig45(down45,:),'Log2FoldChange','ascend');
numSigGenesDow45 = sum(down45);



%% 0.54 Analysis

% Consider the mean
meanUntreated = mean(normCounts(:,1:2),2);
meanTreated = mean(normCounts(:,22:23),2);

% Consider the dispersion
dispUntreated = std(normCounts(:,1:2),0,2) ./meanUntreated;
dispTreated = std(normCounts(:,22:23),0,2) ./meanTreated;


% Plot on loglog scale
%{
figure
loglog(meanTreated, dispTreated, 'r.');
hold on;
loglog(meanUntreated, dispUntreated, 'b.');
xlabel('log2(Mean)');
ylabel('log2(Dispersion)');
legend('Treated','Untreated','Location','southwest')
%}


% Compute the mean and the log2FC
meanBase = (meanTreated + meanUntreated)/2;
foldChange = meanTreated ./ meanUntreated;
log2FC = log2(foldChange);

% Plot mean vs. fold change (MA plot)
%{
mairplot(normCounts(:,3), meanUntreated,'Type','MA','Plotonly',true);
set(get(gca,'Xlabel'),'String','mean of normalized counts')
set(get(gca,'Ylabel'),'String','log2(fold change)')
%}


% Create a table with statistics about each gene
geneTable54 = table(gene_names, meanBase, meanTreated, meanUntreated, foldChange, log2FC);
geneTable54.Properties.RowNames = geneCountTable.gene_id;


diffTableLocal54 = rnaseqde(geneCountTable,["negative-1","negative-2"],["P54-1","P54-2"],"VarianceLink","local",IDColumns="gene_id");


sig = diffTableLocal54.AdjustedPValue < sig_cutoff;
diffTableLocalSig54 = diffTableLocal54(sig,:);
diffTableLocalSig54 = sortrows(diffTableLocalSig54,'AdjustedPValue');

% Find all DE genes
de54 = abs(diffTableLocalSig54.Log2FoldChange) >= decutoff;
deGenes54 = sortrows(diffTableLocalSig54(de54,:),'Log2FoldChange','descend');
numSigGenesDE54 = sum(de54);

% Find up-regulated genes
up54 = diffTableLocalSig54.Log2FoldChange >= decutoff;
upGenes54 = sortrows(diffTableLocalSig54(up54,:),'Log2FoldChange','descend');
numSigGenesUp54 = sum(up54);

% Display the top 10 up-regulated genes
top10GenesUp54 = upGenes54(14:33,:);


% find down-regulated genes
down54 = diffTableLocalSig54.Log2FoldChange <= -decutoff;
downGenes54 = sortrows(diffTableLocalSig54(down54,:),'Log2FoldChange','ascend');
numSigGenesDow54 = sum(down54);

% find top 10 down-regulated genes
top10GenesDown54 = downGenes54(66:85,:);


% Volcano plot for 0.54 vs negative controls
%{
sig_level = 0.05;
figure
patch([-8 -1 -1 -8],[0 0 -log10(sig_level) -log10(sig_level)], [0.7 0.7 0.7])
hold on
patch([1 8 8 1],[0 0 -log10(sig_level) -log10(sig_level)], [0.7 0.7 0.7])
patch([-1 1 1 -1],[0 0 80 80], [0.7 0.7 0.7])
scatter(diffTableLocal54.Log2FoldChange, -log10(diffTableLocal54.AdjustedPValue), 'filled')
plot([-1 -1], [0 80], '--r')
plot([1 1], [0 80], '--r')
plot([-8 8], [-log10(sig_level) -log10(sig_level)], '--r')
xlabel('Log2(Fold Change)')
ylabel('-log10(Padj)')
%}




%% 0.63 Analysis

% Consider the mean
meanUntreated = mean(normCounts(:,1:2),2);
meanTreated = mean(normCounts(:,24:25),2);

% Consider the dispersion
dispUntreated = std(normCounts(:,1:2),0,2) ./meanUntreated;
dispTreated = std(normCounts(:,24:25),0,2) ./meanTreated;


% Plot on loglog scale
%{
figure
loglog(meanTreated, dispTreated, 'r.');
hold on;
loglog(meanUntreated, dispUntreated, 'b.');
xlabel('log2(Mean)');
ylabel('log2(Dispersion)');
legend('Treated','Untreated','Location','southwest')
%}


% Compute the mean and the log2FC
meanBase = (meanTreated + meanUntreated)/2;
foldChange = meanTreated ./ meanUntreated;
log2FC = log2(foldChange);

% Plot mean vs. fold change (MA plot)
%{
mairplot(normCounts(:,3), meanUntreated,'Type','MA','Plotonly',true);
set(get(gca,'Xlabel'),'String','mean of normalized counts')
set(get(gca,'Ylabel'),'String','log2(fold change)')
%}


% Create a table with statistics about each gene
geneTable63 = table(gene_names, meanBase, meanTreated, meanUntreated, foldChange, log2FC);
geneTable63.Properties.RowNames = geneCountTable.gene_id;


diffTableLocal63 = rnaseqde(geneCountTable,["negative-1","negative-2"],["P63-1","P63-2"],"VarianceLink","local",IDColumns="gene_id");


sig = diffTableLocal63.AdjustedPValue < sig_cutoff;
diffTableLocalSig63 = diffTableLocal63(sig,:);
diffTableLocalSig63 = sortrows(diffTableLocalSig63,'AdjustedPValue');

% Find all DE genes
de63 = abs(diffTableLocalSig63.Log2FoldChange) >= decutoff;
deGenes63 = sortrows(diffTableLocalSig63(de63,:),'Log2FoldChange','descend');
numSigGenesDE63 = sum(de63);

% Find up-regulated genes
up63 = diffTableLocalSig63.Log2FoldChange >= decutoff;
upGenes63 = sortrows(diffTableLocalSig63(up63,:),'Log2FoldChange','descend');
numSigGenesUp63 = sum(up63);


% find down-regulated genes
down63 = diffTableLocalSig63.Log2FoldChange <= -decutoff;
downGenes63 = sortrows(diffTableLocalSig63(down63,:),'Log2FoldChange','ascend');
numSigGenesDow63 = sum(down63);



%% 0.72 Analysis

% Consider the mean
meanUntreated = mean(normCounts(:,1:2),2);
meanTreated = mean(normCounts(:,26:27),2);

% Consider the dispersion
dispUntreated = std(normCounts(:,1:2),0,2) ./meanUntreated;
dispTreated = std(normCounts(:,26:27),0,2) ./meanTreated;


% Plot on loglog scale
%{
figure
loglog(meanTreated, dispTreated, 'r.');
hold on;
loglog(meanUntreated, dispUntreated, 'b.');
xlabel('log2(Mean)');
ylabel('log2(Dispersion)');
legend('Treated','Untreated','Location','southwest')
%}


% Compute the mean and the log2FC
meanBase = (meanTreated + meanUntreated)/2;
foldChange = meanTreated ./ meanUntreated;
log2FC = log2(foldChange);

% Plot mean vs. fold change (MA plot)
%{
mairplot(normCounts(:,3), meanUntreated,'Type','MA','Plotonly',true);
set(get(gca,'Xlabel'),'String','mean of normalized counts')
set(get(gca,'Ylabel'),'String','log2(fold change)')
%}


% Create a table with statistics about each gene
geneTable72 = table(gene_names, meanBase, meanTreated, meanUntreated, foldChange, log2FC);
geneTable72.Properties.RowNames = geneCountTable.gene_id;


diffTableLocal72 = rnaseqde(geneCountTable,["negative-1","negative-2"],["P72-1","P72-2"],"VarianceLink","local",IDColumns="gene_id");


sig = diffTableLocal72.AdjustedPValue < sig_cutoff;
diffTableLocalSig72 = diffTableLocal72(sig,:);
diffTableLocalSig72 = sortrows(diffTableLocalSig72,'AdjustedPValue');

% Find all DE genes
de72 = abs(diffTableLocalSig72.Log2FoldChange) >= decutoff;
deGenes72 = sortrows(diffTableLocalSig72(de72,:),'Log2FoldChange','descend');
numSigGenesDE72 = sum(de72);

% Find up-regulated genes
up72 = diffTableLocalSig72.Log2FoldChange >= decutoff;
upGenes72 = sortrows(diffTableLocalSig72(up72,:),'Log2FoldChange','descend');
numSigGenesUp72 = sum(up72);


% find down-regulated genes
down72 = diffTableLocalSig72.Log2FoldChange <= -decutoff;
downGenes72 = sortrows(diffTableLocalSig72(down72,:),'Log2FoldChange','ascend');
numSigGenesDow72 = sum(down72);




%% 0.81 Analysis

% Consider the mean
meanUntreated = mean(normCounts(:,1:2),2);
meanTreated = mean(normCounts(:,28:29),2);

% Consider the dispersion
dispUntreated = std(normCounts(:,1:2),0,2) ./meanUntreated;
dispTreated = std(normCounts(:,28:29),0,2) ./meanTreated;


% Plot on loglog scale
%{
figure
loglog(meanTreated, dispTreated, 'r.');
hold on;
loglog(meanUntreated, dispUntreated, 'b.');
xlabel('log2(Mean)');
ylabel('log2(Dispersion)');
legend('Treated','Untreated','Location','southwest')
%}


% Compute the mean and the log2FC
meanBase = (meanTreated + meanUntreated)/2;
foldChange = meanTreated ./ meanUntreated;
log2FC = log2(foldChange);

% Plot mean vs. fold change (MA plot)
%{
mairplot(normCounts(:,3), meanUntreated,'Type','MA','Plotonly',true);
set(get(gca,'Xlabel'),'String','mean of normalized counts')
set(get(gca,'Ylabel'),'String','log2(fold change)')
%}


% Create a table with statistics about each gene
geneTable81 = table(gene_names, meanBase, meanTreated, meanUntreated, foldChange, log2FC);
geneTable81.Properties.RowNames = geneCountTable.gene_id;


diffTableLocal81 = rnaseqde(geneCountTable,["negative-1","negative-2"],["P81-1","P81-2"],"VarianceLink","local",IDColumns="gene_id");


sig = diffTableLocal81.AdjustedPValue < sig_cutoff;
diffTableLocalSig81 = diffTableLocal81(sig,:);
diffTableLocalSig81 = sortrows(diffTableLocalSig81,'AdjustedPValue');

% Find all DE genes
de81 = abs(diffTableLocalSig81.Log2FoldChange) >= decutoff;
deGenes81 = sortrows(diffTableLocalSig81(de81,:),'Log2FoldChange','descend');
numSigGenesDE81 = sum(de81);

% Find up-regulated genes
up81 = diffTableLocalSig81.Log2FoldChange >= decutoff;
upGenes81 = sortrows(diffTableLocalSig81(up81,:),'Log2FoldChange','descend');
numSigGenesUp81 = sum(up81);


% find down-regulated genes
down81 = diffTableLocalSig81.Log2FoldChange <= -decutoff;
downGenes81 = sortrows(diffTableLocalSig81(down81,:),'Log2FoldChange','ascend');
numSigGenesDow81 = sum(down81);




%% 0.90 Analysis

% Consider the mean
meanUntreated = mean(normCounts(:,1:2),2);
meanTreated = mean(normCounts(:,30:31),2);

% Consider the dispersion
dispUntreated = std(normCounts(:,1:2),0,2) ./meanUntreated;
dispTreated = std(normCounts(:,30:31),0,2) ./meanTreated;


% Plot on loglog scale
%{
figure
loglog(meanTreated, dispTreated, 'r.');
hold on;
loglog(meanUntreated, dispUntreated, 'b.');
xlabel('log2(Mean)');
ylabel('log2(Dispersion)');
legend('Treated','Untreated','Location','southwest')
%}


% Compute the mean and the log2FC
meanBase = (meanTreated + meanUntreated)/2;
foldChange = meanTreated ./ meanUntreated;
log2FC = log2(foldChange);

% Plot mean vs. fold change (MA plot)
%{
mairplot(normCounts(:,3), meanUntreated,'Type','MA','Plotonly',true);
set(get(gca,'Xlabel'),'String','mean of normalized counts')
set(get(gca,'Ylabel'),'String','log2(fold change)')
%}


% Create a table with statistics about each gene
geneTable90 = table(gene_names, meanBase, meanTreated, meanUntreated, foldChange, log2FC);
geneTable90.Properties.RowNames = geneCountTable.gene_id;


diffTableLocal90 = rnaseqde(geneCountTable,["negative-1","negative-2"],["P90-1","P90-2"],"VarianceLink","local",IDColumns="gene_id");


sig = diffTableLocal90.AdjustedPValue < sig_cutoff;
diffTableLocalSig90 = diffTableLocal90(sig,:);
diffTableLocalSig90 = sortrows(diffTableLocalSig90,'AdjustedPValue');

% Find all DE genes
de90 = abs(diffTableLocalSig90.Log2FoldChange) >= decutoff;
deGenes90 = sortrows(diffTableLocalSig90(de90,:),'Log2FoldChange','descend');
numSigGenesDE90 = sum(de90);

% Find up-regulated genes
up90 = diffTableLocalSig90.Log2FoldChange >= decutoff;
upGenes90 = sortrows(diffTableLocalSig90(up90,:),'Log2FoldChange','descend');
numSigGenesUp90 = sum(up90);


% find down-regulated genes
down90 = diffTableLocalSig90.Log2FoldChange <= -decutoff;
downGenes90 = sortrows(diffTableLocalSig90(down90,:),'Log2FoldChange','ascend');
numSigGenesDow90 = sum(down90);



%% Find unique genes
DEC20Table2 = diffTableLocalC20;
DEC40Table2 = diffTableLocalC40;
DE36SBTable2 = diffTableLocal36SB;
DE09Table2 = diffTableLocal09;
DE18Table2 = diffTableLocal18;
DE27Table2 = diffTableLocal27;
DE36Table2 = diffTableLocal36;
DE45Table2 = diffTableLocal45;
DE54Table2 = diffTableLocal54;
DE63Table2 = diffTableLocal63;
DE72Table2 = diffTableLocal72;
DE81Table2 = diffTableLocal81;
DE90Table2 = diffTableLocal90;

DEC20Table2.Properties.VariableNames = ["gene_id", "negMeanC20", "expMeanC20", "Log2FCC20", "PValueC20", "PadjC20"];
DEC40Table2.Properties.VariableNames = ["gene_id", "negMeanC40", "expMeanC40", "Log2FCC40", "PValueC40", "PadjC40"];
DE36SBTable2.Properties.VariableNames = ["gene_id", "negMean36SB", "expMean36SB", "Log2FC36SB", "PValue36SB", "Padj36SB"];
DE09Table2.Properties.VariableNames = ["gene_id", "negMean09", "expMean09", "Log2FC09", "PValue09", "Padj09"];
DE18Table2.Properties.VariableNames = ["gene_id", "negMean18", "expMean18", "Log2FC18", "PValue18", "Padj18"];
DE27Table2.Properties.VariableNames = ["gene_id", "negMean27", "expMean27", "Log2FC27", "PValue27", "Padj27"];
DE36Table2.Properties.VariableNames = ["gene_id", "negMean36", "expMean36", "Log2FC36", "PValue36", "Padj36"];
DE45Table2.Properties.VariableNames = ["gene_id", "negMean45", "expMean45", "Log2FC45", "PValue45", "Padj45"];
DE54Table2.Properties.VariableNames = ["gene_id", "negMean54", "expMean54", "Log2FC54", "PValue54", "Padj54"];
DE63Table2.Properties.VariableNames = ["gene_id", "negMean63", "expMean63", "Log2FC63", "PValue63", "Padj63"];
DE72Table2.Properties.VariableNames = ["gene_id", "negMean72", "expMean72", "Log2FC72", "PValue72", "Padj72"];
DE81Table2.Properties.VariableNames = ["gene_id", "negMean81", "expMean81", "Log2FC81", "PValue81", "Padj81"];
DE90Table2.Properties.VariableNames = ["gene_id", "negMean90", "expMean90", "Log2FC90", "PValue90", "Padj90"];


DEGenes3 = innerjoin(DEC20Table2, DEC40Table2);
DEGenes3 = innerjoin(DEGenes3, DE36SBTable2);
DEGenes3 = innerjoin(DEGenes3, DE09Table2);
DEGenes3 = innerjoin(DEGenes3, DE18Table2);
DEGenes3 = innerjoin(DEGenes3, DE27Table2);
DEGenes3 = innerjoin(DEGenes3, DE36Table2);
DEGenes3 = innerjoin(DEGenes3, DE45Table2);
DEGenes3 = innerjoin(DEGenes3, DE54Table2);
DEGenes3 = innerjoin(DEGenes3, DE63Table2);
DEGenes3 = innerjoin(DEGenes3, DE72Table2);
DEGenes3 = innerjoin(DEGenes3, DE81Table2);
DEGenes3 = innerjoin(DEGenes3, DE90Table2);

DEGenes3 = addvars(DEGenes3, gene_names, 'After', 1);
DEGenesClean3 = rmmissing(DEGenes3);

sig_cutoff = 0.05;
FC_cutoff = log2(2);

sig09 = zeros(size(DEGenesClean3,1),1);
sig18 = zeros(size(DEGenesClean3,1),1);
sig27 = zeros(size(DEGenesClean3,1),1);
sig36 = zeros(size(DEGenesClean3,1),1);
sig45 = zeros(size(DEGenesClean3,1),1);
sig54 = zeros(size(DEGenesClean3,1),1);
sig63 = zeros(size(DEGenesClean3,1),1);
sig72 = zeros(size(DEGenesClean3,1),1);
sig81 = zeros(size(DEGenesClean3,1),1);
sig90 = zeros(size(DEGenesClean3,1),1);


for i=1:size(DEGenesClean3,1)

    if ((DEGenesClean3.Padj09(i) < sig_cutoff) & (abs(DEGenesClean3.Log2FC09(i)) > FC_cutoff))
        sig09(i) = 1;
    end
    if ((DEGenesClean3.Padj18(i) < sig_cutoff) & (abs(DEGenesClean3.Log2FC18(i)) > FC_cutoff))
        sig18(i) = 1;
    end
    if ((DEGenesClean3.Padj27(i) < sig_cutoff) & (abs(DEGenesClean3.Log2FC27(i)) > FC_cutoff))
        sig27(i) = 1;
    end
    if ((DEGenesClean3.Padj36(i) < sig_cutoff) & (abs(DEGenesClean3.Log2FC36(i)) > FC_cutoff))
        sig36(i) = 1;
    end
    if ((DEGenesClean3.Padj45(i) < sig_cutoff) & (abs(DEGenesClean3.Log2FC45(i)) > FC_cutoff))
        sig45(i) = 1;
    end
    if ((DEGenesClean3.Padj54(i) < sig_cutoff) & (abs(DEGenesClean3.Log2FC54(i)) > FC_cutoff))
        sig54(i) = 1;
    end
    if ((DEGenesClean3.Padj63(i) < sig_cutoff) & (abs(DEGenesClean3.Log2FC63(i)) > FC_cutoff))
        sig63(i) = 1;
    end
    if ((DEGenesClean3.Padj72(i) < sig_cutoff) & (abs(DEGenesClean3.Log2FC72(i)) > FC_cutoff))
        sig72(i) = 1;
    end
    if ((DEGenesClean3.Padj81(i) < sig_cutoff) & (abs(DEGenesClean3.Log2FC81(i)) > FC_cutoff))
        sig81(i) = 1;
    end
    if ((DEGenesClean3.Padj90(i) < sig_cutoff) & (abs(DEGenesClean3.Log2FC90(i)) > FC_cutoff))
        sig90(i) = 1;
    end

    bin_code{i} = sprintf('%d', [sig09(i), sig18(i), sig27(i), sig36(i), sig45(i), sig54(i),...
        sig63(i), sig72(i), sig81(i), sig90(i)]);

end

bin_code = string(bin_code);
DEGenesClean3 = addvars(DEGenesClean3, bin_code', 'After', 2);

DEGenesClean3 = renamevars(DEGenesClean3, 'Var3', 'bin_code');

Only09 = DEGenesClean3.bin_code == '1000000000';
Only18 = DEGenesClean3.bin_code == '0100000000';
Only27 = DEGenesClean3.bin_code == '0010000000';
Only36 = DEGenesClean3.bin_code == '0001000000';
Only45 = DEGenesClean3.bin_code == '0000100000';
Only54 = DEGenesClean3.bin_code == '0000010000';
Only63 = DEGenesClean3.bin_code == '0000001000';
Only72 = DEGenesClean3.bin_code == '0000000100';
Only81 = DEGenesClean3.bin_code == '0000000010';
Only90 = DEGenesClean3.bin_code == '0000000001';


Only09Table = DEGenesClean3(Only09,:);
Only09Table = sortrows(Only09Table,'Log2FC09','descend');
Only18Table = DEGenesClean3(Only18,:);
Only18Table = sortrows(Only18Table,'Log2FC18','descend');
Only27Table = DEGenesClean3(Only27,:);
Only27Table = sortrows(Only27Table,'Log2FC27','descend');
Only36Table = DEGenesClean3(Only36,:);
Only36Table = sortrows(Only36Table,'Log2FC36','descend');
Only45Table = DEGenesClean3(Only45,:);
Only45Table = sortrows(Only45Table,'Log2FC45','descend');
Only54Table = DEGenesClean3(Only54,:);
Only54Table = sortrows(Only54Table,'Log2FC54','descend');
Only63Table = DEGenesClean3(Only63,:);
Only63Table = sortrows(Only63Table,'Log2FC63','descend');
Only72Table = DEGenesClean3(Only72,:);
Only72Table = sortrows(Only72Table,'Log2FC72','descend');
Only81Table = DEGenesClean3(Only81,:);
Only81Table = sortrows(Only81Table,'Log2FC81','descend');
Only90Table = DEGenesClean3(Only90,:);
Only90Table = sortrows(Only90Table,'Log2FC90','descend');


char_bin_code = char(DEGenesClean3.bin_code);
Any09 = char_bin_code(:,1) == '1';
Any09Table = DEGenesClean3(Any09,:);
Any09Table = sortrows(Any09Table,'Log2FC09','descend');

Any18 = char_bin_code(:,2) == '1';
Any18Table = DEGenesClean3(Any18,:);
Any18Table = sortrows(Any18Table,'Log2FC18','descend');

Any27 = char_bin_code(:,3) == '1';
Any27Table = DEGenesClean3(Any27,:);
Any27Table = sortrows(Any27Table,'Log2FC27','descend');

Any36 = char_bin_code(:,4) == '1';
Any36Table = DEGenesClean3(Any36,:);
Any36Table = sortrows(Any36Table,'Log2FC36','descend');

Any45 = char_bin_code(:,5) == '1';
Any45Table = DEGenesClean3(Any45,:);
Any45Table = sortrows(Any45Table,'Log2FC45','descend');

Any54 = char_bin_code(:,6) == '1';
Any54Table = DEGenesClean3(Any54,:);
Any54Table = sortrows(Any54Table,'Log2FC54','descend');

Any63 = char_bin_code(:,7) == '1';
Any63Table = DEGenesClean3(Any63,:);
Any63Table = sortrows(Any63Table,'Log2FC63','descend');

Any72 = char_bin_code(:,8) == '1';
Any72Table = DEGenesClean3(Any72,:);
Any72Table = sortrows(Any72Table,'Log2FC72','descend');

Any81 = char_bin_code(:,9) == '1';
Any81Table = DEGenesClean3(Any81,:);
Any81Table = sortrows(Any81Table,'Log2FC81','descend');

Any90 = char_bin_code(:,10) == '1';
Any90Table = DEGenesClean3(Any90,:);
Any90Table = sortrows(Any90Table,'Log2FC90','descend');

% Uncomment to save things to file
%{
% Write DE Genes at each frequency to an Excel file
table_file = 'C:\Users\Tom\Documents\iCloudDrive\Documents\My Papers\p38\Data\RNAseq\GENE_TABLES\DE_Genes_Tables005.xlsx';

writetable(Any09Table, table_file, 'Sheet', '0.09 hr-1');
writetable(Any18Table, table_file, 'Sheet', '0.18 hr-1');
writetable(Any27Table, table_file, 'Sheet', '0.27 hr-1');
writetable(Any36Table, table_file, 'Sheet', '0.36 hr-1');
writetable(Any45Table, table_file, 'Sheet', '0.45 hr-1');
writetable(Any54Table, table_file, 'Sheet', '0.54 hr-1');
writetable(Any63Table, table_file, 'Sheet', '0.63 hr-1');
writetable(Any72Table, table_file, 'Sheet', '0.72 hr-1');
writetable(Any81Table, table_file, 'Sheet', '0.81 hr-1');
writetable(Any90Table, table_file, 'Sheet', '0.90 hr-1');
%}


%% Clustergram of Unique Genes

DEGenes4 = innerjoin(DEC20Table2, DEC40Table2);
DEGenes4 = innerjoin(DEGenes4, DE36SBTable2);
DEGenes4 = innerjoin(DEGenes4, DE09Table2);
DEGenes4 = innerjoin(DEGenes4, DE18Table2);
DEGenes4 = innerjoin(DEGenes4, DE27Table2);
DEGenes4 = innerjoin(DEGenes4, DE36Table2);
DEGenes4 = innerjoin(DEGenes4, DE45Table2);
DEGenes4 = innerjoin(DEGenes4, DE54Table2);
DEGenes4 = innerjoin(DEGenes4, DE63Table2);
DEGenes4 = innerjoin(DEGenes4, DE72Table2);
DEGenes4 = innerjoin(DEGenes4, DE81Table2);
DEGenes4 = innerjoin(DEGenes4, DE90Table2);

DEGenes4 = addvars(DEGenes4, gene_names, 'After', 1);
DEGenesClean4 = rmmissing(DEGenes4);

sig_cutoff = 0.01;
FC_cutoff = log2(2);
FC_cutoff2 = log2(2);

sigSB = zeros(size(DEGenesClean4,1),1);
sig09 = zeros(size(DEGenesClean4,1),1);
sig18 = zeros(size(DEGenesClean4,1),1);
sig27 = zeros(size(DEGenesClean4,1),1);
sig36 = zeros(size(DEGenesClean4,1),1);
sig45 = zeros(size(DEGenesClean4,1),1);
sig54 = zeros(size(DEGenesClean4,1),1);
sig63 = zeros(size(DEGenesClean4,1),1);
sig72 = zeros(size(DEGenesClean4,1),1);
sig81 = zeros(size(DEGenesClean4,1),1);
sig90 = zeros(size(DEGenesClean4,1),1);


for i=1:size(DEGenesClean4,1)

    if ((DEGenesClean4.Padj36SB(i) < sig_cutoff) & (abs(DEGenesClean4.Log2FC36SB(i)) > FC_cutoff))
        sigSB(i) = 1;
    end
    if ((DEGenesClean4.Padj09(i) < sig_cutoff) & (abs(DEGenesClean4.Log2FC09(i)) > FC_cutoff))
        sig09(i) = 1;
    end
    if ((DEGenesClean4.Padj18(i) < sig_cutoff) & (abs(DEGenesClean4.Log2FC18(i)) > FC_cutoff))
        sig18(i) = 1;
    end
    if ((DEGenesClean4.Padj27(i) < sig_cutoff) & (abs(DEGenesClean4.Log2FC27(i)) > FC_cutoff))
        sig27(i) = 1;
    end
    if ((DEGenesClean4.Padj36(i) < sig_cutoff) & (abs(DEGenesClean4.Log2FC36(i)) > FC_cutoff))
        sig36(i) = 1;
    end
    if ((DEGenesClean4.Padj45(i) < sig_cutoff) & (abs(DEGenesClean4.Log2FC45(i)) > FC_cutoff))
        sig45(i) = 1;
    end
    if ((DEGenesClean4.Padj54(i) < sig_cutoff) & (abs(DEGenesClean4.Log2FC54(i)) > FC_cutoff))
        sig54(i) = 1;
    end
    if ((DEGenesClean4.Padj63(i) < sig_cutoff) & (abs(DEGenesClean4.Log2FC63(i)) > FC_cutoff))
        sig63(i) = 1;
    end
    if ((DEGenesClean4.Padj72(i) < sig_cutoff) & (abs(DEGenesClean4.Log2FC72(i)) > FC_cutoff))
        sig72(i) = 1;
    end
    if ((DEGenesClean4.Padj81(i) < sig_cutoff) & (abs(DEGenesClean4.Log2FC81(i)) > FC_cutoff))
        sig81(i) = 1;
    end
    if ((DEGenesClean4.Padj90(i) < sig_cutoff) & (abs(DEGenesClean4.Log2FC90(i)) > FC_cutoff))
        sig90(i) = 1;
    end

    bin_code2{i} = sprintf('%d', [sigSB(i), sig09(i), sig18(i), sig27(i), sig36(i), sig45(i), sig54(i),...
        sig63(i), sig72(i), sig81(i), sig90(i)]);

end


bin_code2 = string(bin_code2);
DEGenesClean4 = addvars(DEGenesClean4, bin_code2', 'After', 2);

DEGenesClean4 = renamevars(DEGenesClean4, 'Var3', 'bin_code');

OnlySB_2 = DEGenesClean4.bin_code == '10000000000';
Only09_2 = DEGenesClean4.bin_code == '01000000000';
Only18_2 = DEGenesClean4.bin_code == '00100000000';
Only27_2 = DEGenesClean4.bin_code == '00010000000';
Only36_2 = DEGenesClean4.bin_code == '00001000000';
Only45_2 = DEGenesClean4.bin_code == '00000100000';
Only54_2 = DEGenesClean4.bin_code == '00000010000';
Only63_2 = DEGenesClean4.bin_code == '00000001000';
Only72_2 = DEGenesClean4.bin_code == '00000000100';
Only81_2 = DEGenesClean4.bin_code == '00000000010';
Only90_2 = DEGenesClean4.bin_code == '00000000001';


OnlySBTable2 = DEGenesClean4(OnlySB_2,:);
OnlySBTable2 = sortrows(OnlySBTable2,'Log2FC36SB','descend');
Only09Table2 = DEGenesClean4(Only09_2,:);
Only09Table2 = sortrows(Only09Table2,'Log2FC09','descend');
Only18Table2 = DEGenesClean4(Only18_2,:);
Only18Table2 = sortrows(Only18Table2,'Log2FC18','descend');
Only27Table2 = DEGenesClean4(Only27_2,:);
Only27Table2 = sortrows(Only27Table2,'Log2FC27','descend');
Only36Table2 = DEGenesClean4(Only36_2,:);
Only36Table2 = sortrows(Only36Table2,'Log2FC36','descend');
Only45Table2 = DEGenesClean4(Only45_2,:);
Only45Table2 = sortrows(Only45Table2,'Log2FC45','descend');
Only54Table2 = DEGenesClean4(Only54_2,:);
Only54Table2 = sortrows(Only54Table2,'Log2FC54','descend');
Only63Table2 = DEGenesClean4(Only63_2,:);
Only63Table2 = sortrows(Only63Table2,'Log2FC63','descend');
Only72Table2 = DEGenesClean4(Only72_2,:);
Only72Table2 = sortrows(Only72Table2,'Log2FC72','descend');
Only81Table2 = DEGenesClean4(Only81_2,:);
Only81Table2 = sortrows(Only81Table2,'Log2FC81','descend');
Only90Table2 = DEGenesClean4(Only90_2,:);
Only90Table2 = sortrows(Only90Table2,'Log2FC90','descend');



DEGenes5 = [Only09Table2; Only18Table2];
DEGenes5 = [DEGenes5; Only27Table2];
DEGenes5 = [DEGenes5; Only36Table2];
DEGenes5 = [DEGenes5; Only45Table2];
DEGenes5 = [DEGenes5; Only54Table2];
DEGenes5 = [DEGenes5; Only63Table2];
DEGenes5 = [DEGenes5; Only72Table2];
DEGenes5 = [DEGenes5; Only81Table2];
DEGenes5 = [DEGenes5; Only90Table2];


cg3 = clustergram([DEGenes5{:,"negMeanC20"} ...
    DEGenes5{:,"expMeanC20"}...
    DEGenes5{:,"expMeanC40"}...
    DEGenes5{:,"expMean36SB"}...
    DEGenes5{:,"expMean09"}...
    DEGenes5{:,"expMean18"}...
    DEGenes5{:,"expMean27"}...
    DEGenes5{:,"expMean36"}...
    DEGenes5{:,"expMean45"}...
    DEGenes5{:,"expMean54"}...
    DEGenes5{:,"expMean63"}...
    DEGenes5{:,"expMean72"}...
    DEGenes5{:,"expMean81"}...
    DEGenes5{:,"expMean90"}...
    ], 'Standardize', 'row');

cg3.RowPDist = 'correlation';
%cg3.OptimalLeafOrder = true;
cg3.Colormap = redbluecmap;
cg3.DisplayRange = 4;
cgAxes = plot(cg3);
set(cgAxes, 'Clim', [-3.5 3.5])
colormap hot


% Export data from clustergram to clusterMatrix before doing this!!
for i=1:size(clusterMatrix.ExprValues,1)
    norm_heat_map(i,:) = clusterMatrix.ExprValues(i,:)./mean(clusterMatrix.ExprValues(i,:));
end


heat_map = clusterMatrix.ExprValues;
ordered_heat_map = [norm_heat_map(:,2) norm_heat_map(:,3) norm_heat_map(:,1) norm_heat_map(:,4) norm_heat_map(:,12) norm_heat_map(:,11) norm_heat_map(:,10) norm_heat_map(:,14) norm_heat_map(:,5) norm_heat_map(:,13) norm_heat_map(:,6) norm_heat_map(:,9) norm_heat_map(:,7) norm_heat_map(:,8)];


labs = cell(size(ordered_heat_map,1),1);
labs(:) = {''};
figure
hm = heatmap(ordered_heat_map);
hm.Colormap = CustomColormap;
hm.YDisplayLabels = labs;
grid(hm, 'off')
cmin = -0.05;
cmax = 2.6;
clim([cmin cmax]);



labs = cell(size(ordered_heat_map,1),1);
labs(:) = {''};
figure
hm = heatmap(log2(ordered_heat_map));
hm.Colormap = CustomColormap;
hm.YDisplayLabels = labs;
grid(hm, 'off')
cmin = -2.25;
cmax = 2.0;
clim([cmin cmax]);


figure
plot(norm_heat_map(:,1))
hold on
plot(norm_heat_map(:,2))
plot(norm_heat_map(:,3))
plot(norm_heat_map(:,4))
plot(norm_heat_map(:,5))
plot(norm_heat_map(:,6))
plot(norm_heat_map(:,7))
plot(norm_heat_map(:,8))
plot(norm_heat_map(:,9))
plot(norm_heat_map(:,10))
plot(norm_heat_map(:,11))
plot(norm_heat_map(:,12))
plot(norm_heat_map(:,13))
plot(norm_heat_map(:,14))




%% PCA
% ALL GENES
PCA_Table = removevars(DEGenes4, [1 2 5 6 7 8 10 11 12 13 15 16 17 18 20 21 22 23 25 26 27 28 30 31 32 33 35 36 37 38 40 41 42 43 45 46 47 48 50 51 52 53 55 56 57 58 60 61 62 63 65 66 67]);
PCA_Matrix = table2array(PCA_Table);


PCA_Table2 = removevars(DEGenes4, [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 20 21 22 23 25 26 27 28 30 31 32 33 35 36 37 38 40 41 42 43 45 46 47 48 50 51 52 53 55 56 57 58 60 61 62 63 65 66 67]);
PCA_Matrix2 = table2array(PCA_Table2);

[coeff, score, latent, tsquared, explained] = pca(PCA_Matrix);
[coeff2, score2, latent2, tsquared2, explained2] = pca(PCA_Matrix2);


figyz = figure;
figyz.Theme = "light";
scatter(coeff(1:3,2), coeff(1:3,3), 100, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0 0 1]);
hold on
scatter(coeff(4,2), coeff(4,3), 100, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.5 0.5 0.5]);
scatter(coeff(5:14,2), coeff(5:14,3), 100, 'MarkerEdgeColor','k', 'MarkerFaceColor',[1 0 0]);
box on
grid on
ax = gca;
ax.LineWidth = 2;
ax.GridLineWidth = 1;
ax.FontSize = 13;
ax.XDir = 'reverse';
xlabel('Second Principal Component')
ylabel('Third Principal Component')
%axis equal
xlim([-0.3 0.75])
ylim([-0.45 0.5])



figxy = figure;
figxy.Theme = "light";
scatter(coeff(1:3,1), coeff(1:3,2), 100, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0 0 1]);
hold on
scatter(coeff(4,1), coeff(4,2), 100, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.5 0.5 0.5]);
scatter(coeff(5:14,1), coeff(5:14,2), 100, 'MarkerEdgeColor','k', 'MarkerFaceColor',[1 0 0]);
box on
grid on
ax = gca;
ax.LineWidth = 2;
ax.GridLineWidth = 1;
ax.FontSize = 13;
xlabel('First Principal Component')
ylabel('Second Principal Component')
%axis equal
xlim([0.2 0.38])
ylim([-0.3 0.75])


figxz = figure;
figxz.Theme = "light";
scatter(coeff(1:3,1), coeff(1:3,3), 100, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0 0 1]);
hold on
scatter(coeff(4,1), coeff(4,3), 100, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.5 0.5 0.5]);
scatter(coeff(5:14,1), coeff(5:14,3), 100, 'MarkerEdgeColor','k', 'MarkerFaceColor',[1 0 0]);
box on
grid on
ax = gca;
ax.LineWidth = 2;
ax.GridLineWidth = 1;
ax.FontSize = 13;
xlabel('First Principal Component')
ylabel('Third Principal Component')
%axis equal
xlim([0.2 0.38])
ylim([-0.45 0.5])


marksize = 300;
fig = figure;
fig.Theme = "light";
bx = scatter3(0.38, coeff(1:3,2), coeff(1:3,3), marksize, 'MarkerEdgeColor',[173 216 230]./256, 'MarkerFaceColor',[173 216 230]./256);
hold on
by = scatter3(coeff(1:3,1), 0.75, coeff(1:3,3), marksize, 'MarkerEdgeColor',[173 216 230]./256, 'MarkerFaceColor',[173 216 230]./256);
bz = scatter3(coeff(1:3,1), coeff(1:3,2), -0.45, marksize, 'MarkerEdgeColor',[173 216 230]./256, 'MarkerFaceColor',[173 216 230]./256);
alpha(bx,0.5);
alpha(by,0.5);
alpha(bz,0.5);
for i=1:3
    plot3([coeff(i,1) 0.4], [coeff(i,2), coeff(i,2)], [coeff(i,3) coeff(i,3)], 'Color', [173 216 230]./256)
    plot3([coeff(i,1), coeff(i,1)], [coeff(i,2), 0.75], [coeff(i,3) coeff(i,3)], 'Color', [173 216 230]./256)
    plot3([coeff(i,1), coeff(i,1)], [coeff(i,2), coeff(i,2)], [-0.45 coeff(i,3)], 'Color', [173 216 230]./256)
end

rx = scatter3(0.38, coeff(5:14,2), coeff(5:14,3), marksize, 'MarkerEdgeColor',[255 204 203]./256, 'MarkerFaceColor',[255 204 203]./256);
ry = scatter3(coeff(5:14,1), 0.75, coeff(5:14,3), marksize, 'MarkerEdgeColor',[255 204 203]./256, 'MarkerFaceColor',[255 204 203]./256);
rz = scatter3(coeff(5:14,1), coeff(5:14,2), -0.45, marksize, 'MarkerEdgeColor',[255 204 203]./256, 'MarkerFaceColor',[255 204 203]./256);
alpha(rx,0.5);
alpha(ry,0.5);
alpha(rz,0.5);
for i=5:14
    rlx = plot3([coeff(i,1) 0.4], [coeff(i,2), coeff(i,2)], [coeff(i,3) coeff(i,3)], 'Color', [255 204 203]./256);
    rly = plot3([coeff(i,1), coeff(i,1)], [coeff(i,2), 0.75], [coeff(i,3) coeff(i,3)], 'Color', [255 204 203]./256);
    rlz = plot3([coeff(i,1), coeff(i,1)], [coeff(i,2), coeff(i,2)], [-0.45 coeff(i,3)], 'Color', [255 204 203]./256);
end

gx = scatter3(0.38, coeff(4,2), coeff(4,3), marksize, 'MarkerEdgeColor',[0.9 0.9 0.9], 'MarkerFaceColor',[0.9 0.9 0.9]);
gy = scatter3(coeff(4,1), 0.75, coeff(4,3), marksize, 'MarkerEdgeColor',[0.9 0.9 0.9], 'MarkerFaceColor',[0.9 0.9 0.9]);
gz = scatter3(coeff(4,1), coeff(4,2), -0.45, marksize, 'MarkerEdgeColor',[0.9 0.9 0.9], 'MarkerFaceColor',[0.9 0.9 0.9]);
alpha(gx,0.5);
alpha(gy,0.5);
alpha(gz,0.5);
glx = plot3([coeff(4,1) 0.4], [coeff(4,2), coeff(4,2)], [coeff(4,3) coeff(4,3)], 'Color', [0.925 0.925 0.925]);
gly = plot3([coeff(4,1), coeff(4,1)], [coeff(4,2), 0.75], [coeff(4,3) coeff(4,3)], 'Color', [0.925 0.925 0.925]);
glz = plot3([coeff(4,1), coeff(4,1)], [coeff(4,2), coeff(4,2)], [-0.45 coeff(4,3)], 'Color', [0.925 0.925 0.925]);

scatter3(coeff(1:3,1), coeff(1:3,2), coeff(1:3,3), marksize, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0 0 1]);
scatter3(coeff(4,1), coeff(4,2), coeff(4,3), marksize, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.5 0.5 0.5]);
scatter3(coeff(5:14,1), coeff(5:14,2), coeff(5:14,3), marksize, 'MarkerEdgeColor','k', 'MarkerFaceColor',[1 0 0]);
box on
grid on
ax = gca;
ax.LineWidth = 2;
ax.GridLineWidth = 1;
ax.FontSize = 13;
xlabel('First Principal Component')
ylabel('Second Principal Component')
zlabel('Third Principal Component')
%axis equal
xlim([0.2 0.38])
ylim([-0.3 0.75])
zlim([-0.45 0.5])







figxy = figure;
figxy.Theme = "light";
scatter(coeff2(1:10,1), coeff2(1:10,2), 100, 'MarkerEdgeColor','k', 'MarkerFaceColor',[1 0 0]);
hold on
box on
grid on
ax = gca;
ax.LineWidth = 2;
ax.GridLineWidth = 1;
ax.FontSize = 13;
ax.YDir = 'reverse';
xlabel('First Principal Component')
ylabel('Second Principal Component')
%axis equal
xlim([0.28 0.36])
ylim([-0.42 0.6])



%% Gene plots

DEGenesClean4 = DEGenesClean3;
DEGenesClean4p = DEGenesClean3;
DEGenesClean4x = DEGenesClean3;
DEGenesClean4.Properties.RowNames = DEGenesClean4{:,1};
DEGenesClean4p.Properties.RowNames = DEGenesClean4{:,1};
DEGenesClean4x.Properties.RowNames = DEGenesClean4{:,1};
exp_columns = [21 26 31 36 41 46 51 56 61 66];
frequencies = [0.09 0.18 0.27 0.36 0.45 0.54 0.63 0.72 0.81 0.90];

DEGenesClean4 = removevars(DEGenesClean4, [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22 23 24 25 27 28 29 30 32 33 34 35 37 38 39 40 42 43 44 45 47 48 49 50 52 53 54 55 57 58 59 60 62 63 64 65 67 68]);
DEGenesClean4p = removevars(DEGenesClean4p, [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 24 25 26 27 29 30 31 32 34 35 36 37 39 40 41 42 44 45 46 47 49 50 51 52 54 55 56 57 59 60 61 62 64 65 66 67]);
DEGenesClean4x = removevars(DEGenesClean4x, [2 3 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 21 22 23 24 26 27 28 29 31 32 33 34 36 37 38 39 41 42 43 44 46 47 48 49 51 52 53 54 56 57 58 59 61 62 63 64 66 67 68]);
DEGenesClean4j = join(DEGenesClean4, DEGenesClean4p);
DEGenesClean4j = removevars(DEGenesClean4j, 1);


% Load Enrichment Data
% PC File
enrichment_file = 'C:\Users\Tom\Documents\iCloudDrive\Documents\My Papers\p38\Data\RNAseq\GENE_TABLES\Enrichment_table_min.xlsx';

% Mac File
%enrichment_file = '/Users/thomaskuhlman/Documents/My Papers/p38/Data/RNAseq/GENE_TABLES/Enrichment_table_min.xlsx';

enrichment_09 = readtable(enrichment_file, 'Sheet', '0.09 hr-1 padj=0.005');
enrichment_18 = readtable(enrichment_file, 'Sheet', '0.18 hr-1 padj=0.005');
enrichment_27 = readtable(enrichment_file, 'Sheet', '0.27 hr-1 padj=0.005');
enrichment_36 = readtable(enrichment_file, 'Sheet', '0.36 hr-1 padj=0.005');
enrichment_45 = readtable(enrichment_file, 'Sheet', '0.45 hr-1 padj=0.005');
enrichment_54 = readtable(enrichment_file, 'Sheet', '0.54 hr-1 padj=0.005');
enrichment_63 = readtable(enrichment_file, 'Sheet', '0.63 hr-1 padj=0.01');
enrichment_72 = readtable(enrichment_file, 'Sheet', '0.72 hr-1 padj=0.005');
enrichment_81 = readtable(enrichment_file, 'Sheet', '0.81 hr-1 padj=0.005');
enrichment_90 = readtable(enrichment_file, 'Sheet', '0.90 hr-1 padj=0.005');


stringsToFind = {'GO:MF','GO:BP','GO:CC','KEGG','REAC','WP','TF','MIRNA','HPA','CORUM','HP'};
excludeStrings = {'REACTOME root term', 'KEGG root term', 'WIKIPATHWAYS','HPA root','CORUM root','HP root'};


sources_09 = enrichment_09(ismember(enrichment_09.source, stringsToFind),:);
enrich_09 = sources_09(~ismember(sources_09.term_name, excludeStrings),:);
radio_09 = table('Size', [1 22], 'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'source', 'term', 'Log2FC09', 'Log2FC18', 'Log2FC27', 'Log2FC36', 'Log2FC45', 'Log2FC54', 'Log2FC63', 'Log2FC72', 'Log2FC81', 'Log2FC90', 'Padj09', 'Padj18', 'Padj27', 'Padj36', 'Padj45', 'Padj54', 'Padj63', 'Padj72', 'Padj81', 'Padj90'});
for i=1:size(enrich_09,1)
    temp_FC = DEGenesClean4j(split(enrich_09.intersections{i},','),:);
    temp_FC.Properties.RowNames = {};

    source = enrich_09.source{i};
    for j=1:size(temp_FC,1)-1
        source = [source; enrich_09.source{i}];
    end
    temp_FC2 = addvars(temp_FC, string(source), 'Before', 'Log2FC09');
    temp_FC2 = renamevars(temp_FC2, 'Var1', 'source');

    term = enrich_09.term_name{i};
    for j=1:size(temp_FC,1)-1
        term = [term; enrich_09.term_name{i}];
    end
    temp_FC3 = addvars(temp_FC2, string(term), 'Before', 'Log2FC09');
    temp_FC3 = renamevars(temp_FC3, 'Var2', 'term');

    radio_09 = [radio_09; temp_FC3];
end
radio_09(1,:) = [];
repeats09a = groupsummary(radio_09, ["Log2FC09","Log2FC18","Log2FC27","Log2FC36","Log2FC45","Log2FC54","Log2FC63","Log2FC72","Log2FC81","Log2FC90","Padj09","Padj18","Padj27","Padj36","Padj45","Padj54","Padj63","Padj72","Padj81","Padj90"]);

sources_18 = enrichment_18(ismember(enrichment_18.source, stringsToFind),:);
enrich_18 = sources_18(~ismember(sources_18.term_name, excludeStrings),:);
radio_18 = table('Size', [1 22], 'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'source', 'term', 'Log2FC09', 'Log2FC18', 'Log2FC27', 'Log2FC36', 'Log2FC45', 'Log2FC54', 'Log2FC63', 'Log2FC72', 'Log2FC81', 'Log2FC90', 'Padj09', 'Padj18', 'Padj27', 'Padj36', 'Padj45', 'Padj54', 'Padj63', 'Padj72', 'Padj81', 'Padj90'});
for i=1:size(enrich_18,1)
    temp_FC = DEGenesClean4j(split(enrich_18.intersections{i},','),:);
    temp_FC.Properties.RowNames = {};

    source = enrich_18.source{i};
    for j=1:size(temp_FC,1)-1
        source = [source; enrich_18.source{i}];
    end
    temp_FC2 = addvars(temp_FC, string(source), 'Before', 'Log2FC09');
    temp_FC2 = renamevars(temp_FC2, 'Var1', 'source');

    term = enrich_18.term_name{i};
    for j=1:size(temp_FC,1)-1
        term = [term; enrich_18.term_name{i}];
    end
    temp_FC3 = addvars(temp_FC2, string(term), 'Before', 'Log2FC09');
    temp_FC3 = renamevars(temp_FC3, 'Var2', 'term');

    radio_18 = [radio_18; temp_FC3];
end
radio_18(1,:) = [];
repeats18a = groupsummary(radio_18, ["Log2FC09","Log2FC18","Log2FC27","Log2FC36","Log2FC45","Log2FC54","Log2FC63","Log2FC72","Log2FC81","Log2FC90","Padj09","Padj18","Padj27","Padj36","Padj45","Padj54","Padj63","Padj72","Padj81","Padj90"]);

sources_27 = enrichment_27(ismember(enrichment_27.source, stringsToFind),:);
enrich_27 = sources_27(~ismember(sources_27.term_name, excludeStrings),:);
radio_27 = table('Size', [1 22], 'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'source', 'term', 'Log2FC09', 'Log2FC18', 'Log2FC27', 'Log2FC36', 'Log2FC45', 'Log2FC54', 'Log2FC63', 'Log2FC72', 'Log2FC81', 'Log2FC90', 'Padj09', 'Padj18', 'Padj27', 'Padj36', 'Padj45', 'Padj54', 'Padj63', 'Padj72', 'Padj81', 'Padj90'});
for i=1:size(enrich_27,1)
    temp_FC = DEGenesClean4j(split(enrich_27.intersections{i},','),:);
    temp_FC.Properties.RowNames = {};

    source = enrich_27.source{i};
    for j=1:size(temp_FC,1)-1
        source = [source; enrich_27.source{i}];
    end
    temp_FC2 = addvars(temp_FC, string(source), 'Before', 'Log2FC09');
    temp_FC2 = renamevars(temp_FC2, 'Var1', 'source');

    term = enrich_27.term_name{i};
    for j=1:size(temp_FC,1)-1
        term = [term; enrich_27.term_name{i}];
    end
    temp_FC3 = addvars(temp_FC2, string(term), 'Before', 'Log2FC09');
    temp_FC3 = renamevars(temp_FC3, 'Var2', 'term');

    radio_27 = [radio_27; temp_FC3];
end
radio_27(1,:) = [];
repeats27a = groupsummary(radio_27, ["Log2FC09","Log2FC18","Log2FC27","Log2FC36","Log2FC45","Log2FC54","Log2FC63","Log2FC72","Log2FC81","Log2FC90","Padj09","Padj18","Padj27","Padj36","Padj45","Padj54","Padj63","Padj72","Padj81","Padj90"]);

sources_36 = enrichment_36(ismember(enrichment_36.source, stringsToFind),:);
enrich_36 = sources_36(~ismember(sources_36.term_name, excludeStrings),:);
radio_36 = table('Size', [1 22], 'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'source', 'term', 'Log2FC09', 'Log2FC18', 'Log2FC27', 'Log2FC36', 'Log2FC45', 'Log2FC54', 'Log2FC63', 'Log2FC72', 'Log2FC81', 'Log2FC90', 'Padj09', 'Padj18', 'Padj27', 'Padj36', 'Padj45', 'Padj54', 'Padj63', 'Padj72', 'Padj81', 'Padj90'});
for i=1:size(enrich_36,1)
    temp_FC = DEGenesClean4j(split(enrich_36.intersections{i},','),:);
    temp_FC.Properties.RowNames = {};

    source = enrich_36.source{i};
    for j=1:size(temp_FC,1)-1
        source = [source; enrich_36.source{i}];
    end
    temp_FC2 = addvars(temp_FC, string(source), 'Before', 'Log2FC09');
    temp_FC2 = renamevars(temp_FC2, 'Var1', 'source');

    term = enrich_36.term_name{i};
    for j=1:size(temp_FC,1)-1
        term = [term; enrich_36.term_name{i}];
    end
    temp_FC3 = addvars(temp_FC2, string(term), 'Before', 'Log2FC09');
    temp_FC3 = renamevars(temp_FC3, 'Var2', 'term');

    radio_36 = [radio_36; temp_FC3];
end
radio_36(1,:) = [];
repeats36a = groupsummary(radio_36, ["Log2FC09","Log2FC18","Log2FC27","Log2FC36","Log2FC45","Log2FC54","Log2FC63","Log2FC72","Log2FC81","Log2FC90","Padj09","Padj18","Padj27","Padj36","Padj45","Padj54","Padj63","Padj72","Padj81","Padj90"]);

sources_45 = enrichment_45(ismember(enrichment_45.source, stringsToFind),:);
enrich_45 = sources_45(~ismember(sources_45.term_name, excludeStrings),:);
radio_45 = table('Size', [1 22], 'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'source', 'term', 'Log2FC09', 'Log2FC18', 'Log2FC27', 'Log2FC36', 'Log2FC45', 'Log2FC54', 'Log2FC63', 'Log2FC72', 'Log2FC81', 'Log2FC90', 'Padj09', 'Padj18', 'Padj27', 'Padj36', 'Padj45', 'Padj54', 'Padj63', 'Padj72', 'Padj81', 'Padj90'});
for i=1:size(enrich_45,1)
    temp_FC = DEGenesClean4j(split(enrich_45.intersections{i},','),:);
    temp_FC.Properties.RowNames = {};

    source = enrich_45.source{i};
    for j=1:size(temp_FC,1)-1
        source = [source; enrich_45.source{i}];
    end
    temp_FC2 = addvars(temp_FC, string(source), 'Before', 'Log2FC09');
    temp_FC2 = renamevars(temp_FC2, 'Var1', 'source');

    term = enrich_45.term_name{i};
    for j=1:size(temp_FC,1)-1
        term = [term; enrich_45.term_name{i}];
    end
    temp_FC3 = addvars(temp_FC2, string(term), 'Before', 'Log2FC09');
    temp_FC3 = renamevars(temp_FC3, 'Var2', 'term');

    radio_45 = [radio_45; temp_FC3];
end
radio_45(1,:) = [];
repeats45a = groupsummary(radio_45, ["Log2FC09","Log2FC18","Log2FC27","Log2FC36","Log2FC45","Log2FC54","Log2FC63","Log2FC72","Log2FC81","Log2FC90","Padj09","Padj18","Padj27","Padj36","Padj45","Padj54","Padj63","Padj72","Padj81","Padj90"]);

sources_54 = enrichment_54(ismember(enrichment_54.source, stringsToFind),:);
enrich_54 = sources_54(~ismember(sources_54.term_name, excludeStrings),:);
radio_54 = table('Size', [1 22], 'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'source', 'term', 'Log2FC09', 'Log2FC18', 'Log2FC27', 'Log2FC36', 'Log2FC45', 'Log2FC54', 'Log2FC63', 'Log2FC72', 'Log2FC81', 'Log2FC90', 'Padj09', 'Padj18', 'Padj27', 'Padj36', 'Padj45', 'Padj54', 'Padj63', 'Padj72', 'Padj81', 'Padj90'});
for i=1:size(enrich_54,1)
    temp_FC = DEGenesClean4j(split(enrich_54.intersections{i},','),:);
    temp_FC.Properties.RowNames = {};

    source = enrich_54.source{i};
    for j=1:size(temp_FC,1)-1
        source = [source; enrich_54.source{i}];
    end
    temp_FC2 = addvars(temp_FC, string(source), 'Before', 'Log2FC09');
    temp_FC2 = renamevars(temp_FC2, 'Var1', 'source');

    term = enrich_54.term_name{i};
    for j=1:size(temp_FC,1)-1
        term = [term; enrich_54.term_name{i}];
    end
    temp_FC3 = addvars(temp_FC2, string(term), 'Before', 'Log2FC09');
    temp_FC3 = renamevars(temp_FC3, 'Var2', 'term');

    radio_54 = [radio_54; temp_FC3];
end
radio_54(1,:) = [];
repeats54a = groupsummary(radio_54, ["Log2FC09","Log2FC18","Log2FC27","Log2FC36","Log2FC45","Log2FC54","Log2FC63","Log2FC72","Log2FC81","Log2FC90","Padj09","Padj18","Padj27","Padj36","Padj45","Padj54","Padj63","Padj72","Padj81","Padj90"]);

sources_63 = enrichment_63(ismember(enrichment_63.source, stringsToFind),:);
enrich_63 = sources_63(~ismember(sources_63.term_name, excludeStrings),:);
radio_63 = table('Size', [1 22], 'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'source', 'term', 'Log2FC09', 'Log2FC18', 'Log2FC27', 'Log2FC36', 'Log2FC45', 'Log2FC54', 'Log2FC63', 'Log2FC72', 'Log2FC81', 'Log2FC90', 'Padj09', 'Padj18', 'Padj27', 'Padj36', 'Padj45', 'Padj54', 'Padj63', 'Padj72', 'Padj81', 'Padj90'});
for i=1:size(enrich_63,1)
    temp_FC = DEGenesClean4j(split(enrich_63.intersections{i},','),:);
    temp_FC.Properties.RowNames = {};

    source = enrich_63.source{i};
    for j=1:size(temp_FC,1)-1
        source = [source; enrich_63.source{i}];
    end
    temp_FC2 = addvars(temp_FC, string(source), 'Before', 'Log2FC09');
    temp_FC2 = renamevars(temp_FC2, 'Var1', 'source');

    term = enrich_63.term_name{i};
    for j=1:size(temp_FC,1)-1
        term = [term; enrich_63.term_name{i}];
    end
    temp_FC3 = addvars(temp_FC2, string(term), 'Before', 'Log2FC09');
    temp_FC3 = renamevars(temp_FC3, 'Var2', 'term');

    radio_63 = [radio_63; temp_FC3];
end
radio_63(1,:) = [];
repeats63a = groupsummary(radio_63, ["Log2FC09","Log2FC18","Log2FC27","Log2FC36","Log2FC45","Log2FC54","Log2FC63","Log2FC72","Log2FC81","Log2FC90","Padj09","Padj18","Padj27","Padj36","Padj45","Padj54","Padj63","Padj72","Padj81","Padj90"]);

sources_72 = enrichment_72(ismember(enrichment_72.source, stringsToFind),:);
enrich_72 = sources_72(~ismember(sources_72.term_name, excludeStrings),:);
radio_72 = table('Size', [1 22], 'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'source', 'term', 'Log2FC09', 'Log2FC18', 'Log2FC27', 'Log2FC36', 'Log2FC45', 'Log2FC54', 'Log2FC63', 'Log2FC72', 'Log2FC81', 'Log2FC90', 'Padj09', 'Padj18', 'Padj27', 'Padj36', 'Padj45', 'Padj54', 'Padj63', 'Padj72', 'Padj81', 'Padj90'});
for i=1:size(enrich_72,1)
    temp_FC = DEGenesClean4j(split(enrich_72.intersections{i},','),:);
    temp_FC.Properties.RowNames = {};

    source = enrich_72.source{i};
    for j=1:size(temp_FC,1)-1
        source = [source; enrich_72.source{i}];
    end
    temp_FC2 = addvars(temp_FC, string(source), 'Before', 'Log2FC09');
    temp_FC2 = renamevars(temp_FC2, 'Var1', 'source');

    term = enrich_72.term_name{i};
    for j=1:size(temp_FC,1)-1
        term = [term; enrich_72.term_name{i}];
    end
    temp_FC3 = addvars(temp_FC2, string(term), 'Before', 'Log2FC09');
    temp_FC3 = renamevars(temp_FC3, 'Var2', 'term');

    radio_72 = [radio_72; temp_FC3];
end
radio_72(1,:) = [];
repeats72a = groupsummary(radio_72, ["Log2FC09","Log2FC18","Log2FC27","Log2FC36","Log2FC45","Log2FC54","Log2FC63","Log2FC72","Log2FC81","Log2FC90","Padj09","Padj18","Padj27","Padj36","Padj45","Padj54","Padj63","Padj72","Padj81","Padj90"]);

sources_81 = enrichment_81(ismember(enrichment_81.source, stringsToFind),:);
enrich_81 = sources_81(~ismember(sources_81.term_name, excludeStrings),:);
radio_81 = table('Size', [1 22], 'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'source', 'term', 'Log2FC09', 'Log2FC18', 'Log2FC27', 'Log2FC36', 'Log2FC45', 'Log2FC54', 'Log2FC63', 'Log2FC72', 'Log2FC81', 'Log2FC90', 'Padj09', 'Padj18', 'Padj27', 'Padj36', 'Padj45', 'Padj54', 'Padj63', 'Padj72', 'Padj81', 'Padj90'});
for i=1:size(enrich_81,1)
    temp_FC = DEGenesClean4j(split(enrich_81.intersections{i},','),:);
    temp_FC.Properties.RowNames = {};

    source = enrich_81.source{i};
    for j=1:size(temp_FC,1)-1
        source = [source; enrich_81.source{i}];
    end
    temp_FC2 = addvars(temp_FC, string(source), 'Before', 'Log2FC09');
    temp_FC2 = renamevars(temp_FC2, 'Var1', 'source');

    term = enrich_81.term_name{i};
    for j=1:size(temp_FC,1)-1
        term = [term; enrich_81.term_name{i}];
    end
    temp_FC3 = addvars(temp_FC2, string(term), 'Before', 'Log2FC09');
    temp_FC3 = renamevars(temp_FC3, 'Var2', 'term');

    radio_81 = [radio_81; temp_FC3];
end
radio_81(1,:) = [];
repeats81a = groupsummary(radio_81, ["Log2FC09","Log2FC18","Log2FC27","Log2FC36","Log2FC45","Log2FC54","Log2FC63","Log2FC72","Log2FC81","Log2FC90","Padj09","Padj18","Padj27","Padj36","Padj45","Padj54","Padj63","Padj72","Padj81","Padj90"]);

sources_90 = enrichment_90(ismember(enrichment_90.source, stringsToFind),:);
enrich_90 = sources_90(~ismember(sources_90.term_name, excludeStrings),:);
radio_90 = table('Size', [1 22], 'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'source', 'term', 'Log2FC09', 'Log2FC18', 'Log2FC27', 'Log2FC36', 'Log2FC45', 'Log2FC54', 'Log2FC63', 'Log2FC72', 'Log2FC81', 'Log2FC90', 'Padj09', 'Padj18', 'Padj27', 'Padj36', 'Padj45', 'Padj54', 'Padj63', 'Padj72', 'Padj81', 'Padj90'});
for i=1:size(enrich_90,1)
    temp_FC = DEGenesClean4j(split(enrich_90.intersections{i},','),:);
    temp_FC.Properties.RowNames = {};

    source = enrich_90.source{i};
    for j=1:size(temp_FC,1)-1
        source = [source; enrich_90.source{i}];
    end
    temp_FC2 = addvars(temp_FC, string(source), 'Before', 'Log2FC09');
    temp_FC2 = renamevars(temp_FC2, 'Var1', 'source');

    term = enrich_90.term_name{i};
    for j=1:size(temp_FC,1)-1
        term = [term; enrich_90.term_name{i}];
    end
    temp_FC3 = addvars(temp_FC2, string(term), 'Before', 'Log2FC09');
    temp_FC3 = renamevars(temp_FC3, 'Var2', 'term');

    radio_90 = [radio_90; temp_FC3];
end
radio_90(1,:) = [];
repeats90a = groupsummary(radio_90, ["Log2FC09","Log2FC18","Log2FC27","Log2FC36","Log2FC45","Log2FC54","Log2FC63","Log2FC72","Log2FC81","Log2FC90","Padj09","Padj18","Padj27","Padj36","Padj45","Padj54","Padj63","Padj72","Padj81","Padj90"]);




range_max = 0.035;
num_genes = size(repeats09a,1) + size(repeats18a,1) + size(repeats27a,1) + size(repeats36a,1) + size(repeats45a,1) + size(repeats54a,1) + size(repeats63a,1) + size(repeats72a,1) + size(repeats81a,1) + size(repeats90a,1);
range = linspace(-1,1,num_genes);
plot_range = range_max.*range;
randomIndices = randperm(length(plot_range));
shuffled_plot_range = plot_range(randomIndices)';

repeats09 = addvars(repeats09a, shuffled_plot_range(1:size(repeats09a,1)));
repeats09 = renamevars(repeats09, 'Var22', 'x');
repeats18 = addvars(repeats18a, shuffled_plot_range((size(repeats09,1)+1):(size(repeats09,1) + size(repeats18a,1)))  );
repeats18 = renamevars(repeats18, 'Var22', 'x');
repeats27 = addvars(repeats27a, shuffled_plot_range( (size(repeats09,1)+size(repeats18,1)+1):(size(repeats09,1)+size(repeats18,1) + size(repeats27a,1))  ));
repeats27 = renamevars(repeats27, 'Var22', 'x');
repeats36 = addvars(repeats36a, shuffled_plot_range( (size(repeats09,1)+size(repeats18,1)+size(repeats27,1)+1):(size(repeats09,1)+size(repeats18,1)+size(repeats27,1) + size(repeats36a,1))  ));
repeats36 = renamevars(repeats36, 'Var22', 'x');
repeats45 = addvars(repeats45a, shuffled_plot_range( (size(repeats09,1)+size(repeats18,1)+size(repeats27,1)+size(repeats36,1)+1):(size(repeats09,1)+size(repeats18,1)+size(repeats27,1)+size(repeats36,1) + size(repeats45a,1))  ));
repeats45 = renamevars(repeats45, 'Var22', 'x');
repeats54 = addvars(repeats54a, shuffled_plot_range( (size(repeats09,1)+size(repeats18,1)+size(repeats27,1)+size(repeats36,1)+size(repeats45,1)+1):(size(repeats09,1)+size(repeats18,1)+size(repeats27,1)+size(repeats36,1)+size(repeats45,1) + size(repeats54a,1))  ));
repeats54 = renamevars(repeats54, 'Var22', 'x');
repeats63 = addvars(repeats63a, shuffled_plot_range( (size(repeats09,1)+size(repeats18,1)+size(repeats27,1)+size(repeats36,1)+size(repeats45,1)+size(repeats54,1)+1):(size(repeats09,1)+size(repeats18,1)+size(repeats27,1)+size(repeats36,1)+size(repeats45,1)+size(repeats54,1) + size(repeats63a,1))  ));
repeats63 = renamevars(repeats63, 'Var22', 'x');
repeats72 = addvars(repeats72a, shuffled_plot_range( (size(repeats09,1)+size(repeats18,1)+size(repeats27,1)+size(repeats36,1)+size(repeats45,1)+size(repeats54,1)+size(repeats63,1)+1):(size(repeats09,1)+size(repeats18,1)+size(repeats27,1)+size(repeats36,1)+size(repeats45,1)+size(repeats54,1)+size(repeats63,1) + size(repeats72a,1))  ));
repeats72 = renamevars(repeats72, 'Var22', 'x');
repeats81 = addvars(repeats81a, shuffled_plot_range( (size(repeats09,1)+size(repeats18,1)+size(repeats27,1)+size(repeats36,1)+size(repeats45,1)+size(repeats54,1)+size(repeats63,1)+size(repeats72,1)+1):(size(repeats09,1)+size(repeats18,1)+size(repeats27,1)+size(repeats36,1)+size(repeats45,1)+size(repeats54,1)+size(repeats63,1)+size(repeats72,1) + size(repeats81a,1))  ));
repeats81 = renamevars(repeats81, 'Var22', 'x');
repeats90 = addvars(repeats90a, shuffled_plot_range( (size(repeats09,1)+size(repeats18,1)+size(repeats27,1)+size(repeats36,1)+size(repeats45,1)+size(repeats54,1)+size(repeats63,1)+size(repeats72,1)+size(repeats81,1)+1):(size(repeats09,1)+size(repeats18,1)+size(repeats27,1)+size(repeats36,1)+size(repeats45,1)+size(repeats54,1)+size(repeats63,1)+size(repeats72,1)+size(repeats81,1) + size(repeats90a,1))  ));
repeats90 = renamevars(repeats90, 'Var22', 'x');




marker = 30;
marker2 = 1.2;
scale = 1./2.5;
linew = 1;
sig_cutoffx = 0.05;
fold_cutoffx = 1;
low_alpha = 0.2;
fig_radio = figure;
fig_radio.Theme = "light";
%yyaxis left;
hold on
% Overlay transparent gray rectangle for not significant region
s2 = patch([0 0.945 0.945 0], [1 1 -1 -1], [0.95 0.95 0.95], 'EdgeColor', 'none');
s2.FaceAlpha = 0.8;

% Plot not significant points in color
for i=1:size(repeats36{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats36{i,(10+j)} <= sig_cutoffx)&(abs(repeats36{i,j}) >= fold_cutoffx)
%            scatter(frequencies(j), repeats36{i,j}, marker, [0.6667 1 0.3333], '.', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
        else
            scatter(frequencies(j) + repeats36{i,22}, repeats36{i,j}, marker + repeats36{i,21}.*scale, [0.6667 1 0.3333], 'o', 'filled', 'MarkerFaceAlpha', low_alpha, 'MarkerEdgeAlpha', low_alpha)
        end    
    end
end
for i=1:size(repeats09{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats09{i,(10+j)} <= sig_cutoffx)&(abs(repeats09{i,j}) >= fold_cutoffx)
%            scatter(frequencies(j), repeats09{i,j}, marker, [1 0.3333 0], '.', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
        else
            scatter(frequencies(j) + repeats09{i,22}, repeats09{i,j}, marker + repeats09{i,21}.*scale, [1 0.3333 0], 'o', 'filled', 'MarkerFaceAlpha', low_alpha, 'MarkerEdgeAlpha', low_alpha)
        end    
    end
end
for i=1:size(repeats18{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats18{i,(10+j)} <= sig_cutoffx)&(abs(repeats18{i,j}) >= fold_cutoffx)
%            scatter(frequencies(j), repeats18{i,j}, marker, [1 0.6667 0], '.', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
        else
            scatter(frequencies(j) + repeats18{i,22}, repeats18{i,j}, marker + repeats18{i,21}.*scale, [1 0.6667 0], 'o', 'filled', 'MarkerFaceAlpha', low_alpha, 'MarkerEdgeAlpha', low_alpha)
        end    
    end
end
for i=1:size(repeats27{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats27{i,(10+j)} <= sig_cutoffx)&(abs(repeats27{i,j}) >= fold_cutoffx)
%            scatter(frequencies(j), repeats27{i,j}, marker, [1 1 0], '.', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
        else
            scatter(frequencies(j) + repeats27{i,22}, repeats27{i,j}, marker + repeats27{i,21}.*scale, [1 1 0], 'o', 'filled', 'MarkerFaceAlpha', low_alpha, 'MarkerEdgeAlpha', low_alpha)
        end    
    end
end
for i=1:size(repeats45{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats45{i,(10+j)} <= sig_cutoffx)&(abs(repeats45{i,j}) >= fold_cutoffx)
%            scatter(frequencies(j), repeats45{i,j}, marker, [0.3333 1 0.6667], '.', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
        else
            scatter(frequencies(j) + repeats45{i,22}, repeats45{i,j}, marker + repeats45{i,21}.*scale, [0.3333 1 0.6667], 'o', 'filled', 'MarkerFaceAlpha', low_alpha, 'MarkerEdgeAlpha', low_alpha)
        end    
    end
end
for i=1:size(repeats54{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats54{i,(10+j)} <= sig_cutoffx)&(abs(repeats54{i,j}) >= fold_cutoffx)
%            scatter(frequencies(j), repeats54{i,j}, marker, [0 1 1], '.', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
        else
            scatter(frequencies(j) + repeats54{i,22}, repeats54{i,j}, marker + repeats54{i,21}.*scale, [0 1 1], 'o', 'filled', 'MarkerFaceAlpha', low_alpha, 'MarkerEdgeAlpha', low_alpha)
        end    
    end
end
for i=1:size(repeats63{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats63{i,(10+j)} <= sig_cutoffx)&(abs(repeats63{i,j}) >= fold_cutoffx)
%            scatter(frequencies(j), repeats63{i,j}, marker, [0 0.6667 1], '.', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
        else
            scatter(frequencies(j) + repeats63{i,22}, repeats63{i,j}, marker + repeats63{i,21}.*scale, [0 0.6667 1], 'o', 'filled', 'MarkerFaceAlpha', low_alpha, 'MarkerEdgeAlpha', low_alpha)
        end    
    end
end
for i=1:size(repeats72{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats72{i,(10+j)} <= sig_cutoffx)&(abs(repeats72{i,j}) >= fold_cutoffx)
%            scatter(frequencies(j), repeats72{i,j}, marker, [0 0.3333 1], '.', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
        else
            scatter(frequencies(j) + repeats72{i,22}, repeats72{i,j}, marker + repeats72{i,21}.*scale, [0 0.3333 1], 'o', 'filled', 'MarkerFaceAlpha', low_alpha, 'MarkerEdgeAlpha', low_alpha)
        end    
    end
end
for i=1:size(repeats81{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats81{i,(10+j)} <= sig_cutoffx)&(abs(repeats81{i,j}) >= fold_cutoffx)
 %           scatter(frequencies(j), repeats81{i,j}, marker, [0 0 1], '.', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
        else
            scatter(frequencies(j) + repeats81{i,22}, repeats81{i,j}, marker + repeats81{i,21}.*scale, [0 0 1], 'o', 'filled', 'MarkerFaceAlpha', low_alpha, 'MarkerEdgeAlpha', low_alpha)
        end    
    end
end
for i=1:size(repeats90{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats90{i,(10+j)} <= sig_cutoffx)&(abs(repeats90{i,j}) >= fold_cutoffx)
%            scatter(frequencies(j), repeats90{i,j}, marker, [0 0 0.6667], '.', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
        else
            scatter(frequencies(j) + repeats90{i,22}, repeats90{i,j}, marker + repeats90{i,21}.*scale, [0 0 0.6667], 'o', 'filled', 'MarkerFaceAlpha', low_alpha, 'MarkerEdgeAlpha', low_alpha)
        end    
    end
end

% Then plot significant points in color
for i=1:size(repeats36{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats36{i,(10+j)} <= sig_cutoffx)&(abs(repeats36{i,j}) >= fold_cutoffx)
            scatter(frequencies(j) + repeats36{i,22}, repeats36{i,j}, marker + repeats36{i,21}.*scale, [0.6667 1 0.3333], 'o', 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'LineWidth', linew, 'MarkerEdgeColor','k')
        end    
    end
end
for i=1:size(repeats09{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats09{i,(10+j)} <= sig_cutoffx)&(abs(repeats09{i,j}) >= fold_cutoffx)
            scatter(frequencies(j) + repeats09{i,22}, repeats09{i,j}, marker + repeats09{i,21}.*scale, [1 0.3333 0], 'o', 'filled','MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'LineWidth', linew, 'MarkerEdgeColor','k')
        end    
    end
end
for i=1:size(repeats18{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats18{i,(10+j)} <= sig_cutoffx)&(abs(repeats18{i,j}) >= fold_cutoffx)
            scatter(frequencies(j) + repeats18{i,22}, repeats18{i,j}, marker + repeats18{i,21}.*scale, [1 0.6667 0], 'o', 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'LineWidth', linew, 'MarkerEdgeColor','k')
        end    
    end
end
for i=1:size(repeats27{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats27{i,(10+j)} <= sig_cutoffx)&(abs(repeats27{i,j}) >= fold_cutoffx)
            scatter(frequencies(j) + repeats27{i,22}, repeats27{i,j}, marker + repeats27{i,21}.*scale, [1 1 0], 'o', 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'LineWidth', linew, 'MarkerEdgeColor','k')
        end    
    end
end
for i=1:size(repeats45{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats45{i,(10+j)} <= sig_cutoffx)&(abs(repeats45{i,j}) >= fold_cutoffx)
            scatter(frequencies(j) + repeats45{i,22}, repeats45{i,j}, marker + repeats45{i,21}.*scale, [0.3333 1 0.6667], 'o', 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'LineWidth', linew, 'MarkerEdgeColor','k')
        end    
    end
end
for i=1:size(repeats54{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats54{i,(10+j)} <= sig_cutoffx)&(abs(repeats54{i,j}) >= fold_cutoffx)
            scatter(frequencies(j) + repeats54{i,22}, repeats54{i,j}, marker + repeats54{i,21}.*scale, [0 1 1], 'o', 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'LineWidth', linew, 'MarkerEdgeColor', 'k')
        end    
    end
end
for i=1:size(repeats63{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats63{i,(10+j)} <= sig_cutoffx)&(abs(repeats63{i,j}) >= fold_cutoffx)
            scatter(frequencies(j) + repeats63{i,22}, repeats63{i,j}, marker + repeats63{i,21}.*scale, [0 0.6667 1], 'o', 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'LineWidth', linew, 'MarkerEdgeColor','k')
        end    
    end
end
for i=1:size(repeats72{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats72{i,(10+j)} <= sig_cutoffx)&(abs(repeats72{i,j}) >= fold_cutoffx)
            scatter(frequencies(j) + repeats72{i,22}, repeats72{i,j}, marker + repeats72{i,21}.*scale, [0 0.3333 1], 'o', 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'LineWidth', linew, 'MarkerEdgeColor','k')
        end    
    end
end
for i=1:size(repeats81{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats81{i,(10+j)} <= sig_cutoffx)&(abs(repeats81{i,j}) >= fold_cutoffx)
            scatter(frequencies(j) + repeats81{i,22}, repeats81{i,j}, marker + repeats81{i,21}.*scale, [0 0 1], 'o', 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'LineWidth', linew, 'MarkerEdgeColor','k')
        end    
    end
end
for i=1:size(repeats90{:,:},1)
    for j=1:size(frequencies,2)
        if (repeats90{i,(10+j)} <= sig_cutoffx)&(abs(repeats90{i,j}) >= fold_cutoffx)
            scatter(frequencies(j) + repeats90{i,22}, repeats90{i,j}, marker + repeats90{i,21}.*scale, [0 0 0.6667], 'o', 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'LineWidth', linew, 'MarkerEdgeColor','k')
        end    
    end
end

box on
set(gca, 'LineWidth', 2);
axis([0.045 0.945 -5 12])
xticks(frequencies);
yticks(linspace(-5, 12, 18))
xlabel('Frequency (hr^{-1})')
ylabel('Log_2(Fold Change)')
fontsize(16, "points")
ax = gca;
ax.YGrid = 'on';
ax.GridLineWidth = 1;
ax.TickLength = [0.0025 0.0025];
pbaspect([5.7 1 1])





num_genes = [size(Only09Table,1) 47 size(repeats09,1);...
                size(Only18Table,1) 42 size(repeats18,1);...
                size(Only27Table,1) 19 size(repeats27,1);...
                size(Only36Table,1) 331 size(repeats36,1);...
                size(Only45Table,1) 65 size(repeats45,1);...
                size(Only54Table,1) 277 size(repeats54,1);...
                size(Only63Table,1) 6 size(repeats63,1);...
                size(Only72Table,1) 28 size(repeats72,1);...
                size(Only81Table,1) 16 size(repeats81,1);...
                size(Only90Table,1) 20 size(repeats90,1)];



freq_offset = 0;
bar_width = 0.08;

fig_bar = figure;
fig_bar.Theme = "light";
hold on
b1 = bar(num_genes(:,1));
b2 = bar(num_genes(:,2));
b3 = bar(num_genes(:,3));

b1.FaceColor = 'flat';
b1.CData = flipud(jet(10) + (1 - jet(10))*0.8);
b1.LineWidth = 0.5;
b2.FaceColor = 'flat';
b2.CData = flipud(jet(10) + (1 - jet(10))*0.5);
b2.LineWidth = 0.5;
b3.FaceColor = 'flat';
b3.CData = flipud(jet(10));
b3.LineWidth = 0.5;
axis([0.2 10.8 1 1300])
box on
set(gca, 'LineWidth', 0.5);

set(gca, 'YScale', 'log')





ylabel('Number of Intersecting DE Genes')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

