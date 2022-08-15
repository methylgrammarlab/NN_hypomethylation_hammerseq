fn = 'all.sorted.bedgraph'; RUNDESC = 'hammerseq by NN score'; timeorder = 1:5;
FLDS = {'chr','start','end','m0','m5','m30','h2','h24','nn', 'solowcgw'};


DATA_COLS = [4:8];
WEIGHT_BY_COUNT = 1; MIN_COUNT = 0; 

EXTREMES = [0 0.00076 0.00176 0.01183 0.02746 0.95698 0.98212 0.99667 0.99817 1];
EXTREMES_PERC = [0 0.005 0.01 0.05 0.1 0.9 0.95 0.99 0.995 1];

if (1)
	tab = readtable(sprintf('%s',fn),'FileType','text','delimiter','\t', ...
		'ReadVariableNames',false, 'ReadRowNames',false,'TreatAsEmpty',{'.','NA','N/A'},...
		'Format','%s%d%d%s%s%s%s%s%f%f'); % somehow it wasn't handling ratio format correctly
	tab.Properties.VariableNames = FLDS;
	
	tab(:,DATA_COLS) = tab(:,DATA_COLS(timeorder));
	
	tab = tab(~isnan(tab.nn),:);
	[num denom] = ratiocolsToArrays(tab, DATA_COLS); % this is slow
end

num_mincount = num;
denom_mincount = denom;

num_mincount(denom<MIN_COUNT) = nan;
denom_mincount(denom<MIN_COUNT) = nan;

ratios = num_mincount./denom_mincount;

xlabs = tab.Properties.VariableNames(DATA_COLS);

means = zeros(length(EXTREMES)-1,length(DATA_COLS));
counts = zeros(length(EXTREMES)-1,length(DATA_COLS));


mean_labs = cell(1,length(EXTREMES)-1);
for i = 2:(length(EXTREMES)+1) % add one to the end for soloWCGW
	
	% Special case for solos
	if (i>length(EXTREMES))
		rows = ~isnan(tab.solowcgw);
	else
		rows = (tab.nn >= EXTREMES(i-1)) & (tab.nn <= EXTREMES(i));
	end

	if (WEIGHT_BY_COUNT == 0)
		% mean of ratios.  Be careful with this if the input data includes
		% denominators of 1
		means(i-1,:) = nanmean(ratios(rows,:));
		counts(i-1,:) = nansum(~isnan(ratios(rows,:)));
	else
		% Average CpGs weighted to number of reads
		counts(i-1,:) = nansum(denom_mincount(rows,:));
		means(i-1,:) = nansum(num_mincount(rows,:)) ./ nansum(denom_mincount(rows,:));
	end
	
	% Special case for solos
	if (i>length(EXTREMES))
		mean_labs{i-1} = sprintf('soloWCGW (%0.0e)',sum(counts(i-1,:)));
	else
		mean_labs{i-1} = sprintf('%0.4f - %0.4f (%0.0e)',EXTREMES(i-1),EXTREMES(i),sum(counts(i-1,:)));
		mean_labs{i-1} = sprintf('<%0.2f ptile (%0.0e)',EXTREMES_PERC(i).*100,sum(counts(i-1,:)));
	end
	
end


if (1)
	ti = sprintf('%s WEIGHTBYCOUNT%d MINCOUNT%d',RUNDESC, WEIGHT_BY_COUNT, MIN_COUNT);
	figure('Position',[500, 500, 250, 250],'Name',ti,'FileName',[ti '.eps']); % 1400,700
	title(ti);
	
	cmap=rev_colormap(redbluecmap(length(mean_labs)-1));
	cmap(length(mean_labs),:) = [0 0 0];
	linestyles = {};
	for i = 1:length(mean_labs)
		linestyles{i} = '-';
	end
	linestyles{length(mean_labs)} = ':'; % soloWCGW
	
	p = plot(means');
	colororder(cmap);
	set(gca,'XTick',1:length(xlabs));
	set(gca,'XTickLabels',xlabs);
	rotateXLabels(gca(),90);
	legend(mean_labs,'Location','southeastoutside');
	ylabel('Maintenance ratio');
	set(gca,'YLim',[0.3 1]);

end

