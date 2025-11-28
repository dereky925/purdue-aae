clear;clc

% stem_leaf_diagram.m
% Builds a stem‐and‐leaf plot using the thousands‐digit as stem and
% the hundreds‐digit as leaf.

% 1) Enter your data as a row vector
data = [...
    1280 5320 4390 2100 1240 3060 4970 ...
    1050  360 3330 3380  340 1000  960 ...
    1320  530 3350  540 3870 1250 2400 ...
     960 1120 2120  450 2250 2320 2400 ...
    3150 5700 5220  500 1850 2460 5750 ...
    2800 2730 1670  100 5770 3150 1990 ...
     510  240  396 1419 2109 ...
];

% 2) Compute stem and leaf
stems  = floor(data/1000);                % thousands digit
leaves = floor(mod(data,1000)/100);       % hundreds digit

% 3) Find unique stems in ascending order
unique_stems = unique(stems);

% 4) Print header
fprintf('Stem | Leaves\n');
fprintf('-----+-----------------------------\n');

% 5) For each stem, collect & sort its leaves, then print
for s = unique_stems
    idx    = (stems == s);
    leaf_s = sort(leaves(idx));
    
    % print stem
    fprintf('%4d | ', s);
    % print all leaves separated by spaces
    fprintf('%d ', leaf_s);
    fprintf('\n');
end


% 2) Define the bin edges
edges = 0:1000:6000;   % creates [0 1000 2000 … 6000]

% 3a) Quick & easy: use built-in histogram()
figure
h = histogram(data, edges);
h.FaceColor = [0.9 0.7 0.4];
h.EdgeColor = 'k';
xlabel('Length');
ylabel('Frequency');
title('Subdivision Lengths Histogram');
grid on


%%
clear;clc

% stem_leaf_specific_gravity.m
% Builds a repeated‐stem stem‐and‐leaf display for specific-gravity data

% 1) Input your data (specific gravities)
data = [ ...
    0.32 0.35 0.36 0.36 0.37 0.38 0.40 0.40 0.40 ...
    0.41 0.41 0.42 0.42 0.42 0.42 0.42 0.43 0.44 ...
    0.45 0.46 0.46 0.47 0.48 0.48 0.49 0.53 0.54 ...
    0.54 0.55 0.58 0.61 0.66 0.66 0.67 0.68 0.78 ...
];

% 2) Scale to integer hundredths
scaled = round(data * 100);    % e.g. 0.32 -> 32

% 3) Extract stem and leaf
stems  = floor(scaled/10);     % tenths digit (3,4,5,6,7,…)
leaves = mod(scaled,10);       % hundredths digit (0–9)

% 4) Prepare to display stems 3 through 7
stem_values = 3:7;

% 5) Print header
fprintf(' Stem | Leaves (leaf = hundredths)\n');
fprintf('------+---------------------------\n');

% 6) Loop over each stem, splitting into Low (L) and High (H)
for s = stem_values
    % Low leaves (0..4)
    low_leaves  = sort(leaves(stems==s & leaves<=4));
    % High leaves (5..9)
    high_leaves = sort(leaves(stems==s & leaves>=5));
    
    % Format leaf lists or 'NONE' if empty
    if isempty(low_leaves)
        low_str = 'NONE';
    else
        low_str = sprintf('%d ', low_leaves);
    end
    if isempty(high_leaves)
        high_str = 'NONE';
    else
        high_str = sprintf('%d ', high_leaves);
    end
    
    % Print the two rows for this stem
    fprintf('  %dL  | %s\n', s, strtrim(low_str));   % e.g. " 3L  | 2"
    fprintf('  %dH  | %s\n', s, strtrim(high_str));  % e.g. " 3H  | 5 6 6 7 8"
end




%%

clear;clc;close all



% endotoxin_analysis.m
% Computes means, medians, compares them, and explains the mean/median gap.

% 1) Input the data
urban = [6.0  5.0 11.0 33.0  4.0  5.0 80.0 18.0 35.0 17.0 23.0];
farm  = [3.0 13.0 12.0  8.0  8.0  6.0  5.0 18.0  4.0  7.6 24.0 8.6 2.0 1.0 0.2];

% 2) Compute means and medians
mu_urban = mean(urban);
mu_farm  = mean(farm);
med_urban = median(urban);
med_farm  = median(farm);

% 3) Print results (rounded to 4 decimal places)
fprintf('Urban mean:  %.4f EU/mg\n', mu_urban);
fprintf('Farm  mean:  %.4f EU/mg\n\n',  mu_farm);

fprintf('Urban median: %.4f EU/mg\n', med_urban);
fprintf('Farm  median: %.4f EU/mg\n\n',  med_farm);

% 4) Compare averages
if      mu_urban > 2*mu_farm
    fprintf('→ The average endotoxin concentration in URBAN homes is more than double that in FARM homes.\n');
elseif  mu_farm  > 2*mu_urban
    fprintf('→ The average endotoxin concentration in FARM homes is more than double that in URBAN homes.\n');
else
    fprintf('→ The average endotoxin concentration is about the same in both URBAN and FARM homes.\n');
end

% 5) Compare medians
if      med_urban > 2*med_farm
    fprintf('→ The median endotoxin concentration in URBAN homes is roughly double that in FARM homes.\n');
elseif  med_farm  > 2*med_urban
    fprintf('→ The median endotoxin concentration in FARM homes is roughly double that in URBAN homes.\n');
else
    fprintf('→ The median endotoxin concentration is about the same in both URBAN and FARM homes.\n');
end

% 6) Explain the mean–median discrepancy
fprintf('\nReason the URBAN mean and median differ so much:\n');
fprintf('  The few very large values in the urban sample raise the mean but have little effect on the median.\n');

%% 
clear;clc

% trimmed_mean_endotoxin.m
% (c) Delete the smallest & largest obs, compute trimmed mean, 
%     trimming percentage, and compare to mean & median.

% 1) Input data
urban = [6.0  5.0 11.0 33.0  4.0  5.0 80.0 18.0 35.0 17.0 23.0];
farm  = [3.0 13.0 12.0  8.0  8.0  6.0  5.0 18.0  4.0  7.6 24.0 8.6 2.0 1.0 0.2];

% 2) Sort each sample
u = sort(urban);
f = sort(farm);

% 3) Remove smallest & largest
u_trim = u(2:end-1);   % removes 1st and last
f_trim = f(2:end-1);

% 4) Compute trimmed means
mu_u_trim = mean(u_trim);
mu_f_trim = mean(f_trim);

% 5) Compute trimming percentages
pct_u = (1/numel(urban))*100;   % two observations out of total
pct_f = (1/numel(farm))*100;

% 6) Original means & medians for comparison
mu_u = mean(urban);   med_u = median(urban);
mu_f = mean(farm);    med_f = median(farm);

% 7) Display results
fprintf('Trimmed means (delete min & max):\n');
fprintf('  Urban: %6.4f EU/mg   (trim = %5.2f%%)\n', mu_u_trim, pct_u);
fprintf('   Farm: %6.4f EU/mg   (trim = %5.2f%%)\n\n', mu_f_trim, pct_f);

fprintf('Original means & medians:\n');
fprintf('  Urban mean = %6.4f, median = %6.4f\n', mu_u, med_u);
fprintf('   Farm mean = %6.4f, median = %6.4f\n\n', mu_f, med_f);

% 8) Compare trimmed mean to original mean & median
compare = @ (trim, orig, name) ...
    fprintf('%s trimmed mean is %s the original %s.\n', ...
            name, ...
            ternary(trim>orig, 'greater than', ...
              ternary(trim<orig,'less than','equal to')), ...
            name);

fprintf('Comparisons (trimmed mean vs original mean/median):\n');
compare(mu_u_trim, mu_u, 'Urban mean');
compare(mu_u_trim, med_u, 'Urban median');
compare(mu_f_trim, mu_f, 'Farm mean');
compare(mu_f_trim, med_f, 'Farm median');


% Helper inline ternary function
    function out = ternary(cond, a, b)
        if cond, out = a; else out = b; end
    end






