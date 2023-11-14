%% Import data from text file
% Script for importing data from the following folder:
%
%    filename: /Users/nicholasbaetge/Box
%    Sync/Phyto_bbp/DATA/FINAL/diel_optics_phyto_cultures/analyses/derivations/
%

clear all
close all

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["tp", "wl", "mean_ap", "sd_ap"];
opts.VariableTypes = ["double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";


% Import the data
A = readtable("//Users/nicholasbaetge/Box Sync/Phyto_bbp/DATA/FINAL/diel_optics_phyto_cultures/analyses/derivations/OL22-3.csv", opts);
B = A(:,1:3);
val_data = unstack(B, "mean_ap", "wl",  'VariableNamingRule', 'preserve');
vals = val_data{:, 2:84};

D = A(:, [1 2 4]);
sd_data = unstack(D, "sd_ap", "wl",  'VariableNamingRule', 'preserve');
sd = sd_data{:, 2:84};

lambda = unique(A.wl);
spec = vals;

%% unsmooth
spec = unsmoothACS(spec, lambda);

%%
wl = lambda;
ap = spec;
wlunc =lambda;
apunc = sd;

[amps, compspec, sumspec] = GaussDecomp(wl,ap,wlunc,apunc)
%[amps, compspec, sumspec] = DiatomGaussDecomp(wl,ap,wlunc,apunc)

%%
figure
plot(lambda,spec(10,:),'k')
hold on
plot(lambda,sumspec(10,:),'r')
plot(lambda,compspec(:, :, 10))

%%
figure
plot(lambda,spec(4,:),'k')
hold on
plot(lambda,sumspec(4,:),'r')
plot(lambda,compspec(:, :, 4))

%%
export_unsmooth = array2table(ap);
writetable(export_unsmooth, "decomposed/OL22-3_unsmoothed.csv");

export_amps = array2table(amps);
writetable(export_amps, "decomposed/OL22-3_amps.csv");

export_sumspec = array2table(sumspec);
writetable(export_sumspec, "decomposed/OL22-3_sumspec.csv");

save('decomposed/OL22-3_decomp_comspec.mat', 'compspec' );

%%

for i = 1:size(compspec,3)
    filename = ['decomposed/OL22-3_compspec',num2str(i),'.csv'] ;
    csvwrite(filename,compspec(:,:, i))
end

