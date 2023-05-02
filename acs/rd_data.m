%% Import data from text file
% Script for importing data from the following directory:
%
%    file path: /Users/nicholasbaetge/Box Sync/Phyto_bbp/DATA/FINAL/2022_Experiments/acs
%
% Re-run this file for each of the bottles of each experiment

clear all
close all

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 12);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["exp", "bottle", "tp", "datetime", "date", "time", "wl", "c", "c_se", "a", "a_se", "b"];
opts.VariableTypes = ["categorical", "categorical", "double", "categorical", "categorical", "datetime", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["bottle", "datetime"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "time", "InputFormat", "HH:mm:ss");

% Import the data
A = readtable("/Users/nicholasbaetge/Box Sync/Phyto_bbp/DATA/FINAL/2022_Experiments/acs/SYN22-5_A.csv", opts);

%% Run correction
clear opts
tp=A{:,3};
AA=unique(tp); %this is the sample number
for i=1:length(AA)
    I=find(tp==AA(i));
    a(i,:)=A{I,10};
    err_a(i,:)=A{I,11};
    c(i,:)=A{I,8};
    err_c(i,:)=A{I,9};
end
wl=A{I,7};
subplot(1,2,1) 
plot(wl,a)
subplot(1,2,2)
plot(wl,c)

opts = optimset('fminsearch');      
opts = optimset(opts,'MaxIter',20000000); 
opts = optimset(opts,'MaxFunEvals',20000);
opts = optimset(opts,'TolX',1e-8);
opts = optimset(opts,'TolFun',1e-8);

wla=wl;
tmp = xlsread('Sullivan_etal_2006_instrumentspecific.xls');
phi_T=interp1(tmp(:,1),tmp(:,2),wl,'linear');

for i=1:15 %change based on sample number
    spectra=a(i,:);
    I = find(wla>=710 & wla<=750);  % spectral range for optimization (710 to 750nm)   
    delT = 0;
    offset = 0;
    [x1] = fminsearch(@f_T, [delT, offset], opts, spectra(I), phi_T(I)');
    spectra=spectra-x1(1).*phi_T';
    a(i,:)=spectra;
    spectra=c(i,:);
    delT = 0;
    offset = 0;
    [x1] = fminsearch(@f_T, [delT, offset], opts, spectra(I), phi_T(I)');
    spectra=spectra-x1(1).*phi_T';
    c(i,:)=spectra;
end

for i=1:15 %change based on sample number
    [ap_corr(i,:), cp_corr(i,:)] = ResidualTemperatureAndScatteringCorrection(a(i,1:80), c(i,1:80), wl(1:80)');
end

figure
subplot(1,2,1) 
plot(wl(1:80),ap_corr)
subplot(1,2,2)
plot(wl(1:80),cp_corr)


%% Save data
export_ap = array2table(ap_corr);
writetable(export_ap, "/Users/nicholasbaetge/Box Sync/Phyto_bbp/DATA/FINAL/2022_Experiments/acs/corrected/SYN22-5_A_ap.csv"); %evaluate these in console

export_cp = array2table(cp_corr);
writetable(export_cp, "/Users/nicholasbaetge/Box Sync/Phyto_bbp/DATA/FINAL/2022_Experiments/acs/corrected/SYN22-5_A_cp.csv");

function costf = f_T(x0, spectra, psiT);
    costf = sum((spectra - psiT.*x0(1) -x0(2)).^2);
end

