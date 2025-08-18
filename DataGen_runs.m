% (sigF, sigM, optM, tMax, n, nLoci, nLociP, nBins, mutRt, rngSeed, condDens)
%
%sigFs = [0.1, 0.5, 1, 5];
sigFs = [5];
sigM = 10;
optM = 0.5;
tMax = 50000;
n = 1000;
nLociT = 100;
nLociP = 100;
nBins = [2 3 4 5 10 100];
%mutRts = [0 0.00001 0.0001 0.001];
mutRts = [0];
condDens = 4;
rngSeed = ceil(sum([sigFs,sigM,optM,tMax,n,nLociT,nLociP,nBins,mutRts,condDens])); % use [] if none

% sigFs = [0.1, 0.5];
% nBinss = [2 3];
% mutRts = [0 0.00001];

results2 = [];
tic
for hh = 1:size(sigFs,2)
    for ii = 1:size(nBinss,2)
        for jj = 1:size(mutRts,2)
            results1 = DataGen(sigFs(hh), sigM, optM, tMax, n, nLociT, nLociP, nBins(ii), mutRts(jj), condDens, rngSeed);
            %results1 = DataGen(sigFs(hh), 10, 0.5, 100, 10, 10, 10, nBins(ii), mutRts(jj), 8122025, 3);
            results2 = [results2; results1];
        end
    end
end
toc

%start 2:15 8/14
% time taken apparently 3524 sec? 

% start with t -> 5000 at 8/15 11:17 am

filename=[];
filename = sprintf('results2_50k_%dx%d_%s.csv', condDens, condDens, datestr(now,'mm-dd-yyyy HH-MM'))
'UNCOMMENT AND EXPORT!'
writematrix(results2,filename)