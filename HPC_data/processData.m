
%Process Simulation Data collected from NYU HPC and save the output to 
% figures2/*csv files as well as plot them for visualizing


close all
clear
wannaplot=1;

% Count number of csv files in current directory
files = dir('*.csv');
nFiles = size(files,1);

densityBL = [0.005,0.01];
densityAP = [50,100,200,300]*10^(-6);%(1:1:10)/10^4;
omegaVal = [0, pi/3];

mu=2;

nBL = length(densityBL);
nAP = length(densityAP);
nO = length(omegaVal);
nAll = nBL*nAP*nO;

data = cell(nFiles,1);
datanew = cell(1,nAP,nBL,nO);

% Read in all files
% Store each file in an entry in 'data'
% (each entry in 'data' has nBL*nAP*nO columns)
for i=1:nFiles
        if (exist(strcat('output',int2str(i),'.csv'))==0)
            warning(strcat('output',int2str(i),'.csv - File does not exist'))
            continue;
        
        else
            data{i}=csvread(strcat('output',int2str(i),'.csv'));
                
            colNum0 = find(data{i}(5,:)==0); 
            NaNarray = isnan(data{i}(2,:));
            
            data{i}(2,NaNarray) = 1./(mu*data{i}(5,NaNarray));
            
            data{i}(3,colNum0) = 1; %when n=0, prob of blockage=1;
            
         
            for ii = 1:(nBL*nAP*nO)
                if(data{i}(5,ii)~=0)
                    datanew{ii} = [datanew{ii},data{i}(1:3,ii)];
                end               
            end
        end
end

% Compute the mean value, standard deviation, and n of
%   avgFreq: Average frequency of simultaneous blockage of all BS. 
%   avgDur: Average duration of simultaneous blockage of all BS.
%   probAllBl: probability of simultaneous blockage of all BS.
% across all the iterations.
% We end up with a 3x(nBL*nAP*nO) matrix for mean and std, 
% and a 1x(nBL*nAP*nO) matrix for n
for ii = 1:nAll
    meanval(:,ii) =  mean(datanew{ii},2);
    stdval(:,ii) =  std(datanew{ii},0,2);
    size_data(ii) = size(datanew{ii},2);
end

% Create a 16-entry list for each value
freqBl = meanval(1, :);
durBl = meanval(2, :);
probBl = meanval(3, :);
% BL density, AP density, and omega for each entry
N = 1:nAll;
omegaIdx = ceil(N/(nAP*nBL));
densBLIdx = ceil(N/nAP) - nO*(ceil(N/(nAP*nBL))-1);
APIs = mod(N , nAP) + 1;
densAPIdx = [ APIs(nAll) APIs(1:nAll-1)];



% densAPIdx = ceil(N / nAP);
% densBLIdx = ceil(N / nBL) - (nBL)*(ceil( N / nAP ) - 1);
% omegas = mod(N, nO) + 1;
% omegaIdx= [omegas(nAll) omegas(1:nAll-1)];

figure
scatter(densityAP(densAPIdx), freqBl, 100, 2*omegaIdx+densBLIdx, 'filled')
title('Mean Blockage Frequency')
xlabel('BS density')