
%Process Simulation Data collected from NYU HPC and save the output to 
% figures2/*csv files as well as plot them for visualizing

close all
clear
wannaplot=1;
files = ls('*.csv');
nFiles = size(files,1);


densityBL = [0.005,0.01,0.02];
densityAP = [50,100,200,300,400,500]*10^(-6);%(1:1:10)/10^4;
omegaVal = [0, pi/3];

mu=2;

nBL = length(densityBL);
nAP = length(densityAP);
nO = length(omegaVal);

data = cell(nFiles,1);
datanew = cell(1,nBL*nAP*nO);


for i=1:nFiles
        if (exist(strcat('output',int2str(i),'.csv'))==0)
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


for ii = 1:(nBL*nAP*nO)
    
    meanval(:,ii) =  mean(datanew{ii},2);
    stdval(:,ii) =  std(datanew{ii},0,2);
    size_data(ii) = size(datanew{ii},2);
end

confidense_int = 1.96*stdval./meshgrid(sqrt(size_data),1:size(stdval,1));
res = [meanval; confidense_int];
res1 = reshape(res,size(res,1),nAP,nBL,nO);



for jj = 1:size(res,1) %3 for freq, dur, pB and 3 for their intrvl
    res2{jj} = squeeze(res1(jj,:,:,:));
    res3{jj} = reshape(res2{jj},nAP,nBL*nO);
    
end
pBCond = [res3{3},res3{6}]; %[mean;intrvl]
freqCond = [res3{1},res3{4}];
durCond = [res3{2},res3{5}]*1000; % convert to ms;

% pB=squeeze(meanVal(3,:,:,:));
% % std_pB = squeeze(stdVal(3,:,:,:));
% freq=squeeze(meanVal(1,:,:,:));
% % std_freq=squeeze(stdVal(1,:,:,:));
% freqCond=squeeze(meanVal0(1,:,:,:));
% durCond=squeeze(meanVal0(2,:,:,:));
% pBCond=squeeze(meanVal0(3,:,:,:));
% freqCond_int = squeeze(confidense_int_new(1,:,:,:));
% durCond_int = squeeze(confidense_int_new(2,:,:,:));
% pBCond_int = squeeze(confidense_int_new(3,:,:,:));
% %%Make it 2D arrays
% pB=reshape(pB(1:end-1,:,[1,2]),length(densityAP),6);
% pBCond=reshape(pBCond(1:end-1,:,[1,2]),length(densityAP),6);
% freq=reshape(freq(1:end-1,:,[1,2]),length(densityAP),6);
% freqCond=reshape(freqCond(1:end-1,:,[1,2]),length(densityAP),6);
% durCond=reshape(durCond(1:end-1,:,[1,2]),length(densityAP),6);





legendArray= {'lamB0.01omega0','lamB0.005omega0',...
    'lamB0.01omega60','lamB0.005omega60'};
colTitle= {'lamT','lamB0.005omega0','lamB0.01omega0',...
    'lamB0.005omega60','lamB0.01omega60',...
    'lamB0.005omega0int','lamB0.01omega0int',...
    'lamB0.005omega60int','lamB0.01omega60int'};




if(wannaplot)
    %     figure(1);
    %     semilogy(densityAP,pB);
    %     set(gca,'YScale','log');
    % %     ylim([1e-4,1]);title('Marginal prob of Blockage')
    %     legend(legendArray);
    figure(2);
    semilogy(densityAP,pBCond(:,1:4)); title('Conditional prob of Bl given n!=0')
    %     ylim([1e-4,1])
    legend(legendArray);
    
    
    %     figure(3);
    %     semilogy(densityAP,freq)
    %     title('Expected Freq of blockage')
    %     legend(legendArray);
    % %     ylim([1e-4,1])
    figure(4);
    semilogy(densityAP,freqCond(:,1:4));
    title('Conditional expectation of freq of bl given n!=0')
    %     ylim([1e-4,1])
    legend(legendArray);
    
    figure(5);
    semilogy(densityAP,durCond(:,1:4))
    title('Conditional expectation of duration of bl given n!=0')
    legend(legendArray);
    
end





