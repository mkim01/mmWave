%processData
%Processing 200 files generated by HPC cluster by running the Simulation
%200 times

close all
clear
nFiles = 200;
densityBL = [0.01,0.1,0.2,0.5,0.65];
densityAP = [1,2,5,10]/10^4;
for i=1:nFiles
    data(:,:,i)=csvread(strcat('Data/output',int2str(i),'.csv'));
end
meanVal = reshape(mean(data,3),6,5,4);
%The 6 columns are [avgFreq,avgDur,probAllBl,th_freqBl,th_durBl,th_probAllBl];
avgFreq = squeeze(meanVal(1,:,:));
avgDur= squeeze(meanVal(2,:,:));
probAllBl= squeeze(meanVal(3,:,:));
% th_freqBl= squeeze(meanVal(4,:,:));
th_durBl= squeeze(meanVal(5,:,:));
th_probAllBl= squeeze(meanVal(6,:,:));

%%We did wrong calculation of theoretical freq of blockage. Correcting it
%%here and the origin code too. It is not \sum\lambda_n, but
%%\sum\lambda_nP_{N-n} which has an average of

V = 1; %velocity m/s
hb = 1.8;
hr = 1.4;
ht = 6;
frac = (hb-hr)/(ht-hr);
simTime = 60*10; %sec Total Simulation time
tstep = 0.0001; %(sec) time step
mu = 2; %Expected bloc dur =1/mu
R = 100; %m Radius
for indB=1:length(densityBL)
    for indT = 1:length(densityAP)
        C(indB) = 2/pi*densityBL(indB)*V*frac;
        a = 1-2*mu./(R*C(indB)) + 2*mu^2./(R^2*C(indB).^2).*log(1+C(indB).*R/mu);
        th_freqBl(indB,indT) = mu*a.*densityAP(indT)*pi*R^2.*exp((a-1).*densityAP(indT)*pi*R^2);
%         th_freqBl(indB,indT)=exp((a(indB)-1).*densityAP(indT)*pi*R^2);
    end
end



% h=figure(2);
% hold on; grid on;
% plot(densityAP, avgFreq,'LineWidth',2)
% plot(densityAP, th_freqBl,'--','LineWidth',2)
% legend('\rho_b=0.01 bl/m^2','\rho_b=0.1 bl/m^2','\rho_b=0.2 bl/m^2',...
%     '\rho_b=0.5 bl/m^2','\rho_b=0.65 bl/m^2' )
% xlabel('\lambda_T (Density of APs per km^2)', 'fontsize',13)
% ylabel('Prob all-blocked given atleast 1 AP in 100m','fontsize',13)
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % print(h,[pwd '/figures/N_block_prob_cond.pdf'],'-dpdf','-r0')
% % print(h,[pwd '/figures/N_block_prob_cond.png'],'-dpng','-r0')
color = {'r','g','b','m','k'};

h=figure(3);
hold on; grid on;
for j=1:5
plot(densityAP, probAllBl(j,:),'color',color{j},'LineWidth',2)
end
for j=1:5
plot(densityAP, th_probAllBl(j,:),'color',color{j},'LineStyle','--','LineWidth',2)
end
g=legend('\rho_b=0.01 bl/m^2','\rho_b=0.1 bl/m^2','\rho_b=0.2 bl/m^2',...
    '\rho_b=0.5 bl/m^2','\rho_b=0.65 bl/m^2' );
set(g,'fontsize',13)
xlabel('\lambda_T (Density of APs per km^2)', 'fontsize',13)
ylabel('Prob all APs blocked within 100m','fontsize',13)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,[pwd '/figures/N_block_prob_cond.pdf'],'-dpdf','-r0')
% print(h,[pwd '/figures/N_block_prob_cond.png'],'-dpng','-r0')

h=figure(4);
hold on; grid on;
for j=1:5
plot(densityAP, avgFreq(j,:),'color',color{j},'LineWidth',2)
end
for j=1:5
plot(densityAP, th_freqBl(j,:),'color',color{j},'LineStyle','--','LineWidth',2)
end
g=legend('\rho_b=0.01 bl/m^2','\rho_b=0.1 bl/m^2','\rho_b=0.2 bl/m^2',...
    '\rho_b=0.5 bl/m^2','\rho_b=0.65 bl/m^2' );
set(g,'fontsize',13)
xlabel('\lambda_T (Density of APs per km^2)', 'fontsize',13)
ylabel('Average frequency of all APs blockage','fontsize',13)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,[pwd '/figures/N_block_prob_cond.pdf'],'-dpdf','-r0')
% print(h,[pwd '/figures/N_block_prob_cond.png'],'-dpng','-r0')
