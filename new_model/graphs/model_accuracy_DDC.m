clear all

%  clc

waves=9;
groups=4;
wealth_q=3
clusters=4
f_t=2;
PI_l=5;
PI_l2=10;
FS=11;
ylim_pi=500
ylim_ic=600
cd('C:\Users\jbueren\Google Drive\JMP\Code\Structural Model\new_model\new_model\new_model');

% Data moments
fileID = fopen('data_moments1.txt');
data_moments=textscan(fileID,'%f');
fclose(fileID);
data_moments=data_moments{1};
data_moments=standardizeMissing(data_moments,-9);
data_moments(data_moments<0) = 0;
data_NW_PI=data_moments(1:waves*groups*PI_l*2);
data_NW_PI=reshape(data_NW_PI,PI_l,waves,groups,2);

% Model moments
fileID = fopen('model_moments.txt');
model_moments=textscan(fileID,'%f');
fclose(fileID);
model_moments=model_moments{1};
model_moments=standardizeMissing(model_moments,-9);
model_moments(model_moments<0) = 0;
model_NW_PI=model_moments(1:waves*groups*PI_l*2);
model_NW_PI=reshape(model_NW_PI,PI_l,waves,groups,2);

data_beq_ic=data_moments(waves*groups*PI_l*2+1:waves*groups*PI_l*2+f_t*PI_l);
data_beq_ic=reshape(data_beq_ic,PI_l,f_t);
model_beq_ic=model_moments(waves*groups*PI_l*2+1:waves*groups*PI_l*2+f_t*PI_l);
model_beq_ic=reshape(model_beq_ic,PI_l,f_t);

data_NW_IC=data_moments(waves*groups*PI_l*2+f_t*PI_l+1:waves*groups*PI_l*2+f_t*PI_l+f_t*waves*groups);
data_NW_IC=reshape(data_NW_IC,f_t,waves,groups);

model_NW_IC=model_moments(waves*groups*PI_l*2+f_t*PI_l+1:waves*groups*PI_l*2+f_t*PI_l+f_t*waves*groups);
model_NW_IC=reshape(model_NW_IC,f_t,waves,groups);

data_hr_PI=data_moments(waves*groups*PI_l*2+f_t*PI_l+f_t*waves*groups+1:waves*groups*PI_l*2+f_t*PI_l+f_t*waves*groups+groups*PI_l);
data_hr_PI=reshape(data_hr_PI,PI_l,groups)

model_hr_PI=model_moments(waves*groups*PI_l*2+f_t*PI_l+f_t*waves*groups+1:waves*groups*PI_l*2+f_t*PI_l+f_t*waves*groups+groups*PI_l);
model_hr_PI=reshape(model_hr_PI,PI_l,groups)

data_hr_IC=data_moments(waves*groups*PI_l*2+f_t*PI_l+f_t*waves*groups+groups*PI_l+1:waves*groups*PI_l*2+f_t*PI_l+f_t*waves*groups+groups*PI_l+groups*f_t);
data_hr_IC=reshape(data_hr_IC,f_t,groups)

model_hr_IC=model_moments(waves*groups*PI_l*2+f_t*PI_l+f_t*waves*groups+groups*PI_l+1:waves*groups*PI_l*2+f_t*PI_l+f_t*waves*groups+groups*PI_l+groups*f_t);
model_hr_IC=reshape(model_hr_IC,f_t,groups)



colors = {[0   0.4470    0.7410]  [0.9290    0.6940    0.1250] [0.4660    0.6740    0.1880] [0.8500    0.3250    0.0980] [0.6350, 0.0780, 0.1840]};
pattern = {'none' 'o' 's' '^' '+'};
width={3 1.5 1.5 1.5 1.5};
waves2=7

% Graph by PI

final_f=figure(3)
set(3,'position',[50    20    500    750])
j=1
h1=subplot(3,1,1)
for m=1:2
for c=1:groups
for p=1:PI_l
    if m==1
        if p>1
            scatter(1+(waves2+1)*(c-1):waves2+(waves2+1)*(c-1),data_NW_PI(p,1:1:waves2,c,j),pattern{p},'MarkerEdgeColor',colors{p},'MarkerFaceColor',colors{p})
        end
        plot(1+(waves2+1)*(c-1):waves2+(waves2+1)*(c-1),data_NW_PI(p,1:1:waves2,c,j),'Color',colors{p},'LineWidth',width{p})        
    else
        plot(1+(waves2+1)*(c-1):waves2+(waves2+1)*(c-1),model_NW_PI(p,1:1:waves2,c,j),':','Color',colors{p},'LineWidth',width{p})
    end
     hold on  
end
end
end
c=1;
p=1;
p1=plot(1+waves2*(c-1):waves2+(waves2)*(c-1),data_NW_PI(p,1:1:waves2,c,j),'Color',colors{p},'LineWidth',width{p});
p=2;
p2=scatter(1+waves2*(c-1):waves2+(waves2)*(c-1),data_NW_PI(p,1:1:waves2,c,j),pattern{p},'MarkerEdgeColor',colors{p},'MarkerFaceColor',colors{p});
p=3;
p3=scatter(1+waves2*(c-1):waves2+(waves2)*(c-1),data_NW_PI(p,1:1:waves2,c,j),pattern{p},'MarkerEdgeColor',colors{p},'MarkerFaceColor',colors{p});
p=4;
p4=scatter(1+waves2*(c-1):waves2+(waves2)*(c-1),data_NW_PI(p,1:1:waves2,c,j),pattern{p},'MarkerEdgeColor',colors{p},'MarkerFaceColor',colors{p});
p=5;
p5=scatter(1+waves2*(c-1):waves2+(waves2)*(c-1),data_NW_PI(p,1:1:waves2,c,j),pattern{p},'MarkerEdgeColor',colors{p},'MarkerFaceColor',colors{p});

set(gcf,'color','w')
xlabel('Age','FontSize',FS)
ylabel('Assets (000s of 2018 dollars)','FontSize',FS)
ylim([-9,ylim_pi])
xticks([1:2:32])
xticklabels({'72','76','80','84',...
             '76','80','84','88',...
             '82','86','90','94',...
             '86','90','94','98'})
xlim([0 32])
alpha=0.6;
line([8 8], [-9 ylim_pi], 'Color', [alpha alpha alpha]);
line([16 16], [-9 ylim_pi], 'Color', [alpha alpha alpha]);
line([24 24], [-9 ylim_pi], 'Color', [alpha alpha alpha]);
text(2.4,ylim_pi-50,'Group 1','FontName','Times New Roman','Fontsize',FS)
text(10.4,ylim_pi-50,'Group 2','FontName','Times New Roman','Fontsize',FS)
text(18.4,ylim_pi-50,'Group 3','FontName','Times New Roman','Fontsize',FS)
text(26.4,ylim_pi-50,'Group 4','FontName','Times New Roman','Fontsize',FS)
I1=legend([p1 p2 p3 p4 p5],'Bottom','Second','Third','Fourth','Top','orientation','horizontal');
legend('boxoff')
I1.FontSize=FS;
newPosition = [0.4 0.87 0.2 0.2];
newUnits = 'normalized';
set(I1,'Position', newPosition,'Units', newUnits);
hold off
set(gca,'FontName','Times New Roman','Fontsize',FS);

h2=subplot(3,1,2)
for m=1:2
for p=2:clusters
    if m==1
        scatter(1:PI_l,data_hr_PI(:,p),pattern{p},'MarkerEdgeColor',colors{p},'MarkerFaceColor',colors{p})                
    else
        plot(1:PI_l,data_hr_PI(:,p),'Color',colors{p},'LineWidth',width{p})
        plot(1:PI_l,model_hr_PI(:,p),':','Color',colors{p},'LineWidth',width{p})
    end
     hold on  
end
end
xticks([1:1:5])
xticklabels({'Bottom','Second','Third','Fourth','Top'})
ylim([0 9])
xlim([0.8 5.2])
I2=legend('Physically Frail','Mentally Frail','Impaired','orientation','horizontal');
legend('boxoff')
I2.FontSize=FS;
newPosition = [0.4 0.55 0.2 0.2];
newUnits = 'normalized';
set(I2,'Position', newPosition,'Units', newUnits);
hold off
set(gca,'FontName','Times New Roman','Fontsize',FS);
xlabel('Permanent Income Quintile','FontSize',FS)
ylabel('Formal Care Hours per Day','FontSize',FS)

h3=subplot(3,1,3)
for m=1:2
for p=2:clusters
    if m==1
        scatter(1:f_t,data_hr_IC(:,p),pattern{p},'MarkerEdgeColor',colors{p},'MarkerFaceColor',colors{p})                
    else
        plot(1:f_t,data_hr_IC(:,p),'Color',colors{p},'LineWidth',width{p})
        plot(1:f_t,model_hr_IC(:,p),':','Color',colors{p},'LineWidth',width{p})
    end
     hold on  
end
end
xticks([1:1:2])
xticklabels({'On your own','Close'})
ylim([0 9])
xlim([0.8 2.2])
I3=legend('Physically Frail','Mentally Frail','Impaired','orientation','horizontal');
legend('boxoff')
I3.FontSize=FS;
newPosition = [0.4 0.25 0.2 0.2];
newUnits = 'normalized';
set(I3,'Position', newPosition,'Units', newUnits);
hold off
set(gca,'FontName','Times New Roman','Fontsize',FS);
xlabel('Family Type','FontSize',FS)
ylabel('Formal Care Hours per Day','FontSize',FS)

figure(4)
set(4,'position',[450    400    700    250])
clrs = [0 0 0; 0.9 0.9 0.9 ];
for m=1:2
    subplot(1,2,m)
    if m==1
        hB=bar(data_beq_ic*100)
    else
        hB=bar(model_beq_ic*100)
    end
set(hB,{'FaceColor'},{clrs(1,:),clrs(2,:)}.')

    I=legend('Distant','Close','Location','NorthWest')
    legend('boxoff')
    I.FontSize=FS
%     newPosition = [0.5 0.925 0.08 0.08];
%     newUnits = 'normalized';
%     set(I,'Position', newPosition,'Units', newUnits);
if m==1
    title('Data')
else
    title('Model')
end
set(gca,'FontName','Times New Roman','Fontsize',FS);
set(gcf,'color','w')
ylim([0 100])
xticks([1 2 3 4 5])
xticklabels({'Bottom','2nd','3rd','4th','Top'})
ylabel('Pr. of leaving any bequest (%)','FontSize',FS)
xlabel('Permanent Income Quintile','FontSize',FS)
end



%% Untargetted moments

fileID = fopen('data_moments_ut.txt');
data_moments_ut=textscan(fileID,'%f');
fclose(fileID);
data_moments_ut=data_moments_ut{1};
data_moments_ut=standardizeMissing(data_moments_ut,-9);
data_moments_ut(data_moments_ut<0) = 0;
data_moments_ut=reshape(data_moments_ut,2,waves);

fileID = fopen('model_NW_h_ut.txt');
model_moments_ut=textscan(fileID,'%f');
fclose(fileID);
model_moments_ut=model_moments_ut{1};
model_moments_ut=standardizeMissing(model_moments_ut,-9);
model_moments_ut(model_moments_ut<0) = 0;
model_moments_ut=reshape(model_moments_ut,2,waves);

figure(6)
set(6,'position',[50    100    500    250])    

    plot(1:7,data_moments_ut(1,1:7),'Color',colors{1},'LineWidth',width{1})
    hold on
    plot(1:7,data_moments_ut(2,1:7),'Color',colors{2},'LineWidth',width{2})
    plot(1:7,model_moments_ut(1,1:7),':','Color',colors{1},'LineWidth',width{1})
    plot(1:7,model_moments_ut(2,1:7),':','Color',colors{2},'LineWidth',width{2})

I2=legend('Healthy','In need of LTC (h>1)','orientation','horizontal');
legend('boxoff')
I2.FontSize=FS;
newPosition = [0.42 0.87 0.2 0.2];
newUnits = 'normalized';
set(I2,'Position', newPosition,'Units', newUnits);
hold off
set(gca,'FontName','Times New Roman','Fontsize',FS);
xlabel('Interview number','FontSize',FS)
ylabel('Assets (000s of 2018 dollars)','FontSize',FS)
set(gcf,'color','w')
ylim([-5 180])
xlim([0.8 7.2])
print(gcf,'-depsc', 'C:\Users\jbueren\Google Drive\JMP\Draft\figures\model_fit_ut.eps')


%% Counterfactuals along wealth distribution
clear all
generations=21
FS=11;
colors = {[0   0.4470    0.7410] [0.4940    0.1840    0.5560]   [0.9290    0.6940    0.1250]  [0.8500    0.3250    0.0980] [0.3010, 0.7450, 0.9330]};
pattern = {'none' 'o' 's' '^'};
width={3 1.5 3.5 1.5 2};
ldash = {'-' '--' ':'  '*'};

figure(6)
set(6,'position',[50    100    600    250])    
for i=1:2
cd('C:\Users\jbueren\Google Drive\JMP\Code\Structural Model\new_model\new_model\new_model');
fileID = fopen('benchmark.txt');
benchmark=textscan(fileID,'%f  %f %f %f %f %f ');
fclose(fileID);
benchmark=benchmark{i};
counterfactual(:,1,:)=benchmark(:,1);

fileID = fopen('noLTC.txt');
noLTC=textscan(fileID,'%f   %f %f %f %f %f ');
fclose(fileID);
noLTC=noLTC{i};
counterfactual(:,2,:)=noLTC;

fileID = fopen('noBeq.txt');
noBeq=textscan(fileID,'%f   %f %f %f %f %f');
fclose(fileID);
noBeq=noBeq{i};
counterfactual(:,3,:)=noBeq;

fileID = fopen('noMed.txt');
noMed=textscan(fileID,'%f  %f %f %f %f %f');
fclose(fileID);
noMed=noMed{i};
counterfactual(:,4,:)=noMed;

for c=1:4
    plot(70:2:100,counterfactual(1:16,c),ldash{c},'Color',colors{c},'LineWidth',width{c})        
    hold on  
end
set(gcf,'color','w')
xlabel('Age','FontSize',FS)
ylabel('Assets (000s of 2018 dollars)','FontSize',FS)
 ylim([-9,500])
%  if i==2
    I=legend('Benchmark','No long-term care','No bequest','No medical expenses','orientation','horizontal');
    legend('boxoff')
    I.FontSize=FS;
    newPosition = [0.4 0.88 0.2 0.2];
    newUnits = 'normalized';
    set(I,'Position', newPosition,'Units', newUnits);
%  end 
set(gca,'FontName','Times New Roman','Fontsize',FS);
end
print(gcf,'-depsc', 'C:\Users\jbueren\Google Drive\JMP\Draft\figures\counterfactuals.eps')
print(gcf,'-depsc', 'C:\Users\jbueren\Google Drive\JMP\Slides\figures\counterfactuals.eps')

%% Experiment: break correlation health and death
clear all
generations=21
FS=11;
colors = {[0   0.4470    0.7410] [0.4940    0.1840    0.5560]   [0.9290    0.6940    0.1250]  [0.8500    0.3250    0.0980] [0.3010, 0.7450, 0.9330]};
pattern = {'none' 'o' 's' '^'};
width={3 1.5 3.5 1.5 2};
ldash = {'-' '--' ':'  '*'};

figure(6)
set(6,'position',[50    100    500    250])    
for i=1:2
cd('C:\Users\jbueren\Google Drive\JMP\Code\Structural Model\DDC');
fileID = fopen('benchmark_h.txt');
benchmark=textscan(fileID,'%f  %f %f %f %f %f ');
fclose(fileID);
benchmark=benchmark{i};
counterfactual(:,1,:)=benchmark(:,1);

fileID = fopen('health_exp.txt');
noLTC=textscan(fileID,'%f   %f %f %f %f %f ');
fclose(fileID);
noLTC=noLTC{i};
counterfactual(:,2,:)=noLTC;

for c=1:2
    plot(70:2:100,counterfactual(1:16,c),ldash{c},'Color',colors{c},'LineWidth',width{c})        
    hold on  
end
set(gcf,'color','w')
xlabel('Age','FontSize',FS)
ylabel('Assets (000s of 2018 dollars)','FontSize',FS)
 ylim([-9,500])
%  if i==2
    I=legend('Benchmark','Experiment','orientation','horizontal');
    legend('boxoff')
    I.FontSize=FS;
    newPosition = [0.42 0.89 0.2 0.2];
    newUnits = 'normalized';
    set(I,'Position', newPosition,'Units', newUnits);
%  end 
set(gca,'FontName','Times New Roman','Fontsize',FS);
end
print(gcf,'-depsc', 'C:\Users\jbueren\Google Drive\JMP\Draft\figures\exp_health.eps')


%% Counterfactuals by IC
clear all
generations=21
FS=11;
colors = {[0   0.4470    0.7410] [0.4940    0.1840    0.5560]   [0.9290    0.6940    0.1250]  [0.8500    0.3250    0.0980] [0.3010, 0.7450, 0.9330]};
pattern = {'none' 'o' 's' '^'};
width={3 1.5 3.5 1.5 2};
ldash = {'-' '--' ':'  '*'};

figure(7)
set(7,'position',[50    100    750    400])    
for i=1:2:3
cd('C:\Users\jbueren\Google Drive\JMP\Code\Structural Model\DDC');
fileID = fopen('benchmark.txt');
benchmark=textscan(fileID,'%f  %f %f %f %f %f ');
fclose(fileID);
benchmark=benchmark{2+i};
counterfactual(:,1,:)=benchmark(:,1);

fileID = fopen('noLTC.txt');
noLTC=textscan(fileID,'%f   %f %f %f %f %f ');
fclose(fileID);
noLTC=noLTC{2+i};
counterfactual(:,2,:)=noLTC;

fileID = fopen('noBeq.txt');
noBeq=textscan(fileID,'%f   %f %f %f %f %f');
fclose(fileID);
noBeq=noBeq{2+i};
counterfactual(:,3,:)=noBeq;

if i==1
    subplot(1,2,1)
    title('On your own')
else
    subplot(1,2,2)
    title('Close Families')
end
hold on 
for c=1:3
    plot(70:2:100,counterfactual(1:16,c),ldash{c},'Color',colors{c},'LineWidth',width{c})        
     
end
set(gcf,'color','w')
xlabel('Age','FontSize',FS)
ylabel('Assets (000s of 2018 dollars)','FontSize',FS)
%  ylim([-9,500])
 if i==1
    I=legend('Benchmark','No long-term care','No bequest','orientation','horizontal');
    legend('boxoff')
    I.FontSize=FS;
    newPosition = [0.42 0.89 0.2 0.2];
    newUnits = 'normalized';
    set(I,'Position', newPosition,'Units', newUnits);
 end 
set(gca,'FontName','Times New Roman','Fontsize',FS);
end
print(gcf,'-depsc', 'C:\Users\jbueren\Google Drive\JMP\Draft\figures\counterfactuals_ic.eps')
