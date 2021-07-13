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
data_moments(data_moments<0)=0;
data_NW_PI=data_moments(1:waves*groups*PI_l);
data_NW_PI=reshape(data_NW_PI,PI_l,waves,groups);

% Model moments
fileID = fopen('model_moments.txt');
model_moments=textscan(fileID,'%f');
fclose(fileID);
model_moments=model_moments{1};
model_moments(isnan(data_moments)) = -9;
model_moments=standardizeMissing(model_moments,-9);
model_moments(model_moments<0) = 0;
model_NW_PI=model_moments(1:waves*groups*PI_l);
model_NW_PI=reshape(model_NW_PI,PI_l,waves,groups);

data_beq_ic=data_moments(waves*groups*PI_l+1:waves*groups*PI_l+f_t*PI_l);
data_beq_ic=reshape(data_beq_ic,PI_l,f_t);
model_beq_ic=model_moments(waves*groups*PI_l+1:waves*groups*PI_l+f_t*PI_l);
model_beq_ic=reshape(model_beq_ic,PI_l,f_t);

data_NW_IC=data_moments(waves*groups*PI_l+f_t*PI_l+1:waves*groups*PI_l+f_t*PI_l+f_t*waves*groups);
data_NW_IC=reshape(data_NW_IC,f_t,waves,groups);
model_NW_IC=model_moments(waves*groups*PI_l+f_t*PI_l+1:waves*groups*PI_l+f_t*PI_l+f_t*waves*groups);
model_NW_IC=reshape(model_NW_IC,f_t,waves,groups);

data_hr_PI=data_moments(waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+1:waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+groups*PI_l);
data_hr_PI=reshape(data_hr_PI,PI_l,groups)
model_hr_PI=model_moments(waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+1:waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+groups*PI_l);
model_hr_PI=reshape(model_hr_PI,PI_l,groups)

data_hr_IC=data_moments(waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+groups*PI_l+1:waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+groups*PI_l+groups*f_t);
data_hr_IC=reshape(data_hr_IC,f_t,groups)
model_hr_IC=model_moments(waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+groups*PI_l+1:waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+groups*PI_l+groups*f_t);
model_hr_IC=reshape(model_hr_IC,f_t,groups)

data_md_IC=data_moments(waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+groups*PI_l+groups*f_t+1:waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+groups*PI_l+groups*f_t+clusters*f_t);
data_md_IC=reshape(data_md_IC,f_t,groups)
model_md_IC=model_moments(waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+groups*PI_l+groups*f_t+1:waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+groups*PI_l+groups*f_t+clusters*f_t);
model_md_IC=reshape(model_md_IC,f_t,groups)

%Untargeted moments
fileID = fopen('data_moments_ut.txt');
data_moments_ut=textscan(fileID,'%f');
fclose(fileID);
data_moments_ut=data_moments_ut{1};
data_moments_ut=standardizeMissing(data_moments_ut,-9);
data_moments_ut(data_moments_ut<0) = 0;
data_moments_ut=reshape(data_moments_ut,f_t,waves,groups)

fileID = fopen('moments_NW_IC1.txt');
model_moments_ut=textscan(fileID,'%f');
fclose(fileID);
model_moments_ut=model_moments_ut{1};
model_moments_ut(isnan(data_moments_ut)) = -9
model_moments_ut=standardizeMissing(model_moments_ut,-9);
model_moments_ut(model_moments_ut<0) = 0;
model_moments_ut=reshape(model_moments_ut,f_t,waves,groups)

fileID = fopen('data_moments_h_ut.txt');
data_moments_ut2=textscan(fileID,'%f');
fclose(fileID);
data_moments_ut2=data_moments_ut2{1};
data_moments_ut2=standardizeMissing(data_moments_ut2,-9);
data_moments_ut2(data_moments_ut2<0) = 0;
data_moments_ut2=reshape(data_moments_ut2,2,f_t,9)

fileID = fopen('model_moments_h_ut.txt');
model_moments_ut2=textscan(fileID,'%f');
fclose(fileID);
model_moments_ut2=model_moments_ut2{1};
model_moments_ut2=standardizeMissing(model_moments_ut2,-9);
model_moments_ut2(model_moments_ut2<0) = 0;
model_moments_ut2=reshape(model_moments_ut2,2,f_t,9)



colors = {[0   0.4470    0.7410]  [0.9290    0.6940    0.1250] [0.4660    0.6740    0.1880] [0.8500    0.3250    0.0980] [0.6350, 0.0780, 0.1840]};
pattern = {'none' 'o' 's' '^' '+'};
width={3 1.5 1.5 1.5 1.5};
waves2=7

% Graph by PI

final_f=figure(1)
set(1,'position',[150 100 700    250])
for m=1:2
for c=1:groups
for p=1:PI_l
    if m==1
        if p>1
            scatter(1+(waves2+1)*(c-1):waves2+(waves2+1)*(c-1),data_NW_PI(p,1:1:waves2,c),pattern{p},'MarkerEdgeColor',colors{p},'MarkerFaceColor',colors{p})
        end
        plot(1+(waves2+1)*(c-1):waves2+(waves2+1)*(c-1),data_NW_PI(p,1:1:waves2,c),'Color',colors{p},'LineWidth',width{p})        
    else
        plot(1+(waves2+1)*(c-1):waves2+(waves2+1)*(c-1),model_NW_PI(p,1:1:waves2,c),':','Color',colors{p},'LineWidth',width{p})
    end
     hold on  
end
end
end
c=1;
p=1;
p1=plot(1+waves2*(c-1):waves2+(waves2)*(c-1),data_NW_PI(p,1:1:waves2,c),'Color',colors{p},'LineWidth',width{p});
p=2;
p2=scatter(1+waves2*(c-1):waves2+(waves2)*(c-1),data_NW_PI(p,1:1:waves2,c),pattern{p},'MarkerEdgeColor',colors{p},'MarkerFaceColor',colors{p});
p=3;
p3=scatter(1+waves2*(c-1):waves2+(waves2)*(c-1),data_NW_PI(p,1:1:waves2,c),pattern{p},'MarkerEdgeColor',colors{p},'MarkerFaceColor',colors{p});
p=4;
p4=scatter(1+waves2*(c-1):waves2+(waves2)*(c-1),data_NW_PI(p,1:1:waves2,c),pattern{p},'MarkerEdgeColor',colors{p},'MarkerFaceColor',colors{p});
p=5;
p5=scatter(1+waves2*(c-1):waves2+(waves2)*(c-1),data_NW_PI(p,1:1:waves2,c),pattern{p},'MarkerEdgeColor',colors{p},'MarkerFaceColor',colors{p});

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
print(gcf,'-depsc', 'C:\Users\jbueren\Google Drive\JMP\Draft\figures\model_fit1.eps')


figure(2)
set(2,'position',[150 100 700    250])
for f_l=1:2
subplot(1,2,f_l)
j=1
for m=1:2
for p=1:2
    if m==1
        if p>1
            scatter(1:6,data_moments_ut2(p,f_l,1:6),pattern{p},'MarkerEdgeColor',colors{p},'MarkerFaceColor',colors{p})
        end
        plot(squeeze(data_moments_ut2(p,f_l,1:6)),'Color',colors{p},'LineWidth',width{p})        
    else
        plot(squeeze(model_moments_ut2(p,f_l,1:6)),':','Color',colors{p},'LineWidth',width{p})
    end
     hold on  
end
end
p=1
p1=plot(squeeze(data_moments_ut2(p,f_l,1:6)),'Color',colors{p},'LineWidth',width{p});
p=2;
p2=scatter(1:6,data_moments_ut2(p,f_l,1:6),pattern{p},'MarkerEdgeColor',colors{p},'MarkerFaceColor',colors{p});
set(gcf,'color','w')
xlabel('Interview','FontSize',FS)
ylabel('Assets (000s of 2018 dollars)','FontSize',FS)
ylim([-9,180])
xticks([1:1:6])
xlim([0 7])
set(gca,'FontName','Times New Roman','Fontsize',FS);
end
I1=legend([p1 p2],'Healthy','In LTC need','location','NorthEast');
legend('boxoff')
I1.FontSize = FS
% I1.FontSize=FS;
% newPosition = [0.41 0.87 0.2 0.2];
% newUnits = 'normalized';
% set(I1,'Position', newPosition,'Units', newUnits);
hold off
for f_l=1:2
    if f_l==1
        subplot(1,2,1)
        title('Distant Family')
    else
        subplot(1,2,2)
        title('Close Family')
    end
end
set(gca,'FontName','Times New Roman','Fontsize',FS);

print(gcf,'-depsc', 'C:\Users\jbueren\Google Drive\JMP\Draft\figures\model_fit3.eps')



figure(3)
set(3,'position',[450    400    700    500])
clrs = [0 0 0; 0.9 0.9 0.9 ];
for m=1:2
    subplot(2,2,m)
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
for m=1:2
    subplot(2,2,m+2)
    if m==1
        hB=bar((data_md_IC*100)')
    else
        hB=bar((model_md_IC*100)')
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
ylim([0 80])
xticks([1 2 3 4])
xticklabels({'Healthy','Physical','Mental','Impaired'})
ylabel('Fraction in Medicaid (%)','FontSize',FS)
xlabel('LTC need','FontSize',FS)
end
print(gcf,'-depsc', 'C:\Users\jbueren\Google Drive\JMP\Draft\figures\model_fit2.eps')

%Table hours matched by PI
          [data_hr_PI([1 3 5],2)';model_hr_PI([1 3 5],2)'; ...
           data_hr_PI([1 3 5],3)';model_hr_PI([1 3 5],3)'; ...
           data_hr_PI([1 3 5],4)';model_hr_PI([1 3 5],4)']
%Table hours matched by PI       
        [data_hr_IC(:,2)';model_hr_IC(:,2)'; ...
           data_hr_IC(:,3)';model_hr_IC(:,3)'; ...
           data_hr_IC(:,4)';model_hr_IC(:,4)']



%% Identification
clear all

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

% Model moments
fileID = fopen('model_moments.txt');
model_moments=textscan(fileID,'%f');
fclose(fileID);
model_moments=model_moments{1};
model_moments=standardizeMissing(model_moments,-9);
model_moments(model_moments<0) = 0;
model_beq_ic=model_moments(waves*groups*PI_l+1:waves*groups*PI_l+f_t*PI_l);
model_beq_ic=reshape(model_beq_ic,PI_l,f_t);
model_hr_IC=model_moments(waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+groups*PI_l+1:waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+groups*PI_l+groups*f_t);
model_hr_IC=reshape(model_hr_IC,f_t,groups)


% Stronger bequest
fileID = fopen('moments_high_beq.txt');
model_moments=textscan(fileID,'%f');
fclose(fileID);
model_moments=model_moments{1};
model_moments=standardizeMissing(model_moments,-9);
model_moments(model_moments<0) = 0;
model_beq_ic_high_beq=model_moments(waves*groups*PI_l+1:waves*groups*PI_l+f_t*PI_l);
model_beq_ic_high_beq=reshape(model_beq_ic_high_beq,PI_l,f_t);
model_hr_IC_high_beq=model_moments(waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+groups*PI_l+1:waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+groups*PI_l+groups*f_t);
model_hr_IC_high_beq=reshape(model_hr_IC_high_beq,f_t,groups)

% Stronger LTC needs
fileID = fopen('moments_high_LTC.txt');
model_moments=textscan(fileID,'%f');
fclose(fileID);
model_moments=model_moments{1};
model_moments=standardizeMissing(model_moments,-9);
model_moments(model_moments<0) = 0;
model_beq_ic_high_ltc=model_moments(waves*groups*PI_l+1:waves*groups*PI_l+f_t*PI_l);
model_beq_ic_high_ltc=reshape(model_beq_ic_high_ltc,PI_l,f_t);
model_hr_IC_high_ltc=model_moments(waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+groups*PI_l+1:waves*groups*PI_l+f_t*PI_l+f_t*waves*groups+groups*PI_l+groups*f_t);
model_hr_IC_high_ltc=reshape(model_hr_IC_high_ltc,f_t,groups)


pattern = {'none' 'o' 's' '^' '+'};
width={3 1.5 1.5 1.5 1.5};
waves2=7
clrs = [0 0 0; 0.5 0.5 0.5 ; 0.9 0.9 0.9 ];

figure(7)
set(7,'position',[450    400    700    250])
p=clusters
hb=bar([model_hr_IC(:,p) model_hr_IC_high_beq(:,p) model_hr_IC_high_ltc(:,p)])  
set(hb,{'FaceColor'},{clrs(1,:),clrs(2,:),clrs(3,:)}.')
xticklabels({'Distant','Close'})
I=legend('Benchmark','High bequest','High LTC','Location','NorthEast')
legend('boxoff')
I.FontSize=FS
ylim([0 7])
hold off
set(gca,'FontName','Times New Roman','Fontsize',FS);
xlabel('Family Type','FontSize',FS)
ylabel('Formal Care Hours per Day','FontSize',FS)
set(gcf,'color','w')
print(gcf,'-depsc', 'C:\Users\jbueren\Google Drive\JMP\Draft\figures\identification1.eps')

figure(8)
set(8,'position',[450    400    700    250])

for m=1:2
    subplot(1,2,m)
    hB=bar([model_beq_ic(:,m)*100 model_beq_ic_high_beq(:,m)*100 model_beq_ic_high_ltc(:,m)*100])
    set(hB,{'FaceColor'},{clrs(1,:),clrs(2,:),clrs(3,:)}.')
    I=legend('Benchmark','High Bequest','High LTC','Location','NorthWest')
    legend('boxoff')
    I.FontSize=FS
if m==1
    title('Distant Family')
else
    title('Close Family')
end
set(gca,'FontName','Times New Roman','Fontsize',FS);
set(gcf,'color','w')
ylim([0 100])
xticks([1 2 3 4 5])
xticklabels({'Bottom','2nd','3rd','4th','Top'})
ylabel('Pr. of leaving any bequest (%)','FontSize',FS)
xlabel('Permanent Income Quintile','FontSize',FS)
end
print(gcf,'-depsc', 'C:\Users\jbueren\Google Drive\JMP\Draft\figures\identification2.eps')


%% Counterfactuals by IC
clear all
generations=21
FS=11;
colors = {[0   0.4470    0.7410] [0.4940    0.1840    0.5560]   [0.9290    0.6940    0.1250]  [0.8500    0.3250    0.0980] [0.3010, 0.7450, 0.9330]};
pattern = {'none' 'o' 's' '^'};
width={3 1.5 3.5 1.5 2};
ldash = {'-' '--' ':'  '*'};

figure(8)
set(8,'position',[50    100    650    350])    
for i=1:2 %family

for p=1:2 % percentile
cd('C:\Users\jbueren\Google Drive\JMP\Code\Structural Model\new_model\new_model\new_model');
fileID = fopen('benchmark.txt');
benchmark=textscan(fileID,'%f  %f %f %f %f %f %f ');
fclose(fileID);
benchmark=benchmark{2+i+(p-1)*2};
counterfactual(:,1)=benchmark(:,1);

fileID = fopen('noLTC.txt');
noLTC=textscan(fileID,'%f %f %f %f %f %f %f ');
fclose(fileID);
noLTC=noLTC{2+i+(p-1)*2};
counterfactual(:,2)=noLTC;

fileID = fopen('noBeq.txt');
noBeq=textscan(fileID,'%f %f %f %f %f %f %f');
fclose(fileID);
noBeq=noBeq{2+i+(p-1)*2};
counterfactual(:,3)=noBeq;

if i==1
    subplot(1,2,1)
    T=title('Distant Families','Position',[85 460])

else
    subplot(1,2,2)
    T=title('Close Families','Position',[85 460])
end
hold on 
for c=1:3
    plot(70:2:100,counterfactual(1:16,c),ldash{c},'Color',colors{c},'LineWidth',width{c})        
end
set(gcf,'color','w')
xlabel('Age','FontSize',FS)
ylabel('Assets (000s of 2018 dollars)','FontSize',FS)
% yticks([0:100:500])
ylim([-9,500])
 if i==1
    I=legend('Benchmark','No long-term care','No bequest','orientation','horizontal');
    legend('boxoff')
    I.FontSize=FS;
    newPosition = [0.42 0.88 0.2 0.2];
    newUnits = 'normalized';
    set(I,'Position', newPosition,'Units', newUnits);
 end 
set(gca,'FontName','Times New Roman','Fontsize',FS);
end
end
print(gcf,'-depsc', 'C:\Users\jbueren\Google Drive\JMP\Draft\figures\counterfactuals_ic.eps')

%% Counterfactuals along no close families
clear all
generations=21
FS=11;
colors = {[0   0.4470    0.7410] [0.4940    0.1840    0.5560]   [0.9290    0.6940    0.1250]  [0.8500    0.3250    0.0980] [0.3010, 0.7450, 0.9330]};
pattern = {'none' 'o' 's' '^'};
width={3 1.5 3.5 1.5 2};
ldash = {'-' '--' ':'  '*'};

figure(7)
set(7,'position',[50    100    600    500])    
subplot(2,1,2)
for i=1:1
cd('C:\Users\jbueren\Google Drive\JMP\Code\Structural Model\new_model\new_model\new_model');
fileID = fopen('benchmark.txt');
benchmark=textscan(fileID,'%f %f %f %f %f %f %f');
fclose(fileID);
benchmark=benchmark{i};
fileID = fopen('noLTC.txt');
noLTC=textscan(fileID,'%f %f %f %f %f %f %f ');
fclose(fileID);
noLTC=noLTC{i};
counterfactual(:,1,:)=benchmark(:,1)-noLTC(:,1);


fileID = fopen('noFemale.txt');
benchmark2=textscan(fileID,'%f %f %f %f %f %f %f');
fclose(fileID);
benchmark2=benchmark2{i};
fileID = fopen('noFemale_noLTC.txt');
noLTC=textscan(fileID,'%f %f %f %f %f %f %f ');
fclose(fileID);
noLTC=noLTC{i};
counterfactual(:,2,:)=benchmark2(:,1)-noLTC(:,1);


for c=1:2
    if c==1
        plot(70:2:90,counterfactual(1:11,c)./benchmark(1:11,1),ldash{c},'Color',colors{c},'LineWidth',width{c})    
    else
        plot(70:2:90,counterfactual(1:11,c)./benchmark2(1:11,1),ldash{c},'Color',colors{c},'LineWidth',width{c}) 
    end
    hold on  
end
set(gcf,'color','w')
xlabel('Age','FontSize',FS)
ylim([-0.01,0.33])
%  yticks([0:0.05:0.3])
xlim([70,90])
%  if i==2
    I=legend('Benchmark','No Close Families','orientation','horizontal');
    legend('boxoff')
    I.FontSize=FS;
    newPosition = [0.4 0.38 0.2 0.2];
    newUnits = 'normalized';
    set(I,'Position', newPosition,'Units', newUnits);
%  end 
set(gca,'FontName','Times New Roman','Fontsize',FS);
end

generations=21
FS=11;
colors = {[0   0.4470    0.7410] [0.4940    0.1840    0.5560]   [0.9290    0.6940    0.1250]  [0.8500    0.3250    0.0980] [0.3010, 0.7450, 0.9330]};
pattern = {'none' 'o' 's' '^'};
width={3 1.5 3.5 1.5 2};
ldash = {'-' '--' ':'  '*'};

subplot(2,1,1)
for i=1:1
cd('C:\Users\jbueren\Google Drive\JMP\Code\Structural Model\new_model\new_model\new_model');

fileID = fopen('benchmark.txt');
benchmark=textscan(fileID,'%f %f %f %f %f %f %f');
fclose(fileID);
benchmark=benchmark{i};

fileID = fopen('noLTC.txt');
noLTC=textscan(fileID,'%f %f %f %f %f %f %f ');
fclose(fileID);
noLTC=noLTC{i};
counterfactual(:,1,:)=benchmark(:,1)-noLTC(:,1);


fileID = fopen('noLTC_physical.txt');
noLTC=textscan(fileID,'%f %f %f %f %f %f %f ');
fclose(fileID);
noLTC=noLTC{i};
counterfactual(:,2,:)=benchmark(:,1)-noLTC(:,1);

fileID = fopen('noLTC_mental.txt');
noLTC=textscan(fileID,'%f %f %f %f %f %f %f ');
fclose(fileID);
noLTC=noLTC{i};
counterfactual(:,3,:)=benchmark(:,1)-noLTC(:,1);

fileID = fopen('noLTC_impaired.txt');
noLTC=textscan(fileID,'%f %f %f %f %f %f %f ');
fclose(fileID);
noLTC=noLTC{i};
counterfactual(:,4,:)=benchmark(:,1)-noLTC(:,1);


for c=1:4
    plot(70:2:90,counterfactual(1:11,c)./benchmark(1:11,1),ldash{c},'Color',colors{c},'LineWidth',width{c})        
    hold on  
end
set(gcf,'color','w')
xlabel('Age','FontSize',FS)
ylim([-0.01,0.33])
%  yticks([0:0.05:0.3])
 xlim([70,90])
%  if i==2
    I=legend('All LTC','No Physical LTC','No Mental LTC','No Impaired LTC','orientation','horizontal');
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
