%% Intensity of bequest, consumption and formal care hours across wealth, family and health status.
clear all
clc   
cd('C:\Users\jbueren\Google Drive\JMP\Code\Structural Model\new_model\new_model\new_model');
% Data moments
fileID = fopen('parameters_prelim.txt');
parameters=textscan(fileID,'%f');
fclose(fileID);
parameters=parameters{1};
%Benchmark
sigma(1)=parameters(1);
nu(1)=parameters(7);
sigma_2(1)=sigma(1);
mu(1,1)=0;
mu(2:4,1)=exp(parameters(3:5));
lambda(1:2,1)=exp(parameters(8:9))*2;
kappa(1)=parameters(6);
om=1;
share_p=parameters(10);
subs_p=parameters(11);

p_fc(1)=8/1000;
theta(1:4,1)=1;

% DFJ
    %2010
    kappa(2)=13.0;
    sigma(2)=2.83;
    sigma_2(2)=sigma(2);
    mu(1:4,2)=0;
    lambda(1:3,2)=39.7;
    theta(1,2)=1;
    theta(2:4,2)=1-0.36;
    p_fc(2)=1/1000;
    
    %2016
%     kappa(2)=13.0;
%     sigma(2)=2.83;
%     mu(1,2)=0;
%     mu(2:4,2)=-13.57+0.2589*85-0.3035*85^2/100-0.2050*85^3/10000-0.2823+0.00653*85;
%     theta(1:4,2)=1;
%     lambda(1:3,2)=39.7;    

% Amerik et al. (WP)
    %Joint Estimation (Asset & SSQ)
    kappa(3)=7.83;
    sigma(3)=5.27;
    sigma_2(3)=sigma(3);
    lambda(1:3,3)=0.67;
    mu(1:4,3)=0;
    theta(1,3)=1;
    theta(2:4,3)=1.09;
    kappa_ltc=-37.44
    %SSQ only
%     kappa(3)=0.34;
%     sigma(3)=3.56;
%     lambda(1:3,3)=0.18;
%     mu(1:4,3)=0;
%     theta(1,3)=1;
%     theta(2:4,3)=0.13;
%     kappa_ltc=-63.97;
    p_fc(3)=0

% Lockwood (2017)
kappa(4)=231.0;
sigma(4)=3.0;
sigma_2(4)=sigma(4);
lambda(1:3,4)=1470;
mu(1:4,4)=0;
theta(1:4,4)=1;
p_fc(4)=0



%Grid of consumption
N_x=500;
a_max=500
x_grid=linspace(0,a_max,N_x);

l_IC(1:3,1:4)=zeros(3,4);
for m_l=1:4
    m_l
    if m_l==1
        %Informal care on their own
        l_IC(1,2:4)=[0.5*365*2+68 0.5*365*2+521 0.5*365*2]; 
        %Informal care close families
        l_IC(2,2:4)=[3.6*365*2+68 8.9*365*2+521 13.6*365*2];
    end
for h_l=1:4;
    h_l
for IC_l=1:2;
for i=1:N_x;
    i;
    u_j(1:i)=-1/0;
    c_j(1:i)=-1/0;
    fc_j(1:i)=-1/0;
    for j=1:i
        a2=x_grid(j);
        cts=x_grid(i)-x_grid(j); %cash to spend
        c_max=cts;
        if mu(h_l,m_l)==0 
            c_try=c_max;
            fc_try=0;
        else
            %Interior solution
            l_fc_min=0;
            l_fc_max=cts/p_fc(m_l);
            crit=2;
            while abs(crit)>0.0000001
                fc_try=(l_fc_min+l_fc_max)/2;
                c_try=(mu(h_l,m_l)/p_fc(m_l))^(-1/sigma(m_l))*share_p^(-1/sigma(m_l))*fc_try^((1-subs_p)/sigma(m_l))*...
                      (share_p*fc_try^subs_p+(1-share_p)*l_IC(IC_l,h_l)^subs_p)^((nu(m_l)+subs_p-1)/(subs_p*sigma(m_l)));
                LHS=c_try+p_fc(m_l)*fc_try;
                crit=LHS-cts;
                if cts>LHS
                    l_fc_min=fc_try;
                else
                    l_fc_max=fc_try;
                end 
            end
        end
        c_j(j)=c_try;
        fc_j(j)=fc_try*p_fc(m_l);
        if m_l==3 && h_l>1 && (c_try>-kappa_ltc)
            u_j(j)=theta(h_l,m_l)*(c_try+kappa_ltc)^(1-sigma(m_l))/(1-sigma(m_l))+lambda(IC_l,m_l)*(a2+kappa(m_l))^(1-sigma_2(m_l))/(1-sigma_2(m_l));
        elseif m_l==3 && h_l>1 && (c_try<-kappa_ltc)
            u_j(j)=-999999;
        elseif mu(h_l,m_l)==0
            u_j(j)=theta(h_l,m_l)*c_try^(1-sigma(m_l))/(1-sigma(m_l))+lambda(IC_l,m_l)*(a2+kappa(m_l))^(1-sigma_2(m_l))/(1-sigma_2(m_l));
        elseif m_l==1
            u_j(j)=theta(h_l,m_l)*c_try^(1-sigma(m_l))/(1-sigma(m_l))+mu(h_l,m_l)*(share_p*fc_try^subs_p+(1-share_p)*l_IC(IC_l,h_l)^subs_p)^((1-nu(m_l))/subs_p)/(1-nu(m_l))+lambda(IC_l,m_l)*(a2+kappa(m_l))^(1-sigma_2(m_l))/(1-sigma_2(m_l));
        else
            u_j(j)=theta(h_l,m_l)*c_try^(1-sigma(m_l))/(1-sigma(m_l))+mu(h_l,m_l)*(fc_try+om*l_IC(IC_l,h_l))^(1-nu(m_l))/(1-nu(m_l))+lambda(IC_l,m_l)*(a2+kappa(m_l))^(1-sigma_2(m_l))/(1-sigma_2(m_l));
        end            
    end
    [val loc]=max(u_j(1:i));
    if val==-999999
        loc=1;
    end
    c(m_l,h_l,IC_l,i)=c_j(loc);
    fc(m_l,h_l,IC_l,i)=fc_j(loc);
    beq(m_l,h_l,IC_l,i)=x_grid(loc);
end
end
end
end 

% My model

colors = {[0    0.4470    0.7410] [0.8500    0.3250    0.0980] [0.9290    0.6940    0.1250] [0.4940    0.1840    0.5560] [0.4660    0.6740    0.1880]};
pattern = {'none' 'x' 's' '^'};
width={4 0.6+2*(4-0.6)/3 0.6+(4-0.6)/3 0.6}
pattern2 = {'-' ':' ':'};
m_l=1
gap=20
figure(1)
for h_l=[1 4]
    for IC_l=[   1 2 ]
        pos=1
        if h_l==4
            pos=2
        end
        subplot(1,2,pos)
        plot(x_grid(1:gap:end),squeeze(beq(m_l,h_l,IC_l,1:gap:end))./x_grid(1:gap:end)','Color',colors{1},'Marker',pattern{1},'MarkerSize',6,'LineWidth',width{1},'LineStyle',pattern2{IC_l})
        hold on
        if h_l>1
            plot(x_grid(1:gap:end),squeeze(fc(m_l,h_l,IC_l,1:gap:end))./x_grid(1:gap:end)','Color',colors{2},'Marker',pattern{2},'MarkerSize',6,'LineWidth',width{2},'LineStyle',pattern2{IC_l})
        end
        plot(x_grid(1:gap:end),squeeze(c(m_l,h_l,IC_l,1:gap:end))./x_grid(1:gap:end)','Color',colors{3},'Marker',pattern{3},'MarkerSize',6,'LineWidth',width{3},'LineStyle',pattern2{IC_l})
        ylim([-0.05 1.05])
        if h_l==4
            legend('Bequest','LTC','Consumption')
        end
    end
end

%% vs others: good health

colors = {[0    0.4470    0.7410] [0.8500    0.3250    0.0980] [0.9290    0.6940    0.1250] [0.4940    0.1840    0.5560] [0.4660    0.6740    0.1880]};
pattern = {'none' 'o' 's' '^' 'x'};
width={4 0.6+2*(4-0.6)/3 0.6+(4-0.6)/3 0.6}
pattern2 = {'-' '-.' ':' '--'};
step=80

ms=90
FS=11

    h_l=1
    figure(h_l+3)
    set(h_l+3,'position',[600    500    650    400*0.75])
    g=subplot(1,2,1)
    % DFJ (2016) & Ameriks et al. (WP) & Lockwood (WP)
    IC_l=1
    p1=plot(x_grid(1:end),squeeze(beq(2,h_l,IC_l,1:end))./x_grid(1:end)','Color',colors{1},'LineWidth',width{3},'LineStyle',pattern2{3})
    p2=plot(x_grid(1:end),squeeze(beq(4,h_l,IC_l,1:end))./x_grid(1:end)','Color',colors{1},'LineWidth',width{3},'LineStyle',pattern2{4})
    hold on
    p3=plot(x_grid(1:end),squeeze(beq(3,h_l,IC_l,1:end))./x_grid(1:end)','-o','Color',colors{1},'LineWidth',width{3},'LineStyle',pattern2{1},'MarkerIndices',1:50:N_x)
    p4=plot(x_grid(1:end),squeeze(beq(3,2,IC_l,1:end))./x_grid(1:end)',':','Color',colors{1},'LineWidth',width{3},'LineStyle',pattern2{1},'MarkerIndices',1:50:N_x)
    
    ylim([-0.03 1.38])
    xlim([0 a_max])
    I=legend([p2,p3,p4],'Lockwood','SSQ, healthty','SSQ, LTC','location','NorthWest');
    legend('boxoff');
    set(I,'FontName','Times New Roman','FontSize',FS);
    ylabel('Allocation','FontSize',FS)
    xlabel('Wealth')
    set(gca,'FontName','Times New Roman','FontSize',FS);
%     newPosition = [0.13 0.87 0.2 0.2];
%     newUnits = 'normalized';
%     set(I,'Position', newPosition,'Units', newUnits);
    g=subplot(1,2,2)
    % My model
    IC_l=1
    p1=plot(x_grid(1:10:end),squeeze(beq(1,1,IC_l,1:10:end))./x_grid(1:10:end)','Color',colors{1},'LineWidth',width{3},'LineStyle',pattern2{1})
    hold on
%     plot(x_grid(1:end),squeeze(beq(1,2,IC_l,1:end))./x_grid(1:end)','Color',colors{1},'LineWidth',width{3},'LineStyle',pattern2{1})
    p2=scatter(x_grid(1:step:end),squeeze(beq(1,2,IC_l,1:step:end))./x_grid(1:step:end)',ms-60,pattern{2},'MarkerEdgeColor',colors{1},'MarkerFaceColor',colors{1})
    p3=plot(x_grid(1:end),squeeze(beq(1,3,IC_l,1:end))./x_grid(1:end)','Color',colors{1},'LineWidth',width{3},'LineStyle',pattern2{4})
    p4=plot(x_grid(1:end),squeeze(beq(1,4,IC_l,1:end))./x_grid(1:end)','Color',colors{1},'LineWidth',width{3},'LineStyle',pattern2{2})
%     plot(100,0.0,'ro')
%     plot(150,0.0,'ro')
%     plot(200,0.5,'ro')
    ylim([-0.03 1.38])
    xlim([0 a_max])
    I=legend([p1,p2,p3,p4],'Healthy','Physically Frail','Mentally Frail','Impaired','location','NorthWest');
    legend('boxoff');
    set(I,'FontName','Times New Roman','FontSize',FS);
%     newPosition = [0.63 0.87 0.2 0.2];
%     newUnits = 'normalized';
%     set(I,'Position', newPosition,'Units', newUnits);
    xlabel('Wealth')
    set(gcf,'color','w')
    set(gca,'FontName','Times New Roman','FontSize',FS);

    print(gcf,'-depsc', 'C:\Users\jbueren\Google Drive\JMP\Draft\figures\spending_death_others_gh.eps')
    print(gcf,'-depsc', 'C:\Users\jbueren\Google Drive\JMP\Slides\figures\spending_death_others_gh.eps')

        
% end

%% Share devoted to fc from 20000 spending

m_l=1

cts=100000; %cash to spend
for h_l=1:4
    c_max=cts;
if mu(h_l,m_l)==0 
    c_try=c_max;
    fc_try=0;
else
    %Interior solution
        l_fc_min=0;
        l_fc_max=cts/p_fc(m_l);
        crit=2;
        while abs(crit)>0.0000001
            fc_try=(l_fc_min+l_fc_max)/2;
            c_try=(mu(h_l,m_l)/p_fc(m_l))^(-1/sigma(m_l))*share_p^(-1/sigma(m_l))*fc_try^((1-subs_p)/sigma(m_l))*...
                  (share_p*fc_try^subs_p+(1-share_p)*l_IC(IC_l,h_l)^subs_p)^((nu(m_l)+subs_p-1)/(subs_p*sigma(m_l)));
            LHS=c_try+p_fc(m_l)*fc_try;
            crit=LHS-cts;
            if cts>LHS
                l_fc_min=fc_try;
            else
                l_fc_max=fc_try;
            end 
        end
end
c_j(h_l)=c_try;
fc_j(h_l)=fc_try*p_fc(m_l);
end
fc_j(1:h_l)/cts





