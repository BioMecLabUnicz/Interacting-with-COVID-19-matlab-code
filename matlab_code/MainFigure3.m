% The script generates the results shown in Figure 3 of the ms. 
% In the script you set the parameter id_country in order to show the
% result for a given country:
% id_country=1, plot the results for Germany, 2 for France, 3 for Italy, 4
% for UK


clear all
close all
clc
global  data_dQR data_QR  n_pop  eta alpha H_0 time1

global beta K_p tau tag_change  

global m w n_h 

global time2  w0 w1  w2 w3 time3 u0 u1 u2 u3

global tag_country

id_country=3; % id for country 1=Germany, 2=France, 3=Italy


% load the fitting results for SIHM2
filename_f=['res_opt/num_country_',num2str(id_country),'_data_dh_time_changes_2_fmincon'];
load(filename_f)
% 
xopt=xfmincon;

if id_country==3 ||  id_country==4
    eta=0.006;
    alpha=0.094;
elseif id_country==1 ||  id_country==2    
    eta=0.006*7;
    alpha=0.094*7;
end

I0=xopt(1)/eta;
beta=xopt(2);
K_p=xopt(3);
n_h=xopt(4);
tau=xopt(5);
time1=xopt(6);
w1=xopt(7);
time2=xopt(8);
w2=xopt(9);        
w0=0;
u0=0;              


%%%%%%%%%%%%%%%%%%%%%
%%%%%% load  data  %%%%%%
%%%%%%%%%%%%%%%%%%%%%

if id_country==1
    run data_dH_Germany
    n_pop = 83e6; %Germany
    
elseif id_country==2
    run data_dH_France.m
    n_pop = 67.1e6; %France
    H_0=H_Fr_0;
elseif id_country==3
    run data_dH_Italy
    n_pop = 60e6; %Italy
    
elseif id_country==4
    run data_dH_UK
    n_pop = 66.6e6; %UK
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% add data to march-april 2021 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if id_country==1
    
    data_dH=[num_dH_Ger;num_dH_Ger_2021];
    H_0=100;
    data_H= cumsum([H_0; data_dH]);
    data_dQR=data_dH(1:60);
    data_QR=data_H(1:61);
    
    time_x=40;
    time3=52.5; % time (week) when alpha variant became dominant
    time_pr=46;      % inital time (week) for prediction simulation
elseif id_country==2

    data_dH=data_dH_week_11_new;
    data_H= cumsum([H_0; data_dH]); % H_0=250, end week 9 H=116, end week 11=400

    data_dQR=data_dH(1:59);
    data_QR=data_H(1:60);

    time_x=39;

    time3=52; % time (week) when alpha variant became dominant
    time_pr=45;      % inital time (week) for prediction simulation

    
elseif id_country==3
    
    data_dH=[data_dH3;data_dH_2021];
    
    H_0=0;
    data_cumsumH=[H_0; cumsum(data_dH)];
    
    data_dQR=data_dH;
    data_QR=data_cumsumH;
    
    time_x=300;
    time3=366;
    time_pr=333;
elseif id_country==4
       return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% introduce the third change for w %%%%%
%%%% for simulating alpha variant effect %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tag_change=3;

% compute w3 (i.e. alpha variant effect)
w3_pred=1.2*(1+xfmincon(9))-1;

w3=w3_pred;
xfmincon1=[xfmincon time3 w3_pred];

size_time_window=length(data_dQR);

if id_country==1 ||id_country==2
    
    tspan1=0:1/7:size_time_window;
    
    idx_tx=find(tspan1>=time_x);
    idx_tpr=find(tspan1>=time_pr);
    idx_t3=find(tspan1>=time3);
   
    if id_country==1
        tx=datetime(2020,3,1) + calweeks(0:size_time_window); 
    else
        tx=datetime(2020,3,8) + calweeks(0:size_time_window);
    end
    
    Tx=tx(1):tx(end);
    
elseif id_country==3
    
    tspan1=0:1/4:size_time_window;


    tx=datetime(2020,2,19) + caldays(0:size_time_window);

    Tx=tx(1):1/4:tx(end);
    
    idx_tx=find(tspan1>=time_x);

    idx_t3=find(tspan1>=time3);
    
    idx_tpr=find(tspan1>=time_pr);
end



u1=0;u2=0;
u3=0; %

I0=xfmincon1(1)/eta;
Y0=zeros(4,1);
Y0(1,1)=n_pop-I0-H_0;
Y0(2,1)=I0;
Y0(3,1)=H_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% simulate the model prediction %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[T0, Ys]=ode45('covid_model_SIHM_distu_opt',tspan1, Y0,[]);

I_pr=Ys(:,2);
data_sim_pr=eta*I_pr;


%%%% set an additional control action to maintain the peak below a
%%%% threshold
if id_country==1
    u3=0.094; % peak 6300
elseif id_country==2
    u3=0.0644; % peak 11k
elseif id_country==3
    u3=0.0419; % peak 1.5k
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% simulate the effects of the additional %%%%%
%%%%%%%%%%% control action %%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[T1, Y1]=ode45('covid_model_SIHM_distu_opt',tspan1, Y0,[]);

I_sim_u3_a=Y1(:,2);
data_sim_u3_a=eta*I_sim_u3_a;


%%%% set another control action to maintain the peak 
%%%% below a threshold
if id_country==1
    u3=0.15; % peak 4900
elseif id_country==2
    u3=0.1016; % peak 10k
elseif id_country==3
    u3=0.1377; % peak 1.1k
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% simulate the effects of the additional %%%%%
%%%%%%%%%%% control action %%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[T2, Y2]=ode45('covid_model_SIHM_distu_opt',tspan1, Y0,[]);

I_sim_u3_b=Y2(:,2);
data_sim_u3_b=eta*I_sim_u3_b;


%%%% plot the model prediction and the effects of the two control actions
type_col_fit=':r';
type_col_pr='-c';
type_col_pr_w='-.c';
type_col_a='-m';
type_col_b='-k';

idx_fig=id_country;
f_o=figure (idx_fig);
semilogy(tx(time_x+1:end),data_dQR(time_x:end),'o','color',[0.5 0.5 0.5],'markersize',8);
hold on
grid on

semilogy(Tx(idx_tx(1):idx_tpr(1)),data_sim_pr(idx_tx(1):idx_tpr(1)),type_col_fit,'Linewidth',3);
semilogy(Tx(idx_tpr(1):end),data_sim_pr(idx_tpr:end),type_col_pr,'Linewidth',3);

% semilogy(Tx(idx_tpr:idx_t3),data_sim_pr(idx_tpr:idx_t3),type_col_pr,'Linewidth',3);
% semilogy(Tx(idx_t3:end),data_sim_pr(idx_t3:end),type_col_pr_w,'Linewidth',3);


semilogy(Tx(idx_t3:end),data_sim_u3_a(idx_t3:end),type_col_a,'Linewidth',3);
semilogy(Tx(idx_t3:end),data_sim_u3_b(idx_t3:end),type_col_b,'Linewidth',3);


if id_country==1
    
    legend('DE dH data','dH fitting', 'dH preditciton ({\it{e_p}}=0.05)','\Delta u_{{opt}}=0.09, peak < 6.5k','\Delta u_{opt}=0.15, peak < 5.5k')

    ylim([2000 2e4])

elseif id_country==2
   
    ylim([1000 2e4])   

    legend('FR dH data','dH fitting','dH prediciton ({\it{e_p}}=0.01)','\Delta u_{{opt}}=0.06, peak < 11k','\Delta u_{{opt}}=0.1, peak < 10k')

elseif id_country==3

    ylim([600 4e3])

    legend('IT dH data','dH fitting','dH prediction ({\it{e_p}}=0.04)','\Delta u_{opt}=0.042, peak < 1.5k','\Delta u_{opt}=0.14, peak < 1.1k')
    
end

ylabel('dH')
T_0=datetime(2020,12,15);
if id_country==1 || id_country==2
    T_end=datetime(2021,4,25);
elseif id_country==3
    T_end=datetime(2021,4,5);
end
set(gca, 'fontsize',13)
xlim([T_0 T_end])
set(f_o,'Position',[10 10 550 200])
