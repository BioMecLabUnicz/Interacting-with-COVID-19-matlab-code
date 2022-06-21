% The script generates the results shown in Figure 2 of the ms. 
clear all
%close all

tag_change=2; % corresponding to k-1 changes of w, with w=w^{k}, so w=w^{3}

% define the vectors for rho_1 and rho_2
val_rho_1=NaN(1,4);
val_rho_2=NaN(1,4);
std_rho_1=NaN(1,4);
std_rho_2=NaN(1,4);

for id_country=1:4

     % load the parameter set for w=w^{3}
     filename_f=['res_opt/num_country_',num2str(id_country),'_data_dh_time_changes_',num2str(tag_change),'_fmincon'];                
     load(filename_f)

    w_val_vect=xfmincon([7,9]);
    % compute rho at time1 and time2 
    rho_val_vect=(1+w_val_vect)./[1 1+w_val_vect(1:end-1)];

    w_std_vect=SEopt([7,9]);
    rho_std_vect=w_std_vect./[1 1+w_val_vect(1:end-1)];

    % store rho for each country
    val_rho_1(id_country)=rho_val_vect(1);
    std_rho_1(id_country)=rho_std_vect(1);
    
    val_rho_2(id_country)=rho_val_vect(2);
    std_rho_2(id_country)=rho_std_vect(2);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% plot rho %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cell{1}={'Germany','France ','UK','Italy'};

sel_nat=[1,2,4,3];
f_w1=figure(5);
errorbar(val_rho_1(sel_nat),std_rho_1(sel_nat),'o','linewidth',2)
hold on
grid on
set(gca, 'XTick',1:4, 'XTickLabel',cell{1})
set(gca,'fontsize',16)
set(f_w1,'Position',[10 10 300 150])
ylim([1.3 2.5])

str_title_w1='{\rho}_{1}';
ylabel(str_title_w1);

f_w2=figure(6);
errorbar(val_rho_2(sel_nat),std_rho_2(sel_nat),'o','linewidth',2)
hold on
grid on
set(gca, 'XTick',1:4, 'XTickLabel',cell{1})
set(gca,'fontsize',16) 
set(f_w2,'Position',[10 10 300 150])
str_title_w2='{\rho}_{2}';
ylabel(str_title_w2);

figure(6)
ylim([1.3 2.5])
