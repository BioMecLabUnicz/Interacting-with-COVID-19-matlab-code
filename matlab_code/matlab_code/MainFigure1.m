% The script generates the results shown in Figure 1 of the ms. 
% In the code, there are flags (do_opt and do_opt_fmincon), 
% which are true if you choose to run GA optimization (do_opt) and fmicnon (do_opt_fmincon),
% otherwise you can load the files with the results of the hybrid optimization (GA+fmincon) 
% within the folder res_opt.
% In the script you set the parameter id_country in order to show the
% results for a given country:
% id_country=1, plot the results for Germany, 2 for France, 3 for Italy, 4 for UK.

clear all
close all
clc
global  data_dQR data_QR  n_pop  eta alpha H_0 time1

global beta K_p tau tag_change  

global m w n_h 

global time2  w1  w2 w3 time3 u1 u2 u3

global tag_country


id_country=1; % id for country 1=Germany, 2=France, 3=Italy, 4=UK

do_opt=0; %  if true do ga optimization otherwise load the saved results

do_opt_fmincon=0; % if true then do fmincon optimization, otherwise load the saved results
save_fmincon=1; % if true then save the fmincon optimization results
use_fmincon_res=1; %  if true then load the fmincon optimization results

n_changes=3; % number of stepwise changes tested for w, i.e. w=w^{\{4\}}

perf_FVAL_ind=NaN(1,n_changes+1);
perf_aic_ind=NaN(1,n_changes+1);
perf_bic_ind=NaN(1,n_changes+1);
perf_fpe_ind=NaN(1,n_changes+1); 

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


if id_country==1    
    data_dH=num_dH_Ger;
    H_0=100;
    data_H= cumsum([H_0; data_dH]);
    data_dQR=data_dH(1:46); % data_dQR represents the dH data for fitting
    data_QR=data_H(1:47);    % data_QR represents the cumulative sum of the hospital admissions dH
elseif id_country==2
    
    data_dH=data_dH_week_11_new;
    
    data_H= cumsum([H_0; data_dH]); 
    data_dQR=data_dH(1:45);  % data_dQR represents the dH data for fitting
    data_QR=data_H(1:46);     % data_QR represents the cumulative sum of the hospital admissions dH
    
elseif id_country==3
    
    
    H_0=0;
    data_cumsumH3=[H_0; cumsum(data_dH3)];
    
    data_dQR=data_dH3(1:333); % data_dQR represents the dH data for fitting
    data_QR=data_cumsumH3(1:334); % data_QR represents the cumulative sum of the hospital admissions dH
    
elseif id_country==4
    
    H_0=cum_sum_H0;
    data_dQR=data_dH2(1:301); % data_dQR represents the dH data for fitting
    data_QR=data_cumsum_H2(1:302);  % data_QR represents the cumulative sum of the hospital admissions dH
    
end


size_time_window=length(data_QR);

v_days=0:size_time_window-1;

%     figure
% 
%     semilogy(v_days(2:end),data_dQR,'or')
%     hold on
%     grid on
%     semilogy(v_days,data_QR,'sb')

if id_country==3 ||  id_country==4
    
    eta=0.006;
    alpha=0.094;
elseif id_country==1 ||  id_country==2
    
    eta=0.006*7;
    alpha=0.094*7;
end



for tag_change=0:n_changes 
    % tag_change=k-1, corresponding to the stepwise changes of w, with w=w^{\{k\}},
    
    % tag_change=0, i.e. w=w^{1}=w0=0 (w0 in the code corresponds to w_1 in the manuscript)
    % in general wi in the code corresponds to w_{i+1}  in the ms
    % tag_change=1, i.e. w=w^{2}, 1 stepwise change of w from w0 to w1 at Tc1 (time1) (w1 corresponds to w_2 in the ms)
    
    % tag_change=2; i.e w=w^{3},  2 stepwise changes, 1) from w0 to w1 at Tc1 (time1) 
                                                                            %  2) from w1 to w2 at Tc2 (time2) (w2 corresponds to w_3 in the ms)
    % tag_change=3; i.e w=w^{4}, 3 stepwise changes,  1) from w0 to w1 at Tc1 (time1) 
                                                                             % 2) from w1 to w2 at Tc2 (time2) 
                                                                             % 3) from w2 to w3 at Tc3 (time3) (w3 corresponds to w_4 in the ms)
                                                                             
    if do_opt  % do the ga optimization
        
        % number of parameters to be identified
        if tag_change==0
            
            n_id_par=5;
            
        elseif tag_change==1
            
            n_id_par=7;
            
        elseif tag_change==2
            
            n_id_par=9;
            
        elseif tag_change==3
            
            n_id_par=11;
            
        end
        init_num_it=1;
        end_num_it=1; % number of ga repetitions
        
        num_it_h_pop=3;
        num_it_init_pop_id=2;
        mtx_par=NaN(end_num_it,n_id_par+1);
        mtx_par_ga=NaN(end_num_it,n_id_par+1);
        
        for num_it=init_num_it:end_num_it
            
            % name of the file for storing the optimization results
            filename=['res_opt/num_country_',num2str(id_country),'_data_dh_time_changes_',num2str(tag_change),'_num_it_new_',num2str(num_it)];
            
            % initialize the procedure
            if tag_change==1
                if id_country==1
                    
                    Init_pop=[ ];
                    
                elseif id_country==3
                    
                    Init_pop=[];     
                
                elseif  id_country==4
                
                    Init_pop=[];
                    
                elseif id_country==2
                    
                    Init_pop=[];
                    
                end
                
            elseif tag_change==2
                
                if id_country==1
                    Init_pop=[ ];
                elseif id_country==4
                    Init_pop=[ ];                  
                elseif id_country==3
                    Init_pop=[];    
                elseif id_country==2
                    Init_pop=[ ];
                end
     
            elseif tag_change==3
                
                
                 if id_country==1
                     Init_pop=[];
                 elseif id_country==2                
                     Init_pop=[];
                 elseif id_country==3
                     Init_pop=[];
                elseif id_country==4
                    Init_pop=[];
                 end
            elseif tag_change==0
                
                
                if id_country==1
                    
                    Init_pop=[ ];
                elseif id_country==3
                    Init_pop=[  ];
                elseif id_country==4
                    Init_pop=[];
                elseif id_country==2
                    Init_pop=[];
                end
                
            end
            
            % define the bounds
            p1_l=0;
            p1_h=5e4;
            
            if tag_change==1
                
                if id_country==1 ||id_country==2
                    
                    LB_pw=[p1_l 0   0   0   0  25  0 ];
                    
                    UB_pw=[p1_h 7  1e5 10  15  35 10];
                    
                elseif id_country==3
                    LB_pw=[p1_l 0   0    0  0  210  0  ];
                    
                    
                    UB_pw=[p1_h 1  1e5  10 120 250 10 ];
                    
                elseif id_country==4
                    
                    
                    LB_pw=[p1_l 0   0    0  0   160  0];
                    
                    UB_pw=[p1_h 1  1e6  10 120 210  10];
                end
                
            elseif tag_change==2
                
                if id_country==1
                    LB_pw=[p1_l 0   0   0   0  25  0 35 0];
                    
                    UB_pw=[p1_h 7  1e5 10 15  35 10 42 10];
                    
                elseif id_country==4
                    
                    LB_pw=[p1_l 0   0    0  0  166  0  240 0];
                    
                    UB_pw=[p1_h 1  1e4  10 120 197 10 280 10];
                    
                elseif id_country==2
                    
                    LB_pw=[p1_l 0   0   0   0  25  0 36 0];
                    
                    
                    UB_pw=[p1_h 7  1e5 10  15  35 10 43 10];
                elseif id_country==3
                    
                    LB_pw=[p1_l  0   0    0  0  230  0 290 0];
                    
                    UB_pw=[p1_h  1  1e4  10 120 245 10 305 10];
                    
                end
            elseif tag_change==3
                
                if id_country==1
                     
                    LB_pw=[p1_l  0   0    0   0   6    0   23  0 35 0];
                    
                    UB_pw=[p1_h 7  1e5 10 15  22 10  34 10 42 10];
                    
                elseif id_country==2
                     
                    LB_pw=[p1_l  0   0    0   0   6    0   22  0 36 0];
                    
                    UB_pw=[p1_h 7  1e5 10 15  21 10  35 10 43 10];   
                 
                elseif id_country==3 
                    
                   LB_pw=[p1_l    0    0     0    0  50    0  180  0 270   0];
                    
                    UB_pw=[p1_h  1  1e4  10 120 170 10 240 10 310 10];
              
                  
                elseif id_country==4
                    LB_pw=[p1_l 0     0    0    0   40  0  165  0 230 0];
                    
                    UB_pw=[p1_h 1  1e6  10 120 120 10 210 10 280 10];
                end
                
            elseif tag_change==0
                
                if id_country==3 ||  id_country==4
                    LB_pw=[p1_l   0   0     0  0    ];
                    
                    UB_pw=[p1_h  1  1e4 10 120 ];
                    
                    
                elseif id_country==1 ||  id_country==2
                    
                    LB_pw=[p1_l 0   0   0   0   ];
                    
                    UB_pw=[p1_h 7  1e5 10  15   ];
                end
                
            end
            
            
            m=(LB_pw+UB_pw)/2;
            w=(UB_pw-LB_pw)/2;
            LB_pw_no=(LB_pw-m)./w;
            UB_pw_no=(UB_pw-m)./w;
            
            if isempty(Init_pop)
                Init_pop_no=[];
            else
                
                Init_pop_no=(Init_pop-m)./w;
            end
            if num_it>num_it_init_pop_id
                Init_pop=mtx_par_ga(1:num_it-1,1:n_id_par);
                
                Init_pop_nor=[Init_pop_no;Init_pop];
            else
                Init_pop_nor=Init_pop_no;
            end
            
            if num_it>num_it_h_pop
                
                myoptions=gaoptimset('InitialPopulation',Init_pop_nor,'PopulationSize',500);
            else
                myoptions=gaoptimset('InitialPopulation',Init_pop_nor);
            end
            
            
            [xoptGA_pw_no,JobjGA_pw,exitflag_pw,output_pw,population_pw,scores_pw]=ga(@my_cost_fun_covid19_model_SIHM_distu_country,length(LB_pw_no),...
                [],[],[],[],LB_pw_no,UB_pw_no,[],myoptions);
            
            xoptGA_pw=xoptGA_pw_no.*w+m;
            
            xopt=xoptGA_pw;
            % Store the identified parameters
            
            mtx_par(num_it,:)=[xoptGA_pw JobjGA_pw];
            mtx_par_ga(num_it,:)=[xoptGA_pw_no JobjGA_pw];
            save(filename,'mtx_par');
            
            eval_single_opt=1;
            if eval_single_opt
                % simulate the model with the optimized parameters            
                I0=xopt(1)/eta;
                beta=xopt(2);                
                K_p=xopt(3);
                n_h=xopt(4);
                tau=xopt(5);
                if tag_change>=1
                    
                    time1=xopt(6);
                    w1=xopt(7);
                end
                if tag_change>=2
                    time2=xopt(8);
                    w2=xopt(9);
                    
                end
                if tag_change==3
                    
                    time3=xopt(10);
                    w3=xopt(11);
                    
                end
                
                
                tspan=v_days;
                
                Y0=zeros(4,1);
                
                Y0(1,1)=n_pop-I0-H_0;
                Y0(2,1)=I0;
                Y0(3,1)=H_0;
                
                [T, Y]=ode45('covid_model_SIHM_distu',tspan, Y0,[]);
                
                S_pop=Y(:,1);
                I_pop=Y(:,2);
                cumH_pop=Y(:,3);
                M_sim=Y(:,4);
                
                
                f_idx=figure (10+id_country);
                semilogy(T,S_pop,'--b',T,I_pop,'--r',T,cumH_pop,'--m',T,M_sim,'-k','Linewidth',1);
                hold on
                grid on
                
                % semilogy(v_days,QR_data,'om')
                semilogy(v_days(2:end),data_dQR','sr')
                semilogy(v_days,eta*I_pop,'--c')
                legend('S','I','H','M','dH data','dH sim')
                
            end
            
        end
        
    else
        
        % load the stored identifications (six) by ga        
        num_it=6;
        filename=['res_opt/num_country_',num2str(id_country),'_data_dh_time_changes_',num2str(tag_change),'_num_it_',num2str(num_it)];
        
        load (filename)
        % select the best fitting solution
        
    end   
    if do_opt_fmincon % perform fmincon optimization
        
        % select the best fitting solution from GA as initial condition for
        % fmincon
        [~,num_id_sel]=min(mtx_par(:,end));
        
        % optimezed parameters
        x=mtx_par(num_id_sel,1:end-1);
        
        % define the bounds for fmincon
        bound=x*.4;
        
        % specific bounds as for the time changes (day or weeks)
        if id_country==1
            
            if tag_change==0
                
                bound_x=x*.2;
                num_idx_x=3;
                
                bound_y=[];
                num_idx_y=[];
                
                bound_z=[];
                num_idx_z=[];
                
            elseif tag_change==1
                
                bound_x=x*.25;
                num_idx_x= 6;
                
                bound_y=[];
                num_idx_y=[];
                
                bound_z=[];
                num_idx_z=[];
                
            elseif tag_change==2
                
                bound_x=x*.1;
                num_idx_x=[1 3];
                
                bound_y=x*.06;
                num_idx_y=[6 8];
                
                bound_z=x*.04;
                num_idx_z=[9];
            else
                bound_x=[];
                num_idx_x=[];
                num_idx_y=[];
                bound_y=[];
                bound_z=[];
                num_idx_z=[];
            end
            
        elseif id_country==2
            
            if tag_change==2
                
                bound_z=x*.3;
                num_idx_z=[2 4];
                
                bound_y=x*.08;
                num_idx_y=[ 7 9];
                
                bound_x=x*.05;
                num_idx_x=[ 6 8];
                
                
            elseif tag_change==1
                
                bound_x=x*.15;
                num_idx_x= 6;
                
                bound_y=x.*425;
                num_idx_y=2;
                
                bound_z=[];
                num_idx_z=[];
                
            else
                
                num_idx_x=[];
                num_idx_y=[];
                bound_x=[];
                bound_y=[];
                bound_z=[];
                num_idx_z=[];
            end
            
        elseif id_country==3
            
            bound_z=[];
            num_idx_z=[];
            
            if tag_change==2
                
                bound_x=x*.025;
                num_idx_x=[  6 8];
                
                bound_y=x*.08;
                num_idx_y=[7 9];
                
            elseif tag_change==1
                
                bound_x=x*.1;
                num_idx_x= [1 2  6];
                
                bound_y=x*.3;
                num_idx_y= 3;
                
                
            else
                
                num_idx_x=[];
                num_idx_y=[];
                bound_x=[];
                bound_y=[];
                
                
            end
            
        elseif id_country==4
            
            if tag_change==2
                
                bound_x=x*.25;
                num_idx_x=1:5;
                
                bound_y=x*.06;
                num_idx_y=[ 8 9];
                
                bound_z=x*.1;
                num_idx_z=[ 6 7];
                
            elseif tag_change==1
                
                bound_x=x*.2;
                num_idx_x= 1:5;
                
                bound_y=x*.06;
                num_idx_y=[6 7];
                
                bound_z=x;
                num_idx_z=[];
            else
                bound_x=x*.2;
                num_idx_x=[1 4];
                
                num_idx_y=[];
                bound_y=[];
                
                bound_z=x;
                num_idx_z=[];
            end
            
        end
        
        LB_pw=x-bound;
        UB_pw=x+bound;
        
        LB_pw(num_idx_x)=x(num_idx_x)-bound_x(num_idx_x);
        UB_pw(num_idx_x)=x(num_idx_x)+bound_x(num_idx_x);
        LB_pw(num_idx_y)=x(num_idx_y)-bound_y(num_idx_y);
        UB_pw(num_idx_y)=x(num_idx_y)+bound_y(num_idx_y);
        LB_pw(num_idx_z)=x(num_idx_z)-bound_z(num_idx_z);
        UB_pw(num_idx_z)=x(num_idx_z)+bound_z(num_idx_z);
        
        %
        
        m=(LB_pw+UB_pw)/2;
        
        w=(UB_pw-LB_pw)/2;
        
        
        LB_pw_no=(LB_pw-m)./w;
        
        UB_pw_no=(UB_pw-m)./w;
        
        xoptGA_pw_no=(x-m)./w;
        
        [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(@my_cost_fun_covid19_model_SIHM_distu_country,xoptGA_pw_no,[],[],[],[],LB_pw_no(1:end),UB_pw_no(1:end));
        SE=sqrt(diag(inv(HESSIAN)));
        SEopt=SE'.*w
        
        xfmincon=X.*w+m
        xopt=xfmincon;
        n_s=length(data_dQR);
        n_p=length(xfmincon);
        
        aic_ind=n_s*log(FVAL/n_s) +2*n_p;

        bic_ind=n_s*log(FVAL/n_s) +n_p*log(n_s);
                
        fpe_ind=(FVAL/n_s)*((1+n_p/n_s)/(1-n_p/n_s));
        
        perf_FVAL_ind(tag_change+1)=FVAL/n_s;
        perf_aic_ind(tag_change+1)=aic_ind;
        perf_bic_ind(tag_change+1)=bic_ind;
         
        perf_fpe_ind(tag_change+1)=fpe_ind;
        
        
        if save_fmincon
            filename_f=['res_opt/num_country_',num2str(id_country),'_data_dh_time_changes_',num2str(tag_change),'_new_fmincon'];
            
            save(filename_f,'xfmincon','FVAL', 'SEopt','aic_ind','bic_ind','fpe_ind');
        end
    else
        
        
        if use_fmincon_res% load the fmincon results
            
            filename_f=['res_opt/num_country_',num2str(id_country),'_data_dh_time_changes_',num2str(tag_change),'_fmincon'];
            
            load(filename_f)
            
            xopt=xfmincon;
        else
            
            num_it=6;
            
            filename=['res_opt/num_country_',num2str(id_country),'_data_dh_time_changes_',num2str(tag_change),'_num_it_',num2str(num_it)];
            
            load (filename)
            [FVAL,num_id_sel]=min(mtx_par(:,end));           
            
            % optimezed parameters
            xopt=mtx_par(num_id_sel,1:end-1);
            
        end
    end
    
    % compute scores
    
    n_s=length(data_dQR);
    n_p=length(xopt);
    
    aic_ind=n_s*log(FVAL/n_s) +2*n_p;

    bic_ind=n_s*log(FVAL/n_s) +n_p*log(n_s);
    
    fpe_ind=(FVAL/n_s)*((1+n_p/n_s)/(1-n_p/n_s));
  
    perf_FVAL_ind(tag_change+1)=FVAL/n_s;
    perf_aic_ind(tag_change+1)=aic_ind;
    perf_bic_ind(tag_change+1)=bic_ind;
    perf_fpe_ind(tag_change+1)=fpe_ind;
    
    
    thr_score=0.05;
    if tag_change>0
        
        % compare the scores
        comp_fval=perf_FVAL_ind(tag_change+1)-(perf_FVAL_ind(tag_change)-abs(perf_FVAL_ind(tag_change)*thr_score));
        comp_aic=perf_aic_ind(tag_change+1)-(perf_aic_ind(tag_change)-abs(perf_aic_ind(tag_change)*thr_score));
        comp_bic=perf_bic_ind(tag_change+1)-(perf_bic_ind(tag_change)-abs(perf_bic_ind(tag_change)*thr_score));
        comp_fpe=perf_fpe_ind(tag_change+1)-(perf_fpe_ind(tag_change)-abs(perf_fpe_ind(tag_change)*thr_score));
        
         if (comp_fval<0)&&(comp_aic<0)&& (comp_bic<0)&& (comp_fpe<0)
            all_scores_improved=1;
            
            
            if tag_change==n_changes
                disp('increase the number of changes')
            end
        else
            all_scores_improved=0;
            %%% found the model with index equal to tag_change-1
            id_sel=tag_change-1;
            if id_sel>1
                text=['best model with ',num2str(id_sel),' changes'];
            else
                text=['best model with ',num2str(id_sel),' change'];
            end
            disp(text);
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%% Plot fitting results %%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            tag_country=id_country;
            idx_fig=tag_country;
            type_col='-k';
                    
            for  tag_change=0:id_sel % plot the results for SIHM with w^{k}, k=tag_change+1 
               
                filename_f=['res_opt/num_country_',num2str(id_country),'_data_dh_time_changes_',num2str(tag_change),'_fmincon'];        
                load(filename_f)
                if tag_change==0
                    plot_data=1;
                else
                    plot_data=0;
                end
                my_cost_fun_covid19_model_SIHM_distu_country_plot(xfmincon,idx_fig,plot_data, type_col)
                
            end
            f=figure(idx_fig);

            set(gca, 'fontsize',13)

            set(f,'Position',[10 10 350 200]);
            ylabel('dH')


            if id_country==1
                legend({'DE dH','w^{\{1\}} (w=w_1)' , 'w^{\{2\}} (w=w_1)','w^{\{2\}} (w=w_2)', 'w^{\{3\}} (w=w_1)','w^{\{3\}} (w=w_2)','w^{\{3\}} (w=w_3)'},'fontsize',8);
                ylim([10 3e4])
            elseif  id_country==2
                ylim([10 4e4])
                legend({'FR dH','w^{\{1\}} (w=w_1)' , 'w^{\{2\}} (w=w_1)','w^{\{2\}} (w=w_2)', 'w^{\{3\}} (w=w_1)','w^{\{3\}} (w=w_2)','w^{\{3\}} (w=w_3)'},'fontsize',8);

            elseif id_country==3
                
                legend({'IT dH','w^{\{1\}} (w=w_1)' , 'w^{\{2\}} (w=w_1)','w^{\{2\}} (w=w_2)', 'w^{\{3\}} (w=w_1)','w^{\{3\}} (w=w_2)','w^{\{3\}} (w=w_3)'},'fontsize',8);
                    ylim([10 1e4])
            elseif id_country==4
                ylim([60 6e3])
                legend({'UK dH','w^{\{1\}} (w=w_1)' , 'w^{\{2\}} (w=w_1)','w^{\{2\}} (w=w_2)', 'w^{\{3\}} (w=w_1)','w^{\{3\}} (w=w_2)','w^{\{3\}} (w=w_3)'},'fontsize',8);
            end
            
            return 
            
         end         
         
    end
   
       
end

    

