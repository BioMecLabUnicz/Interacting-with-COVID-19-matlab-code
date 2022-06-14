
function f_cost=my_cost_fun_covid19_model_SIHM_distu_country_plot(x,idx_f,plot_data,type_col)

% function used for plotting the behavior of SIHM model for a given
% parameter set , where x is the parameter vector, idx_f is the figure number,
% plot data is a tag that is true for plotting the data, and type_col is a
% string for defining the line and color of the simulated SIHM_0 model
% for example,type_col='-k';

global  data_dQR data_QR

global   eta  n_pop H_0 beta K_p tau  n_h

global  tag_change

global time1 time2 time3  time4 tag_kp_tau

global w1 w2 w3 w4 u0 w0


global tag_country 


if tag_change==1
    
    
    eta_I0=x(1);
    I0=eta_I0/eta;
    beta=x(2);
    K_p=x(3);
    n_h=x(4);
    tau=x(5);
    time1=x(6);
    w1=x(7);
    
    
    
elseif tag_change==2
    
    eta_I0=x(1);
    I0=eta_I0/eta;
    beta=x(2);
    K_p=x(3);
    n_h=x(4);
    tau=x(5);
    time1=x(6);
    w1=x(7);
    time2=x(8);
    w2=x(9);
    
elseif tag_change==3
    
    
    eta_I0=x(1);
    I0=eta_I0/eta;
    beta=x(2);
    K_p=x(3);
    n_h=x(4);
    tau=x(5);
    time1=x(6);
    w1=x(7);
    time2=x(8);
    w2=x(9);
    time3=x(10);
    w3=x(11);
    
elseif tag_change==4
    
    
    eta_I0=x(1);
    I0=eta_I0/eta;
    beta=x(2);
    K_p=x(3);
    n_h=x(4);
    tau=x(5);
    time1=x(6);
    w1=x(7);
    time2=x(8);
    w2=x(9);
    time3=x(10);
    w3=x(11);
    
    time4=x(12);
    w4=x(13);
    
    
    
    
else
    eta_I0=x(1);
    I0=eta_I0/eta;
    beta=x(2);
    K_p=x(3);
    n_h=x(4);
    tau=x(5);
    
end

    
dQR_data=data_dQR; % in our case corresponds to dH
QR_data=data_QR;

QR_0=H_0;
u0=0;
w0=0;

size_time_window=length(dQR_data);
v_days=0:size_time_window;



tspan=v_days;

Y0=zeros(4,1);
Y0(1,1)=n_pop-I0-QR_0;
Y0(2,1)=I0;
Y0(3,1)=QR_0;


[T, Y]=ode45('covid_model_SIHM_distu',tspan, Y0,[]);

QR_sim=Y(:,3);
I_sim=Y(:,2);
data_sim=eta*I_sim;

if tag_country==2 || tag_country==1 
    tspan1=0:1/7:size_time_window;
elseif tag_country==3 || tag_country==4 
     tspan1=0:1/4:size_time_window;
end


[T1, Y1]=ode45('covid_model_SIHM_distu',tspan1, Y0,[]);

I_sim1=Y1(:,2);

data_sim1=eta*I_sim1;


if  sum(data_sim<=0)==0
 
    if length(QR_sim)==length(QR_data)
        
        f_cost_dQR_data=0;


        for idx=1:length(dQR_data)
            
            if dQR_data(idx)<=0 || isnan(dQR_data(idx))
                
                err_sample=0;
            else
                err_sample=(log(data_sim(idx+1))-log(dQR_data(idx))).^2;
            
            end
                
         
                
                
                
            
            f_cost_dQR_data=f_cost_dQR_data+err_sample;
           
        end
        
        
        
            
           
        
    else
        
        f_cost_dQR_data=Inf;
       
    end
    
  
    
else
    
    f_cost_dQR_data=Inf;
    
end

f_cost=f_cost_dQR_data;
figure (idx_f)
subplot(1,1,1)

%set(gca, 'YScale', 'log')

% semilogy(v_days,QR_data,'om')
%semilogy(T(2:end),dQR_data,'ok');
%t=datetime(2020,3,19) + caldays(1:length(v_days));

if tag_country==3
    t=datetime(2020,2,19) + caldays(0:length(v_days)-1);
    dH_sim=data_sim1;
    t1=t(1):1/4:t(end);
elseif tag_country==2 
    %t=datetime(2020,3,8) + calweeks(0:48);
    t=datetime(2020,3,8) + calweeks(0:size_time_window);
    dH_sim=data_sim1;
    t1=t(1):t(end);
elseif tag_country==1 
    %t=datetime(2020,3,8) + calweeks(0:48);
    t=datetime(2020,3,1) + calweeks(0:size_time_window);
    dH_sim=data_sim1;
    t1=t(1):t(end);

elseif tag_country==4 
    t=datetime(2020,3,22) + caldays(0:length(v_days)-1);
    dH_sim=data_sim1;
    t1=t(1):1/4:t(end);
    

    
end


if plot_data
    semilogy(t(2:end),data_dQR,'o','color',[0.5 0.5 0.5],'markersize',6);
  
    hold on
    grid on
end

if tag_change==2
    type_col1='-r';
    type_col2='-.r';
    type_col3=':r';
%     type_col2='-r';
%     type_col3='-g';
%     type_col4='-c';
%     type_col5='-m';
elseif tag_change==1
    type_col1='-b';
    type_col2='-.b';
%      type_col1='-c';
%     type_col2='-m';

elseif tag_change>=3
    type_col1='-g';
    type_col2='-.g';
    type_col3=':g';
    type_col4='--g';
    type_col5='-m';
end
Tx=t1;

if tag_change>=1
    if tag_country==2 || tag_country==1 
         idx_t1=find(T1>time1-1/7);

         if tag_change==2
             idx_t2=find(T1>time2-1/7);
         end

         if tag_change==3
             idx_t2=find(T1>time2-1/7);
             idx_t3=find(T1>time3-1/7);
         end
         if tag_change==4
             idx_t2=find(T1>time2-1/7);
             idx_t3=find(T1>time3-1/7);
             idx_t4=find(T1>time4-1/7);
         end
    elseif tag_country==3 || tag_country==4
        idx_t1=find(T1>time1-1/4);
        if tag_change==2
            idx_t2=find(T1>time2-1/4);
        end

         if tag_change==3
             idx_t2=find(T1>time2-1/4);
             idx_t3=find(T1>time3-1/4);
         end
         
    end
else
   idx_t1=[]; 
end
if isempty(idx_t1)
    idx_t1=length(Tx);
    
    plot_ph_1=1;
else
    plot_ph_1=0;
end

plot_ph_3=0;
plot_ph_2=0;
plot_ph_4=0;

if tag_change==2 
    
    if isempty(idx_t2)
        idx_t2=length(T1);
        plot_ph_2=1;
    end
elseif tag_change==3 
    
    if isempty(idx_t2)
        idx_t2=length(T1);
        plot_ph_2=1;
    end
    if isempty(idx_t3)
        idx_t3=length(T1);
        plot_ph_3=1;
    end
    
    
    
elseif tag_change==4
    
    if isempty(idx_t2)
        idx_t2=length(T1);
        plot_ph_2=1;
    end
    if isempty(idx_t3)
        idx_t3=length(T1);
        plot_ph_3=1;
    end
    if isempty(idx_t4)
        idx_t4=length(T1);
        plot_ph_4=1;
    end
else
    idx_t3=length(T1);
    idx_t2=length(T1);
    idx_t4=length(T1);
end




if plot_ph_1
    semilogy(Tx,dH_sim,type_col,'Linewidth',3);
    

elseif  plot_ph_2
        semilogy(Tx(idx_t1:end),dH_sim(idx_t1:end),type_col,'Linewidth',3);
        
elseif  plot_ph_3
        semilogy(Tx(idx_t2:end),dH_sim(idx_t2:end),type_col,'Linewidth',3);
elseif  plot_ph_4
        semilogy(Tx(idx_t3:end),dH_sim(idx_t3:end),type_col,'Linewidth',3);
        

elseif tag_change==2
         
    semilogy(Tx(1:idx_t1),dH_sim(1:idx_t1),type_col1,...
        Tx(idx_t1:idx_t2),dH_sim(idx_t1:idx_t2),type_col2,...
         Tx(idx_t2:end),dH_sim(idx_t2:end),type_col3,'Linewidth',3);
   
    
elseif tag_change==3
         
    semilogy(Tx(1:idx_t1),dH_sim(1:idx_t1),type_col1,...
        Tx(idx_t1:idx_t2),dH_sim(idx_t1:idx_t2),type_col2,...
        Tx(idx_t2:idx_t3),dH_sim(idx_t2:idx_t3),type_col3,...
         Tx(idx_t3:end),dH_sim(idx_t3:end),type_col4,'Linewidth',3);
     

elseif tag_change==4
         
    semilogy(Tx(1:idx_t1),dH_sim(1:idx_t1),type_col1,...
        Tx(idx_t1:idx_t2),dH_sim(idx_t1:idx_t2),type_col2,...
        Tx(idx_t2:idx_t3),dH_sim(idx_t2:idx_t3),type_col3,...
         Tx(idx_t3:idx_t4),dH_sim(idx_t3:idx_t4),type_col4,...
          Tx(idx_t4:end),dH_sim(idx_t4:end),type_col5,'Linewidth',3);
     
     
elseif tag_change==1 || tag_kp_tau==1
    
    semilogy(Tx(1:idx_t1),dH_sim(1:idx_t1),type_col1,...
        Tx(idx_t1:end),dH_sim(idx_t1:end),type_col2,...
       'Linewidth',3);
end
   

hold on
grid on

if plot_data
    xtickformat('MMM')
    xlim([t(1) t(end)])
end



