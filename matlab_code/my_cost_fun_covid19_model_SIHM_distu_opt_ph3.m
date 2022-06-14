
function f_cost=my_cost_fun_covid19_model_SIHM_distu_opt_ph3(x_no)
%  cost function for implementing a control strategy with the aim to keep the
%  peak value below the desired threshold

global data_dQR data_QR   w m 

global eta  n_pop H_0  u0 I0
 
global tag_change  

global time3  data_week data_peak

global u1 u2 u3 u4 u5 ka_1 ka_2



% from normalized parameter set to real one
x=x_no.*w+m;

%%% for our case tag_change==3;% u0=0; % w0=0;


if tag_change==1
    
        u1=u0;
       

elseif tag_change==2
    
        u1=u0;
        u2=u0;
        
            
elseif tag_change==3
    
    
        u1=u0;
        u2=u0;
        u3=x(1);
        
elseif tag_change==4
    
    
        u1=x(1);
        u2=x(2);
        u3=x(3);
        u4=x(4);
       
elseif tag_change==5
    
    
        u1=x(1);
        u2=x(2);
        u3=x(3);
        u4=x(4);
        u5=x(5);
end
   

dQR_data=data_dQR; % in our case corresponds to dH
QR_data=data_QR;


QR_0=H_0;

size_time_window=length(dQR_data);

if data_week==1
    v_days=0:1/4:size_time_window;
elseif data_week==7
    v_days=0:1/14:size_time_window;
end
    


tspan=v_days;

Y0=zeros(4,1);
Y0(1,1)=n_pop-I0-QR_0;
Y0(2,1)=I0;
Y0(3,1)=QR_0;

[T, Y]=ode45('covid_model_SIHM_distu_opt',tspan, Y0,[]);


QR_sim=Y(:,3);
I_sim=Y(:,2);

idx_ts=find(T>=time3);
data_sim=eta*I_sim(idx_ts(1):end);
u_v=zeros(length(T)-idx_ts(1)+1,1);
for idx_t=idx_ts:length(T)

        if T(idx_t)>=time3            
            u_v(idx_t)=u3;
        else
            u_v(idx_t)=u0;
        end

end

u_v=u_v(idx_ts(1):end);

%%%% Italy %%%%%%
% data_week=1;
% data_peak=1600;
% ka_2=(1/(data_week*(data_peak-100)-data_week*data_peak)^2)^-1;
% x_k=ka_2/900^2; ka_1=(1-x_k)/(0.25^2);

%%%% Germany %%%%
% data_week=7;
% data_peak=1000;
% ka_2=(1/(data_week*(data_peak-100)-data_week*data_peak)^2)^-1;
% x_k=ka_2/(data_week*500)^2; ka_1=(1-x_k)/(0.25^2);

%%%% France %%%%
% data_week=7;
% data_peak=1500;
% ka_2=(1/(data_week*(data_peak-100)-data_week*data_peak)^2)^-1;
% x_k=ka_2/(data_week*500)^2; ka_1=(1-x_k)/(0.25^2);

n_sim=length(data_sim);
data_sim_th=zeros(n_sim,1);
for idx_d=1:length(data_sim)
    data_sim_th(idx_d)=min([data_sim(idx_d)-data_week*data_peak,0]);
    
end
f_cost=sum(ka_1*(u_v).^2+ka_2./(data_sim_th).^2); 
