
function f_cost=my_cost_fun_covid19_model_SIHM_distu_country(x_no)
% cost function for fitting procedure with the aim to find the parameter set 
% that minimize the distance between the data and  
% the simulated ones by the devised SIHM model

global  data_dQR data_QR   w m 

global   eta  n_pop H_0 beta K_p tau n_h  w0 u0
 
global w1 w2 time1  tag_change  

global time2 time3  w3  


% from normalized parameter set to real one
x=x_no.*w+m;
    
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
