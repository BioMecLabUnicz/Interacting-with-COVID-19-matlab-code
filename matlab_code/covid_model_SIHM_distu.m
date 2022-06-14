function xdot=covid_model_SIHM_distu(t,x,flag)


global alpha beta  eta n_pop tau tag_change  time3 w4  time4  K_p w0 u0

global w1 w2 w3 time1  time2 n_h 
s=x(1);

I=x(2);

q_plus_r=x(3);

m=x(4);

tau_val=tau;

if tag_change==1
    if t>=time1
        w_val=w1;
    else
        w_val=w0;
    end
        
elseif tag_change==2
    
    if t>=time1 && t<time2
        w_val=w1;
        
    elseif t>=time2
        w_val=w2;
                
    else
        w_val=w0;
        
    end
    
elseif tag_change==3
    
    
    if t>=time1 && t<time2
        w_val=w1;
       
    elseif t>=time2 && t<time3
        w_val=w2;
    elseif t>=time3
        w_val=w3;       
    else
        w_val=w0;
    end
    
elseif tag_change==4
    
    if t>=time1 && t<time2
        w_val=w1;
    
    elseif t>=time2 && t<time3
        w_val=w2;
        
    elseif t>=time3 && t<time4
        w_val=w3;
              
    elseif t>=time4
        w_val=w4;
               
    else
        w_val=w0;
             
    end
else
    w_val=w0;
    
end

beta_val=beta;

u_val=u0;


beta_f=beta_val/(1+(m/K_p)^n_h)*(1+w_val)*(1-u_val);
s_dot= -beta_f*s*I/n_pop;

i_dot= beta_f*s*I/n_pop - (alpha+eta)*I;

q_plus_r_dot= eta*I ;

m_dot=(eta*I - m)/tau_val;

xdot=[s_dot;i_dot;q_plus_r_dot;m_dot ];

