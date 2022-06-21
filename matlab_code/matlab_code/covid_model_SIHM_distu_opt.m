function xdot=covid_model_SIHM_distu_opt(t,x,flag)

% function used for defining the ode system implementing 
% the SIHM model with beta_f=u*(1+w)*(1-Delta_u)

global alpha beta eta n_pop tag_change time3 w4  time4 time5 w5 K_p w0 u0

global w1 w2 w3 time1 time2  n_h  tau

global u1 u2 u3 u4 u5

s=x(1);

I=x(2);

q_plus_r=x(3);

m=x(4);

% tag_change=k-1, corresponding to the stepwise changes of w, with w=w^{k},
% in the tested case tag_change=3, i.e. w=w^{4} with 4-1=3 stepwise changes, 
% note that wi in the code corresponds to w_{i+1}  in the ms
% u1=u2=0 and u3=Delta_u is the
% additional control action to be implemented in response to w3 with the
% aim to keep the peak value below the desired threshold

if tag_change==1
    if t>=time1
        w_val=w1;
        u_val=u1;
    
    else
        w_val=w0;
        u_val=u0;
     end
    
    
    
elseif tag_change==2
    
    if t>=time1 && t<time2
        w_val=w1;
        u_val=u1;
        
        
        
    elseif t>=time2
        w_val=w2;
        u_val=u2;
        
        
    else
        w_val=w0;
        u_val=u0;
        
    end
    
    
elseif tag_change==3
    
    
    if t>=time1 && t<time2
        w_val=w1;
        u_val=u1;
       
    elseif t>=time2 && t<time3
        w_val=w2;
        u_val=u2;
       
    elseif t>=time3
        w_val=w3;
        u_val=u3;
        
    else
        w_val=w0;
        u_val=u0;
        
    end
    
elseif tag_change==4
    
    
    if t>=time1 && t<time2
        w_val=w1;
        u_val=u1;
    
    elseif t>=time2 && t<time3
        w_val=w2;
        u_val=u2;
        
    elseif t>=time3 && t<time4
        w_val=w3;
        u_val=u3;
        
        
    elseif t>=time4
        w_val=w4;
        u_val=u4;
        
        
    else
        w_val=w0;
        u_val=u0;
        
    end
    
elseif tag_change==5
    
    
    if t>=time1 && t<time2
        w_val=w1;
        u_val=u1;
    
    elseif t>=time2 && t<time3
        w_val=w2;
        u_val=u2;
        
    elseif t>=time3 && t<time4
        w_val=w3;
        u_val=u3;
        
        
    elseif t>=time4 && t<time5
        w_val=w4;
        u_val=u4;
        
    elseif t>=time5
        w_val=w5;
        u_val=u5;
            
    else
        w_val=w0;
        u_val=u0;
        
    end    
else
    w_val=w0;
    u_val=u0;
    
end

beta_val=beta;
tau_val=tau;

beta_f=beta_val/(1+(m/K_p)^n_h)*(1+w_val)*(1-u_val);
s_dot= -beta_f*s*I/n_pop;

i_dot= beta_f*s*I/n_pop - (alpha+eta)*I;

q_plus_r_dot= eta*I ;

m_dot=(eta*I - m)/tau_val; 

xdot=[s_dot;i_dot;q_plus_r_dot;m_dot ];

