function [ qchfgroen ] = qchfgroen(ichan,jlev,z,de,p,g,xaeq,L,rhof,rhog)
%Module Groeneveld
xe_lim=[-5.00,-0.50,-0.40,-0.30,-0.20,-0.15,-0.10,-0.05,0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.60,0.70,0.80,0.90,1.0,5.0];
xe_interval=[1, 2, 2, 3, 3, 4, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24];
g_lim=[-50,0,50,100,300,500,750,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000, 40000000];
g_interval=[1, 2, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9,   9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11,  11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13,13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 22]
p_lim=[0,100,300,500,1000,2000,3000,5000,7000,10000,12000,14000,16000,18000,20000,21000,40000000]
p_interval=[1, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16]
%Setting equilibrium quality limits for the reduced 2x2x2 n1 = xe_interval(max(1,min(size(xe_interval),ceiling(xaeq/0.05)+11))) xelow = xelim(n1) xeup = xelim(n1+1) ! Setting mass-flux limits for the reduced 2x2x2 n2 = g_interval(min(size(g_interval),ceiling(gsi/50)+1)) glow = glim(n2) gup = glim(n2+1) ! Setting pressure limits for the reduced 2x2x2 n3 = p_interval(min(size(p_interval),ceiling(psi/100))) plow = plim(n3) pup = plim(n3+1) ! Develops 2x2x2 matrix from original 17x23x25 using pressure, xe and ! mass flux boundaries Gro_red = real(Gro(n2:n2+1,n1:n1+1,n3:n3+1))
n1=xe_interval(max(1,min(size(xe_interval),real(xaeq/0.05)+11)))
xelow = zeros(n1)
xeup = zeros(n1+1)

%Setting mass-flux limits for the reduced 2x2x2
n2 = g_interval(min(size(g_interval),real(gsi/50)+1))
glow = zeros(n2)
gup = zeros(n2+1)

%Setting pressure limits for the reduced 2x2x2
n3 = p_interval((min(size(p_interval))),floor(psi/100))
plow = zeros(n3)
pup = zeros(n3+1)

Gro_red = real(Gro(n2:n2+1,n1:n1+1,n3:n3+1))
Gro_red=zeros(2,2,2);
Gro_redp=zeros(2,2);
Gro_redpg=zeros(2);
Gro=zeros(23,25,17);
Gro=reshape(chf_array(23,25,17),[2,1,3]);

%Conversion units
t_lbm_kg=0.453592;
t_hr_s=3600;
t_ft_m=0.3048;
t_psi_MPa=0.00689476;
t_MPa_kPa=1000;

gsi = abs(g)*t_lbm_kg/(t_hr_s*t_ft_m^2) %[kg/m^2s] 
psi = p*t_psi_MPa*t_MPa_kPa %[kPa]
desi = de*t_ft_m %[m] 
zsi = z*t_ft_m %[m] 
Lsi = L*t_ft_m %[m]

%Develops Matrix from 17x23x25 using pressure , xe and mass flux boundaries
%Setting equilibrium quality limits for the reduced
n1=xe_interval(max(1,min(size(xe_interval),ge(xaeq/0.05)+11)));
    xelow=xelim(n1);
    xeup=xelim(n1+1);
    
%Setting mass-flux limits for the reduced
n2 = g_interval(min(size(g_interval),ge(gsi/50)+1));
glow = glim(n2);
gup = glim(n2+1);
%Setting pressure limits for the reduced
n3 = p_interval(min(size(p_interval),ge(psi/100))) 
plow = plim(n3)
pup = plim(n3+1)

Gro_red = real(Gro(n2:n2+1,n1:n1+1,n3:n3+1))
%Develops 1x2x2 matrix reduced
Gro_redp(1,1) = (Gro_red(1,1,2)-Gro_red(1,1,1))/(pup-plow)*(psi-plow)+Gro_red(1,1,1) 
Gro_redp(1,2) = (Gro_red(1,2,2)-Gro_red(1,2,1))/(pup-plow)*(psi-plow)+Gro_red(1,2,1) 
Gro_redp(2,1) = (Gro_red(2,1,2)-Gro_red(2,1,1))/(pup-plow)*(psi-plow)+Gro_red(2,1,1)
Gro_redp(2,2) = (Gro_red(2,2,2)-Gro_red(2,2,1))/(pup-plow)*(psi-plow)+Gro_red(2,2,1)

%Develops 1x2 matrix reduced/interpolated in pressure and mass flux
Gro_redpg(1) = (Gro_redp(2,1)-Gro_redp(1,1))/(gup-glow)*(gsi-glow)+ Gro_redp(1,1) 
Gro_redpg(2) = (Gro_redp(2,2)-Gro_redp(1,2))/(gup-glow)*(gsi-glow)+ Gro_redp(1,2)

%Final interpolated chf value, reduced/interpolated in pressure, ! mass flux and xe 
Gro_int = (Gro_redpg(2)-Gro_redpg(1))/(xeup-xelow)*(xaeq-xelow)+ Gro_redpg(1)
qint = 1000*Gro_int*((1./t_btu_J)/(1./t_hr_s))/((1./t_ft_m)^2)

%Subchannel cross-section geometry factor 
    if (desi > 0.003) & (desi < 0.025) ;
    K1 = (0.008/desi)^0.5 
    elseif (desi >= 0.025)
    K1 = 0.57 
    else
    K1=1.633
    end 
    % Heated length factor
    if (xaeq>0);
        alpha=xaeq*rhof/(xaeq*rhof+(1-xaeq)*rhog);
    else 
        alpha=0;
    end
    if (Lsi/desi>5);
        K4=exp(desi/Lsi)*exp(2*alpha);
    else
        K4=1.0;
    end
    for j=1:jmax;
        if (xe(j)>=0);
            xe(j)=xe(j);
        else
            xe(j)=0;
        end
    end
        K5=1;
        if  (osb_exists)==True;
            for jh=osb:j:jmax;
                sum_q=0;
                num_sum=0;
                for j_int = osb_j:jh;
                    sum_q=sum_q+q(j_int);
                    num_sum=num_sum+1;
                end
                if (jh~=osb_j)
                    qbla=sum_q/num_sum;
                else
                    qbla=sum_q;
                end
                if q(jh)>0;
                    K5(jh)=q(jh)/qbla;
                else
                    K5(jh)=1.0;
                end
            end
        end
         qbla=K1*K4*K5*qint;
end