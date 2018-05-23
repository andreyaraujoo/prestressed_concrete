tic
clear
%Homework: Prestressed Concrete Design
%Student:Andrey Araújo dos Santos
%Profesor: Dr. Ivan Araújo
%Course: Prestressed Concrete
%Degree: Civil Engineering
%University: Universidade Federal da Integração Latino Americana
%This work is licensed under the Creative Commons Attribution-ShareAlike 3.0 Unported License. 
%To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/3.0/ or send a letter to
%Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Geometria   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_b1 =50; g_b2 =30; g_b3 =60;
g_h1 =15; g_h2 =80; g_h3 =15; g_h4 =20; g_ht =g_h4+g_h3+g_h2+g_h1;
%%%%% Secao 1 %%%%%
g_a1=((g_b1+g_b2)*g_h1)/2;
g_y1=(((g_h1/3)*(g_b1+2*g_b2))/(g_b1+g_b2));
g_a1y1=g_a1*g_y1;
g_i01=((g_h1^3)*(g_b2^2+4*g_b2*g_b1+g_b1^2))/(36*(g_b1+g_b2));
%%%%% secao 2 %%%%%
g_a2=(g_b2*g_h2);
g_y2=g_h1+g_h2/2;
g_a2y2=g_a2*g_y2;
g_i02=(g_b2*g_h2^3)/12;
%%%%% secao 3 %%%%%
g_a3=((g_b2+g_b3)*g_h3)/2;
g_y3=g_h1+g_h2+(((g_h3/3)*(g_b2+2*g_b3))/(g_b2+g_b3));
g_a3y3=g_a3*g_y3;
g_i03=((g_h3^3)*(g_b2^2+4*g_b2*g_b3+g_b3^2))/(36*(g_b3+g_b2));
%%%%% secao 4 %%%%%
g_a4=(g_b3*g_h4);
g_y4=g_h1+g_h2+g_h3+g_h4/2;
g_a4y4=g_a4*g_y4;
g_i04=(g_b3*g_h4^3)/12;
%%%%%%%%%%% Total %%%%%%%%%%%
g_at=g_a1+g_a2+g_a3+g_a4;
g_ayt=g_a4y4+g_a3y3+g_a2y2+g_a1y1;
g_yinf=(g_a1*(g_h4+g_h3+g_h2+g_h1/2)+g_a2*(g_h4+g_h3+g_h2/2)+g_a3*(g_h4+g_h3/2)+g_a4*(g_h4/2))/g_at;
g_ysup=g_h1+g_h2+g_h3+g_h4-g_yinf;
g_ad1=g_a1*((g_y1)-g_ysup)^2;
g_ad2=g_a2*(g_h1+(g_h2/2)-g_ysup)^2;
g_ad3=g_a3*(g_h1+g_h2+((((g_h3/3)*(g_b2+2*g_b3))/(g_b2+g_b3)))-g_ysup)^2;
g_ad4=g_a4*(g_h1+g_h2+g_h3+(g_h4/2)-g_ysup)^2;
g_i0t=g_i04+g_i03+g_i02+g_i01+g_ad1+g_ad2+g_ad3+g_ad4;
wicsup=g_i0t/g_ysup;
wicinf=g_i0t/g_yinf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Coordenadas Iniciais   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0cb1=0;    x0cb2=0;    x0cb3=1;  x0cb4=3.0;
y0cb1=0.3;  y0cb2=0.9;  y0cb3=1.3;  y0cb4=1.3;
xmcb1=14.5; xmcb2=14.5; xmcb3=14.5; xmcb4=14.5;
ymcb1=0.08;  ymcb2=0.08;  ymcb3=0.2;  ymcb4=0.2;
acb1 = ((y0cb1-ymcb1)/((x0cb1-xmcb1)^2));
acb2 = ((y0cb2-ymcb2)/((x0cb2-xmcb2)^2));
acb3 = ((y0cb3-ymcb3)/((x0cb3-xmcb3)^2));
acb4 = ((y0cb4-ymcb4)/((x0cb4-xmcb4)^2));
xcb1 = [x0cb1,1,2,3,4,5,6,7,8,9,10,11,12,13,14.5];
xcb2 = [x0cb2,1,2,3,4,5,6,7,8,9,10,11,12,13,14.5];
xcb3 = [x0cb3 ,2,3,4,5,6,7,8,9,10,11,12,13,14.5];
xcb4 = [x0cb4 ,4,5,6,7,8,9,10,11,12,13,14.5];
ycb1= acb1* ((xcb1-xmcb1).^2) +ymcb1;
ycb2= acb2* ((xcb2-xmcb2).^2) +ymcb2;
ycb3= acb3* ((xcb3-xmcb3).^2) +ymcb3;
ycb4= acb4* ((xcb4-xmcb4).^2) +ymcb4;
for i=1:length(xcb1)
    yln(i)=g_yinf/100;
    ex_cb1(i)=(g_yinf/100-ycb1(i));
end
for i=1:length(xcb2)
    ex_cb2(i)=(g_yinf/100-ycb2(i));
end
for i=1:length(xcb3)
    ex_cb3(i)=(g_yinf/100-ycb3(i));
end
for i=1:length(xcb4)
    ex_cb4(i)=(g_yinf/100-ycb4(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   Calculo de perdas   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   Propriedades dos materiais   %%%%%%%%%%%%%%%%%%%%%%%
fptk =  1900000;    %kPa
fpyk =  1710000;    %kPa
ep   =  200000000;  %kPa
fck  =  30000;       %kPa
sigmapipre = min(0.77*fptk,0.85*fpyk); %Kpa(Pré-tração)
sigmapipos = min(0.74*fptk,0.82*fpyk); %Kpa(Pós-tração)
mi = 0.05; %cordoalha engraxada
k = 0.01*mi;
deltaw =2;
cpiiis=0.38;
n=2;
t12=14;t34=28;%dias
beta1_14=exp(cpiiis*(1-(28/t12)^0.5));
beta1_28=exp(cpiiis*(1-(28/t34)^0.5));
eci_14=1000*5600*((fck/1000)*beta1_14)^0.5;%kpa
eci_28=1000*5600*((fck/1000)*beta1_28)^0.5;%kpa
eci=(2*eci_14 + 2*eci_28)/4;
alfap=ep/eci;
ep_gp1 = g_yinf-ymcb1;
ep_gp2 = g_yinf-ymcb3;
y0=(ymcb1+ymcb2+ymcb3+ymcb4)/(4);%m
excp=(g_yinf/100)-y0;%m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  Resumo cabos %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncordcb1 = 7;ncordcb2 = 7;ncordcb3 = 7; ncordcb4 = 7;
acordcb1 = 1.4;
acordcb2 = 1.4;
acordcb3 = 1.4;
acordcb4 = 1.4;
ascb1= acordcb1*ncordcb1;ascb2= acordcb2*ncordcb2;ascb3= acordcb3*ncordcb3;ascb4= acordcb4*ncordcb4;
p0cb1 = sigmapipos*ascb1/100^2; p0cb2 = sigmapipos*ascb2/100^2;p0cb3 = sigmapipos*ascb3/100^2;p0cb4 = sigmapipos*ascb4/100^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Perdas por atrito     %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sumalphacb1 = zeros(1,15);
sumalphacb2 = zeros(1,15);
sumalphacb3 = zeros(1,14);
sumalphacb4 = zeros(1,12);
picb1 = zeros(1,15);
picb2 = zeros(1,15);
picb3 = zeros(1,14);
picb4 = zeros(1,12);
for i=2:1:15
    sumalphacb1(1)=0;
    sumalphacb1(i) = sumalphacb1(i-1) + (2*(ycb1(i-1)-ycb1(i))/(xcb1(i)-xcb1(i-1)));
    picb1(1)= p0cb1 * exp(-(mi*sumalphacb1(1)+k*xcb1(1)));
    picb1(i)= p0cb1 * exp(-(mi*sumalphacb1(i)+k*xcb1(i)));
end
for i=2:1:15
    sumalphacb2(1)=0;
    sumalphacb2(i) = sumalphacb2(i-1) + (2*(ycb2(i-1)-ycb2(i))/(xcb2(i)-xcb2(i-1)));
    picb2(1)= p0cb2 * exp(-(mi*sumalphacb2(1)+k*xcb2(1)));
    picb2(i)= p0cb2 * exp(-(mi*sumalphacb2(i)+k*xcb2(i)));
end
for i=2:1:14
    sumalphacb3(1)=0;
    sumalphacb3(i) = sumalphacb3(i-1) + (2*(ycb3(i-1)-ycb3(i))/(xcb3(i)-xcb3(i-1)));
    picb3(1)= p0cb3 * exp(-(mi*sumalphacb3(1)+k*xcb3(1)));
    picb3(i)= p0cb3 * exp(-(mi*sumalphacb3(i)+k*xcb3(i)));
end

for i=2:1:12
    sumalphacb4(1)=0;
    sumalphacb4(i) = sumalphacb4(i-1) + (2*(ycb4(i-1)-ycb4(i))/(xcb4(i)-xcb4(i-1)));
    picb4(1)= p0cb4 * exp(-(mi*sumalphacb4(1)+k*xcb4(1)));
    picb4(i)= p0cb4 * exp(-(mi*sumalphacb4(i)+k*xcb4(i)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Perdas por ancoragem %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stopcb1 = (deltaw/1000)*ep*(ascb1/10000);
stopcb2 = (deltaw/1000)*ep*(ascb2/10000);
stopcb3 = (deltaw/1000)*ep*(ascb3/10000);
stopcb4 = (deltaw/1000)*ep*(ascb4/10000);
areapacb1=cumtrapz(picb1);
areapacb2=cumtrapz(picb2);
areapacb3=cumtrapz(picb3);
areapacb4=cumtrapz(picb4);
%%%%%%   descontar (pi+1)*deltax %%%%%%%%%
for i=2:+15
    mcb1(i-1)=(picb1(i-1)-picb1(i))/(xcb1(i)-xcb1(i-1));mcb1(15)=0;
    mcb2(i-1)=(picb2(i-1)-picb2(i))/(xcb2(i)-xcb2(i-1));mcb2(15)=0;
end
for i=2:+14
    mcb3(i-1)=(picb3(i-1)-picb3(i))/(xcb3(i)-xcb3(i-1));mcb3(14)=0;
end
for i=2:+12
    mcb4(i-1)=(picb4(i-1)-picb4(i))/(xcb4(i)-xcb4(i-1));mcb4(12)=0;
end
%%%%%%%%%%% Calcular curvas de ancoragem %%%%%%%%%%%%%%%%%%
areaanccb1(1)=0;
icb1=2;
while areaanccb1<stopcb1
    icb1=icb1+1;
    areaanccb1(icb1-1)=(areapacb1(icb1)-picb1(icb1)*(xcb1(icb1)-xcb1(1)))*2;
    ccb1(icb1-1) = areaanccb1(icb1-1)-stopcb1;
    bcb1(icb1-1) = 2*mcb1(icb1-1)*(xcb1(icb1)-xcb1(1));
    acb1(icb1-1) = mcb1(icb1-1);
end
pcb1 = [acb1(icb1-1) bcb1(icb1-1) ccb1(icb1-1)];
rpcb1 = roots(pcb1);
for i=1:icb1-1
    x_anccb1(i)=+i-1;
    p_anccb1(i) = picb1(i) - 2*(picb1(i)-picb1(icb1-1))+mcb1(icb1-1)*rpcb1(2);
end
areaanccb2(1)=0;
icb2=2;
while areaanccb2<stopcb2
    icb2=icb2+1;
    areaanccb2(icb2-1)=(areapacb2(icb2)-picb2(icb2)*(xcb2(icb2)-xcb2(1)))*2;
    ccb2(icb2-1) = areaanccb2(icb2-1)-stopcb2;
    bcb2(icb2-1) = 2*mcb2(icb2-1)*(xcb2(icb2-1)-xcb2(1));
    acb2(icb2-1) = mcb2(icb2-1);
end
pcb2 = [acb2(icb2-1) bcb2(icb2-1) ccb2(icb2-1)];
rpcb2 = roots(pcb2);
for i=1:icb2-1
    x_anccb2(i)=+i-1;
    p_anccb2(i) = picb2(i) - 2*(picb2(i)-picb2(icb2-1))+mcb2(icb2-1)*rpcb2(2);
end
areaanccb3(1)=0;
icb3=2;
while areaanccb3<stopcb3
    icb3=icb3+1;
    areaanccb3(icb3-1)=(areapacb3(icb3)-picb3(icb3)*(xcb3(icb3)-xcb3(1)))*2;
    ccb3(icb3-1) = areaanccb3(icb3-1)-stopcb3;
    bcb3(icb3-1) = 2*mcb3(icb3-1)*(xcb3(icb3-1)-xcb3(1));
    acb3(icb3-1) = mcb3(icb3-1);
end
pcb3 = [acb3(icb3-1) bcb3(icb3-1) ccb3(icb3-1)] ;
rpcb3 = roots(pcb3);
for i=1:icb3-1
    x_anccb3(i)=+i;
    p_anccb3(i) = picb3(i) - 2*(picb3(i)-picb3(icb3-1))+mcb3(icb3-1)*rpcb3(2);
end
areaanccb4(1)=0;
icb4=2;
while areaanccb4<stopcb4
    icb4=icb4+1;
    areaanccb4(icb4-1)=(areapacb4(icb4)-picb4(icb4)*(xcb4(icb4)-xcb4(1)))*2;
    ccb4(icb4-1) = areaanccb4(icb4-1)-stopcb4;
    bcb4(icb4-1) = 2*mcb4(icb4-1)*(xcb4(icb4-1)-xcb4(1));
    acb4(icb4-1) = mcb4(icb4-1);
end
pcb4 = [acb4(icb4-1) bcb4(icb4-1) ccb4(icb4-1)];
rpcb4 = roots(pcb4);
for i=1:icb4-1
    x_anccb4(i)=+i+2;
    p_anccb4(i) = picb3(i) - 2*(picb4(i)-picb4(icb4-1))+mcb4(icb4-1)*rpcb4(2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Perdas por encurtamento %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%agrupar valores
for i=1:length(xcb1)
    if (i<(icb1-1))
        pi_enc1(i)=p_anccb1(i);
    else
        pi_enc1(i)=picb1(i);
    end
end
for i=1:length(xcb2)
    if (i<(icb2-1))
        pi_enc2(i)=p_anccb2(i);
    else
        pi_enc2(i)=picb2(i);
    end
end
for i=1:length(xcb3)
    if (i<(icb3-1))
        pi_enc3(i)=p_anccb3(i);
    else
        pi_enc3(i)=picb3(i);
    end
end
for i=1:length(xcb4)
    if (i<(icb4-1))
        pi_enc4(i)=p_anccb4(i);
    else
        pi_enc4(i)=picb4(i);
    end
end

for i=1:length(xcb1)
    enc_m(i)=-6.1*xcb1(i)^2+181.17*xcb1(i);
end

for i=1:length(xcb1)
    asexcab1(i)=ascb1*ex_cb1(i);
    asexcab2(i)=ascb2*ex_cb2(i);
    asexcab3(1)=0;
    
    if i>1
        
        asexcab3(i)=ascb3*ex_cb3(i-1);
        
    end
    asexcab4(1)=0;
    asexcab4(2)=0;
    asexcab4(3)=0;
    if i>3
        asexcab4(i)=ascb4*ex_cb4(i-3);
    end
    
    sumasep(i)=(asexcab1(i)+asexcab2(i)+asexcab3(i)+asexcab4(i))/(ascb1+ascb2+ascb3+ascb4);
end

for i=1:length(xcb1)
    sigmacp(i)=(-pi_enc1(i)-pi_enc2(i))/(g_at*100^-2)+((-pi_enc1(i)*ex_cb1(i)-pi_enc2(i)*ex_cb2(i))/g_i0t*100^-4)*(sumasep(i));
    if i==3
        sigmacp(i)=(-pi_enc1(i)-pi_enc2(i)-pi_enc3(i-1))/(g_at*100^-2)+((-pi_enc1(i)*ex_cb1(i)-pi_enc2(i)*ex_cb2(i)-pi_enc3(i-1)*ex_cb3(i-1))/g_i0t*100^-4)*(sumasep(i));
    elseif i>3
        sigmacp(i)=(-pi_enc1(i)-pi_enc2(i)-pi_enc3(i-1)-pi_enc4(i-3))/(g_at*100^-2)+((-pi_enc1(i)*ex_cb1(i)-pi_enc2(i)*ex_cb2(i)-pi_enc3(i-1)*ex_cb3(i-1)-pi_enc4(i-3)*ex_cb4(i-3))/g_i0t*100^-4)*(sumasep(i));
    end
    sigmacg(i)=enc_m(i)*sumasep(i)/(g_i0t*100^-4);
    deltasigma(i)=(n-1)/(2*n)*alfap*(sigmacp(i)-sigmacg(i));
end
for i=1:length(xcb1)
   enc_deltapcb1(i)=deltasigma(i)*ascb1*100^-2;
   pf_enc1(i)=pi_enc1(i)+enc_deltapcb1(i);
end
for i=1:length(xcb2)
       enc_deltapcb2(i)=deltasigma(i)*ascb2*100^-2;
          pf_enc2(i)=pi_enc2(i)+enc_deltapcb2(i);

end
for i=1:length(xcb3)
    enc_deltapcb3(i)=deltasigma(i+1)*ascb3*100^-2;
       pf_enc3(i)=pi_enc3(i)+enc_deltapcb3(i);

end
for i=1:length(xcb4)
        enc_deltapcb4(i)=deltasigma(i+2)*ascb4*100^-2;
           pf_enc4(i)=pi_enc4(i)+enc_deltapcb4(i);


end




%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Gráfico    %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
subplot(3,1,1)
plot(xcb1,yln,"k",xcb1, ycb1,"g", xcb2, ycb2,"b", xcb3, ycb3,"m", xcb4, ycb4,"r" )
axis([0 14.5 0 1.3])
label = ['Traçado dos Cabos'];
title(label)
ylabel('h (m)')
xlabel('x (m)')
legend('Y_{cinf}','Cabo 1','Cabo 2', 'Cabo 3','Cabo 4')
lgd.NumColumns = 2;
subplot(3,2,3)
plot(xcb1, picb1,"g", x_anccb1,p_anccb1,'--g',xcb1,pf_enc1,"-.g")
% axis([0 14.5 0 2500])
label = ['Perdas imediatas cabo 1'];
title(label)
ylabel('P0(kN)')
xlabel('x (m)')
subplot(3,2,4)
plot(xcb2, picb2,"b", x_anccb2,p_anccb2,'--b',xcb2,pf_enc2,"-.b")
% axis([0 14.5 0 2500])
label = ['Perdas imediatas cabo 2'];
title(label)
ylabel('P0(kN)')
xlabel('x (m)')
subplot(3,2,5)
plot(xcb3, picb3,"m", x_anccb3,p_anccb3,'--m',xcb3,pf_enc3,"-.m")
% axis([0 14.5 0 2500])
label = ['Perdas imediatas cabo 3'];
title(label)
ylabel('P0(kN)')
xlabel('x (m)')
subplot(3,2,6)
plot(xcb4, picb4,"r", x_anccb4,p_anccb4,'--r',xcb4,pf_enc4,"-.r")
% axis([0 14.5 0 2500])
label = ['Perdas imediatas cabo 4'];
title(label)
ylabel('P0(kN)')
xlabel('x (m)')

toc