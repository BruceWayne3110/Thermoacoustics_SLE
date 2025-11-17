% Code to drive stuart landau oscillator using mmpp
%% code to write the timeseries into a file.
writematrix(qdotprimemodel,'hrrfcombnoise.txt');
%% code to solve stuart landau kicked oscillator with mmpp governing the
% parameter mu. 
clc;clear;
set(0, 'defaultFigureWindowState', 'maximized');
set(groot, 'defaultAxesBox', 'on');
set(groot, 'defaultTextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultColorbarTickLabelinterpreter','latex'); 
set(groot, 'defaultLineLineWidth', 1.5);
sigma=1;    
amp=0.0035; 
mu=0.4;     
lambda_array=[105.9843,86.39,68.80,53.21,39.62]./2;
gr=1;    
w1=1;  
gi=0;    
T=1200;   h=10^-4;
[u,t,v,mmpp]=rk4forstuartlandau(T,mu,sigma,gr,gi,w1,amp,lambda_array,h,0.0001,0.0001);
clearvars -except u t v mu;
% Amplitude=u.^2+v.^2;
% qdotprimemodel=3*(u)./5+0.2;
qdotprimemodel=(u-mean(u));
tsmodel=t./1000;
figure();
plot(tsmodel,qdotprimemodel,'k');
figure();
% plot(tsmodel,Amplitude);
% xlabel('$t$','FontSize',22,'interpreter','latex')
% ylabel('$u$','FontSize',22,'interpreter','latex')
% title(mu);
% set(gca,'linewidth',2,'FontSize',15)
% saveas(gcf,'fixed point unf','fig');
% saveas(gcf,'fixed point unf','png');                                                                                                                                                                                                          
%% Bifurcation plot function
% clc;clear;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          c;clear;
N=51;
sigma=1;  %4
amp=0.0035; %0.25
mu0=0.02;  %1
lambda_array=[105.9843,86.39,68.80,53.21,39.62]./2;
gr=1;    %4
w1=1;  %100
gi=0;    %4
T=1200;h=10^-4;
mu_arr=linspace(-0.5,1,N); 
mu_arr2=fliplr(mu_arr);
u_all = zeros(10000001, length(mu_arr));
u_all2 = zeros(10000001, length(mu_arr));
% kicks=zeros(T/h,length(mu_arr));
parpool(15);
parfor i=1:N   
   [u,t,v,mmpp]=rk4forstuartlandau(T,mu_arr(i),sigma,gr,gi,w1,amp,lambda_array,h,0.0001,0.0001);
   qdotprimemodel=u-mean(u);
    tsmodel=t./1000;
   u_all(:, i) = qdotprimemodel(end-10000000:end);
   u_all2(:, i) = u(end-10000000:end);
   % kicks(:,i)=mmpp;
end
% parfor i=1:N
%    [u,t,v]=rk4forstuartlandau(T,mu_arr(i),sigma,gr,gi,w1,0,lambda_array,tpm,h,0.0001,0.0001);
%    qdotprimemodel=1.5*(1+u+mu_arr(i)/2)./4;
%     tsmodel=t./1000;
%    u_all2(:, i) = qdotprimemodel(end-10000000:end);
% end
u_amp=rms(u_all,1);
u_amp2=rms(u_all2,1);
figure();
scatter(mu_arr,u_amp,'k');
hold on 
scatter(mu_arr,u_amp2,'k');
% numberofkicks=sum(kicks);
% figure()
% plot(mu_arr,numberofkicks);
xlabel('$\mu$','FontSize',22,'interpreter','latex')
ylabel('$u_{rms}$','FontSize',22,'interpreter','latex')
set(gca,'linewidth',2,'FontSize',15)
 % save('continuoustpm_exp3.mat',"-v7.3")
% hold on;
% scatter(mu_arr2,u_amp2,'k');
% saveas(gcf,'bif_rk4_mmpp0.05','fig');
% saveas(gcf,'bif_rk4_mmpp0.05','png');
% u2=[mu_arr',u_all'];
% u3=[mu_arr2',u_all2'];
% writematrix(u2,'bifurcation_rk4_mmpp_fwd0.05.txt');
% writematrix(u3,'bifurcation_rk4_mmpp_back0.05.txt');
%% Function for generating kicks according to MMPP
function newtimearray=MMPP(k,time,lambda_array,h) 
%calculation of stationary state probability.
[V,D,W]=eig(k);
P0=W(:,1)/sum(W(:,1));
Q=k.*lambda_array'-diag(lambda_array); % Infinitesimal generator matrix


n=size(k,2); % order of the matrix

MMPP_state_array_size = 0;
% MMPP_samples_time_array = 0;
%MMPP_samples_state_array = 0;
initial_state_prob = 0;
initial_state_prob = 1/n; 
initial_state = 0;
i = 1;
for j = 1:n
    initial_state(i,j)=j;
    initial_state(i+1,j)=initial_state_prob;
end

current_state = randsrc (1,1,initial_state);
current_time = 0;
next_interrupt_time = 0;
next_sojourn_time=0;
i = 1;
%Calculating the sojourn times and the time of arrivals according to it
sojourn_times(i)=exprnd(-1/Q(current_state,current_state));
MMPP_samples_time_array (i) = current_time;
j=1;
while (next_sojourn_time<=time)
    i=i+1;
    current_sojourn= next_sojourn_time;
    current_state = randsrc (1,1,[initial_state(1,:);k(current_state,:)]);
    MMPP_samples_state_array (i) = current_state;
    while (next_interrupt_time <=next_sojourn_time|| next_interrupt_time<=time)
           j=j+1;
           current_time = next_interrupt_time;
           MMPP_samples_time_array (j) = current_time;
           next_interrupt_time = current_time + exprnd(1/lambda_array(current_state));
    end
    
    next_sojourn_time=next_sojourn_time+exprnd(-1/Q(current_state,current_state));
    soujorn_times(i)=sojourn_times(end)+next_sojourn_time;
end
i=1;
newtimearray(i)=0;
time_steps(i)=0;
while time_steps(end)<=time
    i=i+1;
    time_steps(i)=time_steps(i-1)+h;
    for j=1:size(MMPP_samples_time_array,2)
        if (MMPP_samples_time_array(j)-time_steps(i)<0) && (MMPP_samples_time_array(j)-time_steps(i-1)>0)
              newtimearray(i)=1;
         % else 
         %       newtimearray(i)=0;
        end    
    end    
end
newtimearray=[newtimearray,zeros(1,time/h-length(newtimearray))];
r=sum(newtimearray==1);
end
function out=f1(mu,sigma,w1,u,v,gr,gi,t)
out=mu*sigma*u-w1*v-gr*u^3-gr*u*v^2;
end
function out=f2(mu,sigma,w1,u,v,gr,gi)
out=w1*u+mu*sigma*v-gr*v^3-gr*u^2*v;
end

function [u,t,v,mmpp]=rk4forstuartlandau(T,mu,sigma,gr,gi,w1,amp,arrival_rates,h,u0,v0)
tpm=[1-exp(-1*(10*mu+5)),exp(-1*(10*mu+5))/4,exp(-1*(10*mu+5))/4,exp(-1*(10*mu+5))/4,exp(-1*(10*mu+5))/4;
    1-exp(-1*(10*mu+5)),exp(-1*(10*mu+5))/4,exp(-1*(10*mu+5))/4,exp(-1*(10*mu+5))/4,exp(-1*(10*mu+5))/4;
    1-exp(-1*(10*mu+5)),exp(-1*(10*mu+5))/4,exp(-1*(10*mu+5))/4,exp(-1*(10*mu+5))/4,exp(-1*(10*mu+5))/4;
    1-exp(-1*(10*mu+5)),exp(-1*(10*mu+5))/4,exp(-1*(10*mu+5))/4,exp(-1*(10*mu+5))/4,exp(-1*(10*mu+5))/4;
    1-exp(-1*(10*mu+5)),exp(-1*(10*mu+5))/4,exp(-1*(10*mu+5))/4,exp(-1*(10*mu+5))/4,exp(-1*(10*mu+5))/4;];
mmpp=MMPP(tpm,T,arrival_rates,h);
N=ceil(T/h);
dW=sqrt(h).*randn(1,N);
u=u0;
v=v0;
t=0;
for i=1:N
      m1=h*f2(mu,sigma,w1,u(i),v(i),gr,gi);
      l1=h*f1(mu,sigma,w1,u(i),v(i),gr,gi,t(i));

      a=u(i)+l1*0.5;
      b=v(i)+m1*0.5;
      t2=t(i)+h*0.5;

      m2=h*f2(mu,sigma,w1,a,b,gr,gi);
      l2=h*f1(mu,sigma,w1,a,b,gr,gi,t(i));

      a=u(i)+l2*0.5;
      b=v(i)+m2*0.5;
      t2=t(i)+h*0.5;

      m3=h*f2(mu,sigma,w1,a,b,gr,gi);
      l3=h*f1(mu,sigma,w1,a,b,gr,gi,t(i));

      a=u(i)+l3;
      b=v(i)+m3;
      t2=t(i)+h;

      m4=h*f2(mu,sigma,w1,a,b,gr,gi);
      l4=h*f1(mu,sigma,w1,a,b,gr,gi,t(i));

      dx=(l1+l2*2+l3*2+l4)/6;
      dy=(m1+m2*2+m3*2+m4)/6;

      u(i+1)=u(i)+dx+mmpp(i)*(amp);
      v(i+1)=v(i)+dy+mmpp(i)*(amp);
      t(i+1)=t(i)+h;

end    
end