clc;
clear all;
axonmain();

function [t,x]=axonmain()
global P % parameters
setUp;
[t,x]=ode15s(@Dyn,P.tD,P.Xo);
plotResults(t,x); % axon
end
% dynamic function
function Dx=Dyn(t,x)
% setProcess: Establish G matrix & other values for the process
%Set the drive matrix as well
global P                                                	% param&I-stim
N=length(x)/4; G=zeros(N); f=zeros(N,1);                	% numberOnodes
V=x(1:N);m=x(N+1:2*N);h=x(2*N+1:3*N);n=x(3*N+1:4*N);    	% extractStates
for i=1:N                                               	% loop eachNode
  ga=P.xSectArea(i)/((P.Ra*P.dx(i)));                   	% axon conduct
  gL=P.surfArea(i)/(P.RL);                              	% chloride/leak
  gNaMax=P.gMax_Na*P.surfArea(i)*P.isActive(i);         	% Na maxConduct
  gKMax =P.gMax_K *P.surfArea(i)*P.isActive(i);         	% K maxConduct
  C=P.C*P.surfArea(i);                                  	% capacitance-F
  gNa=m(i)^3*h(i)*gNaMax;                               	% conductWgates
  gK =n(i)^4*gKMax;                                     	% conductWgates
  a = 0.8e-4;
  Rout = 24e4;
  gout = 1/(Rout/P.xSectArea(i));
  SUMgChan=gNa+gK+gL;                                   	% sum memb g's
  f(i)=(gNa*P.ENa)+(gK*P.EK)+(gL*P.ECl);                	% sum drive
  if N==1
  	G = -SUMgCHan;
  	f = f+interp1(P.Io(:,1),P.Io(:,2),t);
else
  	if i == 1
     	G(1, 1) = - ga - SUMgChan ;
     	G(1, 2) = ga;
     	f(1) = f(1)+interp1(P.Io(:,1),P.Io(:,2),t);
  	elseif i == N
      	G(N,N) = -(SUMgChan+ga);
      	G(N,N-1)= ga;
  	else
      	G(i,i) = -(2*ga+SUMgChan);
      	G(i,i-1) = ga;
      	G(i,i+1) = ga;
  	end
  end
 
 
end                                                    	% endForEachNod
DV=(G*V+f)./C;                                          	% curentBalance
[mA,nA,hA,mB,nB,hB]=getRates(V,P.Vrest);                	% alpha & betas
rn = randi(4,i,1);
if (rn(i)==2)
	nB = 0;
end
Dm=(-(mA+mB).*m+mA);Dh=(-(hA+hB).*h+hA);Dn=(-(nA+nB).*n+nA);% gates probs
Dx=[DV;Dm;Dh;Dn];                                       	% asembleStates

set(P.h,'yData',V);title(num2str(t));drawnow;           	% updateAnims
                                             	% process vars
end
% setUp : intialize variables
function P=setUp(P)
% messages
close all;                                               	% clear
fprintf('No Branching Axon');	% announce
% setup variables
global P
P.tMax=0.004;                               	% simu time sec
P.tD=0:P.tMax/100:P.tMax;                   	% desired output time steps
P.Vrest = -60e-3;                           	% rest voltage
P.C = 1e-6;                                 	% bulk capcitance F/cm^2
P.RL = 15e3;                                	% bulk Cl Resist Ohms*cm^2
P.Ra = 3e3;                                 	% bulk axon resist ohms*cm -changeBak2 3e3
P.gMax_K  = 60e-3;                          	% bulk K conductance S/cm^2
P.gMax_Na = 120e-3;                         	% bulk Na conduct. S/cm^2
P.ENa =  55e-3;                             	% Nernst poten. Na volts
P.EK  = -75e-3;                             	% Nernst poten. K  volts
P.ECl = -40e-3;                            	% Nernst poten. Cl volts -changeBak2 -403e-3
P.N=5;                                     	% # compartments/nodes/Sect
P.L=.007;                                   	% overall length (cm)
% VECTOR for EACH:
Z=ones(P.N,1); z=zeros(P.N,1);              	% utility vectors 1's & 0's
P.dx=Z*(P.L/P.N);                           	% subsection length (cm)
P.a=Z*3e-5;                                 	% radii of all nodes (cm)
P.surfArea=2*pi.*P.a.*P.dx;                 	% cm^2 -VECTOR for EACH
P.xSectArea=pi.*P.a.^2;                     	% cm^2 -VECTOR for EACH
P.isActive=Z;                               	% 1 if active node, else 0
% initial conditions
Vo=P.Vrest*Z; mo=0.1; ho=0.5; no=0.3; % initial gate probability
P.Xo=[Vo;mo*Z;ho*Z;no*Z];                   	% stack initial states ->Xo
% Injected current pulse (Istim = Io):
I=1e-10; t1=1e-4; dur=1e-4; t2=t1+dur; miu=1e-9;% pulse start & length time
P.Io= [0 t1-miu  t1  t2  t2+miu   2*P.tMax; ... % assemble time and
  	[0 0   	1   1   0     	0]*I]';  	% current mag (pico=E-12)
setUpPlots(100,'Axon with 35 Nodes');  % Seup plots ahead of sim
end
% setup figures & animation
function setUpPlots(figNum,titleText)
global P
S=get(0,'ScreenSize');                              	% ask 4 screen size
if ~exist('figNum','var'), figure,                  	% if fig# not input
else figure(figNum); end
Vo=P.Xo(1:P.N);mo=P.Xo(P.N+1:2*P.N);                	% extractInitStates
ho=P.Xo(2*P.N+1:3*P.N);no=P.Xo(3*P.N+1:4*P.N);      	% extractInitStates
% place summary figure on screen:
set(gcf,'name',titleText)                           	% window title
set(gcf,'Position',[.68*S(3) 0 .3*S(3) .9*S(4)]);   	% put figure
% setup Stim (Io) window
subplot(4,1,1);
plot(P.Io(:,1),P.Io(:,2),'.-','linewidth',2);       	% plot Io
hold on; xlim([0 P.tMax]);
xlabel('Sec'); ylabel('Amp'); title('Io');
% setup volt vs node window
subplot(4,1,2);
P.h=plot(1:P.N,P.Xo(1:P.N),'.-','markersize',12);   	% P.h points 2 plot
hold on; xlabel('Node #'); ylabel('Volts');         	%
set(gca,'ytick',[-.3 -.2 -.1 P.Vrest -.025 0 .1 ]); 	%
ylim(.15*[-1 1]);                                   	%
% setup gates vs time window
subplot(4,1,3);
plot(0,mo,'b.-',0,ho,'g.-',0,no,'r.-'); hold on;    	% all gates
xlabel('sec');ylabel('Prob.');ylim([0 1]);xlim([0 P.tMax])%
text(0,.96,'  m','fontsize',6,'fontweight','bold','color','b')
text(0,.90,'  h','fontsize',6,'fontweight','bold','color','g')
text(0,.84,'  n','fontsize',6,'fontweight','bold','color','r')
% setup 3d figure window
subplot(4,1,4);
plot3(zeros(P.N,1),1:P.N,Vo,'.-','markersize',7);   	%
ylabel('Node #'); xlabel('time(s)'); zlabel('Volts');   %
set(gca,'ztick',[-.3 -.2 -.1 P.Vrest -.025 0 .1 ]); 	% display
grid on; xlim([0 P.tMax]); zlim(.15*[-1 1]); hold on;   %
% shrink fonts and turn off box on all:
for i=1:4,
  subplot(4,1,i);
  set(gca,'box','off','fontsize',5);
end
subplot(4,1,2);
fprintf('...Plots & variables set up...Simulating...');
drawnow;
end
function [alpha_m,alpha_n,alpha_h,beta_m,beta_n,beta_h]=getRates(Vm,Vm_rest,plotIt)
V = Vm-Vm_rest; V=V*1e3;                          	% convert2 millivolts
%% rates:
%Input Rate Equation Functions Here
alpha_n = 0.01*((-V+10)./(exp((-V+10)/10)-1))*1000;
beta_n = 0.125*exp(-V/80)*1000;
alpha_m = 0.1*((-V+25)./(exp((-V+25)/10)-1))*1000;
beta_m = 4*exp(-V/18)*1000;
alpha_h = 0.07*exp(-V/20)*1000;
beta_h = 1./(exp((-V+30)/10)+1)*1000;
%% fix near singularities:
if(abs(V-25)<1e-2), alpha_m=1*1e3;   end;         	% @ v=25
if(abs(V-10)<1e-2), alpha_n=0.1*1e3; end;         	% @ v=10
%% plot if asked:
if exist('plotIt','var')                          	% PLOT if variable in
  plot(Vm,alpha_m,'b.:',Vm,beta_m,'b.-', ...
   	Vm,alpha_n,'g.:',Vm,beta_n,'g.-', ...
   	Vm,alpha_h,'r.:',Vm,beta_h,'r.-')
  text(0,800,'m','fontsize',8,'fontweight','bold','color','b')
  text(0,800,'  	h','fontsize',8,'fontweight','bold','color','r')
  text(0,800,'      	n','fontsize',8,'fontweight','bold','color','g')
  title('Hodgekin and Huxley Gate Kinetic Rates');
  ylabel('(1/sec)'); xlabel('Volts')
end
end
% plotResults function: wrap up by ploting final results
function plotResults(t,x);
global P;
% get variables
N=size(x,2)/4; nt=size(x,1); v1=ones(nt,1);             	% #S time,nodes
V=x(:,1:N);m=x(:,N+1:2*N);h=x(:,2*N+1:3*N);n=x(:,3*N+1:4*N);% extractStates
% V,m,h,n
% make plots
 subplot(4,1,2); ;                       	% volt VS node
 plot(1:N,V);                          	% prob VS time
 subplot(4,1,3);                        	% 3D plot t,N,V
 plot(t,m,'-b',t,h,'-g',t,n,'-r')
 subplot(4,1,4);
 for i= 1:N
 	plot3(t,i*v1,x(:,i))
 end
axis auto
sound(sin(1:1000),1e4);           	% annouyingBeep
end
