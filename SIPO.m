%% SIPO algorithm MATLAB Code
clc;
clear all;
close all ;
%%

prompt   = {'Please enter the number of run:'};
title    = 'SIPO Algorithm';
dims     = [1 45];
nline    = 1;
definput = {'100','a'};
answer   = inputdlg(prompt,title,dims,definput)
Run_Num  = answer(1,:);Run_Num=str2num(Run_Num{:});

prompt   = {'maxt','npop' ,'F','Beta','c','m_Ratio'};
title    = 'SIPO parameters';
nline    = 1;
dims     = [1 45];
definput = {'500','25' '1','0.8','2','0.2','a'};
answer   = inputdlg(prompt,title,dims,definput);
maxt     = answer(1,:);maxt=str2num(maxt{:});
npop     = answer(2,:);npop=str2num(npop{:});
F        = answer(3,:);F=str2num(F{:});
Beta     = answer(4,:);Beta=str2num(Beta{:});
c        = answer(5,:);c=str2num(c{:});
m_Ratio  = answer(6,:);m_Ratio=str2num(m_Ratio{:});

n  = 0;
for n = 1:Run_Num
tic
n
Function_name               = 'F21'
[lb,ub,dim,fobj]            = Get_Functions_details(Function_name);
costfunction                = fobj;     

nvar    = dim;
varsize = [1 nvar];
varmin  = lb;
varmax  = ub;

%% 
empty_ball.position     =[];
empty_ball.cost         =[];
empty_ball.velocity     =[];
empty_ball.acceleration =[];

ball             = repmat(empty_ball,npop,1);
globalbest.cost  = inf;

for i = 1:npop
  
    
    ball(i).position       = unifrnd(varmin,varmax,varsize);
    ball(i).velocity       = zeros(varsize);
    ball(i).Acceleration   = zeros(varsize);
    ball(i).sbetter        = zeros(varsize);
    ball(i).mean           = zeros(varsize);
    ball(i).cost           = costfunction(ball(i).position);
    if ball(i).cost < globalbest.cost
       globalbest.position = ball(i).position;
       globalbest.cost     = ball(i).cost;    
   end
end

bests = zeros(maxt,1);
T     = m_Ratio.*maxt;  
%% 
for t = 1:maxt
     
    sumcost = 0;
    s       = 1;
    for i= 1:npop
        ball(i).sbetter = ball(i).position;
            for j= 1:npop
            df = ball(j).cost - ball(i).cost;
            if df < 0
               ball(i).sbetter = ball(i).sbetter + ball(j).position;
               s               = s+1;
            end
        end
        ball(i).mean         = ((ball(i).sbetter) ./ s);
       
        
        P_MEAN = F.*(maxt./t);                % Eq. (15) in the paper
        k1     = (1./t)^(Beta) ;              % Eq. (16) in the paper
        k2     = c ./ (1 + exp( - (t-T)));    % Eq. (17) in the paper
        
        ball(i).velocity    = globalbest.position-ball(i).position;       % Eq. (14) in the paper
        ball(i).Acceleration = P_MEAN .* ball(i).mean - ball(i).position; % Eq. (17) in the paper
           
        ball(i).position = ball(i).position + ...
                          k1 .* (ball(i).Acceleration) .* rand(varsize)+...
                          k2 .* ball(i).velocity .* rand(varsize);  % Eq. (5) in the paper
        
        ball(i).position = min(max(ball(i).position,varmin),varmax);
        ball(i).cost     = costfunction(ball(i).position);
      
        if ball(i).cost < globalbest.cost
           globalbest.position  = ball(i).position;
           globalbest.cost      = ball(i).cost; 
        end
       bests(t) = globalbest.cost; 
       sumcost   = sumcost+ball(i).cost;
    end
    
disp(['Iteration' num2str(t) ':bestcost=' num2str(bests(t))]);
meanfits(t) = sumcost/npop;

end

Bests(n) = bests(t);
RunTime(n)=toc;
end

        disp([' ']);
        disp(['                   SIPO                         ']);
        disp(['-----------------------------------------------']);
        disp(['Number of run     = ' num2str(Run_Num)]);
        disp([' ']);
        disp([' ']);
        disp(['****************    Statistical indexes : Time    ****************']);
        disp(['------------------------------------------------']);
        disp(['Per run            = ' num2str(RunTime)]);
        disp(['Average            = ' num2str(mean(RunTime))]);
        disp(['Standard deviation = ' num2str(std(RunTime))]);
        disp(['Maximum            = ' num2str(max(RunTime))]);
        disp(['Minimum            = ' num2str(min(RunTime))]);
        
        disp([' ']);
        disp(['*****************   Statistical indexes : Fitness    ****************']);
        disp(['-----------------------------------------------']);
        disp(['Number of run      = ' num2str(Run_Num)]);
        disp(['Best cost per run  = ' num2str(Bests)]);
        disp(['Average            = ' num2str(mean(Bests))]);
        disp(['Standard deviation = ' num2str(std(Bests))]);
        disp(['Maximum            = ' num2str(max(Bests))]);
        disp(['Minimum            = ' num2str(min(Bests))]);
    




% semilogy(bests,'color','r','linewidth',2);
% hold on
% semilogy(meanfits,'-.','color','k','linewidth',2);
% % title('Version.2')
% xlabel('Iteration');
% ylabel('');
% axis tight
% grid on
% box on
% legend('Bests','Meanfits')


figure(1);
plot(bests,'color','r','linewidth',2);
% hold on
% plot(meanfits,'-.','color','k','linewidth',2);
% title('Version.2')
xlabel('Iteration');
ylabel('Fitness Value');
% axis tight
% grid on
% box on
legend('SIPO')
% title('Convergence Curve')
% figure(1);
% plot(P_MEAN,'k-','LineWidth',1.5)

