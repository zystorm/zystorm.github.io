clear
close all
clc
tic

%% 初始化
G = 200;
m = 35;
Q = 100; %信息素
Rho = 0.2;
Alpha = 1.2; % 1.4 2.2
Beta = 4.2; %启发式信息
addpath("Dataset-master/TSPLIB");
load eil51.txt;
CitysData = eil51;

Rhomax = 0.3;
Rhomin = 0.3;
TauMax = 50;
TauMin = 0.5;

CityNum = length(CitysData);
CitysData = CitysData(:,2:3);
DisTable = zeros(CityNum);
for ii = 1:CityNum
    for jj = ii:CityNum
        if ii == jj
            DisTable(ii,jj) = eps;
        else
            DisTable(ii,jj) = round(norm(CitysData(ii,:)-CitysData(jj,:)));
            DisTable(jj,ii) = DisTable(ii,jj);
        end
    end
end

%% 进化
% LbestAll = nan(round((4-1)/0.2)*round((5-0)/0.2),3);
% RbestAll = nan(round((4-1)/0.2)*round((5-0)/0.2),CityNum);

% for a = 1:round((4.2-1)/0.2)
% Alpha = 0.8+0.2*a;
% for b = 1:round((5.2-0)/0.2)
% Beta = -0.2+0.2*b;

LbA = zeros(1,100);
RbA = zeros(1,CityNum);
for dd = 1:100
%周游循环
Tau = zeros(CityNum); %信息素-迭代
Eta = 1./DisTable; %启发式-不变
Rbest = zeros(G,CityNum);
Lbest = zeros(G,1);
CityList = 1:1:CityNum;

for n = 1:G
    if n == 1
        Tau = 0.5 + 10./DisTable;
%         Rho = Rhomax;
    else
        Tau = (1-Rho)*Tau + DeltaTau;
%         if Rho >= Rhomin
%             Rho = Rho*0.95;
%         else
%             Rho = Rhomin;
%         end
        %%%%% 修改1 %%%%%%%
        if mod(n,5) == 0
            Tau(Tau>TauMax) = TauMax;
            Tau(Tau<TauMin) = TauMin;
        end
    end
    % Ant
    Tabu = zeros(m,CityNum);
    Lall = zeros(m,1);
    DeltaTau = zeros(CityNum);
   for k = 1:m
        CityNext = randi(CityNum);
        Tabu(k,1) = CityNext;
        for c = 2:CityNum
            CityTemp = CityNext;

            p1 = (Tau(CityTemp,:).^Alpha).*(Eta(CityTemp,:).^Beta);
            p1(Tabu(k,1:(c-1))) = 0;
            P = [CityList;p1/sum(p1)];
            P(:,Tabu(k,1:(c-1))) = [];
            
            % 轮盘赌法
            Pcum = P(2,:);
            Pcum = cumsum(Pcum);
            Select = find(Pcum>=rand);
            CityNext = P(1,Select(1));
            Tabu(k,c) = CityNext;

            Lall(k) = Lall(k)+DisTable(Tabu(k,c-1),Tabu(k,c));
        end
        Lall(k) = Lall(k)+DisTable(Tabu(k,end),Tabu(k,1));
%         Begin = Tabu(k,1:end)';
%         End = Tabu(k,2:end)';
%         End = [End;Begin(1)];
%         Index = (End-1)*CityNum + Begin;
%         DeltaTau(Index) = DeltaTau(Index) + Q/Lall(k);

        %%%%%%%% 修改3 %%%%%%%%%
        
   end
   LTabu = [Lall,Tabu];
   LTabu = sortrows(LTabu,1);
   Tabu1 = LTabu(:,2:end);
   Lall1 = LTabu(:,1);
   for k = 1:10
        Begin = Tabu1(k,1:end)';
        End = Tabu1(k,2:end)';
        End = [End;Begin(1)];
        Index = (End-1)*CityNum + Begin;
        DeltaTau(Index) = DeltaTau(Index) + (11-k)*Q/Lall1(k);
   end
    [Lbest(n),IndexMin] = min(Lall);
    Rbest(n,:) = Tabu(IndexMin,:);
    %%%%%% 修改2 %%%%%%%
%     if Lbest(n) < min(Lbest(1:m))+100
%         Begin = Rbest(n,1:end)';
%         End = Rbest(n,2:end)';
%         End = [End;Begin(1)];
%         Index = (End-1)*CityNum + Begin;
%         DeltaTau(Index) = DeltaTau(Index) + 10*Q/Lbest(n);
%     elseif n>1 && Lbest(n) > Lbest(n-1)+20 
%         Begin = Rbest(n,1:end)';
%         End = Rbest(n,2:end)';
%         End = [End;Begin(1)];
%         Index = (End-1)*CityNum + Begin;
%         DeltaTau(Index) = DeltaTau(Index)*0.8;
%     end
end
[LbA(dd),index] = min(Lbest);
RbA = Rbest(index,:);
LbA(dd)
end
mean(LbA)
min(LbA)
% LbestAll((a-1)*25+b,:) = [Alpha,Beta,Lbest(end)];
% end
% end
toc

%% plot
figure(1)
[L,index] = min(Lbest);
R = Rbest(index,:);
scatter(CitysData(R,1),CitysData(R,2));
hold on
plot(CitysData(R,1),CitysData(R,2));
xlabel('x');
ylabel('y');

figure(2)
plot(1:1:G,Lbest);
xlabel('代');
ylabel('最小值');

% LbestAll = sortrows(LbestAll,3);
% c = linspace(1,10,401);
% s = scatter3(LbestAll(:,1),LbestAll(:,2),LbestAll(:,3),15,c,'filled');
% xlabel('\alpha');
% ylabel('\beta');
% zlabel('最短路径');