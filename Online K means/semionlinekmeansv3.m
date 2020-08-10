close all
clear all

n = 100000;

load('4QAM_SNR_5.mat');
I = real(receivedSignal);
Q = imag(receivedSignal);


msize = numel(I);
idx = randperm(msize);


%scatter(I(idx(1:1000)),Q(idx(1:1000)),'.')

X = [I(idx(1:n)),Q(idx(1:n))];
rng(1); % For reproducibility
[index,C,~,d] = kmeans(X,4,'Replicates',4);

w = sum((min(d,[],2)).^2);


%w = 78.9055;  %15 SNR 4qam
%w = 795.5980;  %4 SNR 4qam

figure;
plot(X(index==1,1),X(index==1,2),'r.','MarkerSize',12)
hold on
plot(X(index==2,1),X(index==2,2),'b.','MarkerSize',12)
hold on
plot(X(index==3,1),X(index==3,2),'y.','MarkerSize',12)
hold on
plot(X(index==4,1),X(index==4,2),'m.','MarkerSize',12)
hold on
plot(X(index==5,1),X(index==5,2),'g.','MarkerSize',12)
hold on
plot(X(index==6,1),X(index==6,2),'c.','MarkerSize',12)
hold on
plot(X(index==7,1),X(index==7,2),'w.','MarkerSize',12)
hold on
plot(X(index==8,1),X(index==8,2),'k.','MarkerSize',12)

plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Centroids',...
       'Location','NW')
title 'Cluster Assignments and Centroids'
hold off


% load('2PSK_SNR_10.mat');
% I = real(receivedSignal);
% Q = imag(receivedSignal);
X = [I,Q];


k=4;
f=w/(log10(n));
clusters = zeros(n,2);
F = [];
Q = [];
R = [];
indexofclusterpoint = [];
indexoflastpoint = [];
    C = [];
    q=0;r=1;
    clusters = zeros(n,2);
for i=1:n
   x = X(i,:);
   if isempty(C)==1
       p = 1;   
   else
       D = zeros(size(C,1),1);
       for j = 1:size(C,1)
        D(j) = sum((x-C(j,:)).^2);
       end
   [~,I] = min(D);
   p = min(sum((x-C(I,:)).^2)/f,1);
   c = C(I,:);
   end
   if binornd(1,p) >= 1
       C = [C; x];
       q = q+1;
       indexofclusterpoint = [indexofclusterpoint; i];
   end
   
      %-----------------------------
   Ddash = zeros(size(C,1),1);
       for in = 1:size(C,1)
        Ddash(in) = sum((x-C(in,:)).^2);
       end
   [~,J] = min(Ddash);
   clusters(i,:) = C(J,:);
   
   
   %-----------------------------
   Q = [Q; q];
   if q >= 0.3*k*(1+log10(n))
       
       r = r+1;q = 0; f = 2*f;
       F = [F; f];
       R = [R; r];
       C = [];
       indexoflastpoint = [indexoflastpoint; i];
   end

end


  figure;
for i=indexoflastpoint(end-1):n
    if clusters(i,:) == C(1,:)
        plot(X(i,1),X(i,2),'r.','MarkerSize',12);
    elseif clusters(i,:) == C(2,:)
        plot(X(i,1),X(i,2),'b.','MarkerSize',12);
    elseif clusters(i,:) == C(3,:)
        plot(X(i,1),X(i,2),'g.','MarkerSize',12);
    
    elseif clusters(i,:) == C(4,:)
        plot(X(i,1),X(i,2),'y.','MarkerSize',12);
    elseif clusters(i,:) == C(5,:)
        plot(X(i,1),X(i,2),'m.','MarkerSize',12);
    %{
        elseif clusters(i,:) == C(6,:)
        plot(X(i,1),X(i,2),'c.','MarkerSize',12);
    
    elseif clusters(i,:) == C(7,:)
        plot(X(i,1),X(i,2),'w.','MarkerSize',12);
    elseif clusters(i,:) == C(8,:)
        plot(X(i,1),X(i,2),'k.','MarkerSize',12);   
     %}
    end
    hold on
    
end
 plot(C(:,1),C(:,2),'kx',...
      'MarkerSize',15,'LineWidth',3) 
 legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Centroids',...
        'Location','NW')
title 'Cluster Assignments and Centroids'
hold off
   

[iC,ia,ic] = unique(clusters);
a_counts = accumarray(ic,1);
value_counts = [iC, a_counts]; 