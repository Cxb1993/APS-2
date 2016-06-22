
clear

N= 64;
n = N+1; k=8; L=1;
u = zeros(N,1);
x = linspace(0,1,n);
u = sin(k*pi*x/L);
plot(x,u,'blue','LineWidth',3);  hold on;


  set(gca,'fontsize',30)
  set(gca,'fontweight','bold')

for iter=1:10
    
    for i=2:N

        u(i) = (u(i+1)+u(i-1))/2;

    end
    
end  
plot(x,u,'red','LineWidth',3);  hold on;
set(gca,'fontsize',30)
 set(gca,'fontweight','bold')
 
xlabel('x/L','FontSize',14,'FontWeight','bold')
ylabel('u','FontSize',14,'FontWeight','bold')

legend('Initial','Final',14,'FontWeight','bold');

% lambda vs k
N=64;
k = linspace(1,60,100);

for j=1:length(k)
lamda(j) = 1 - (sin(k(j)*pi*0.5/N)).^2;
end

plot(k,lamda,'blue','LineWidth',3);
xlabel('wavenumber, k','FontSize',14,'FontWeight','bold')
ylabel('eigenvalue, \lambda','FontSize',14,'FontWeight','bold')

 set(gca,'fontsize',30)
 set(gca,'fontweight','bold')
 
 % K/N ratio
 
k=1;
N = linspace(2,20,65);
temp= 1./N;
lamda = 1 - (sin(k*pi*0.5*temp)).^2;

plot(temp,lamda,'blue','LineWidth',3);
xlabel('ratio, k/N','FontSize',14,'FontWeight','bold')
ylabel('eigenvalue, \lambda','FontSize',14,'FontWeight','bold')
axis([0 0.5 0.5 1]);
 set(gca,'fontsize',30)
 set(gca,'fontweight','bold')
% N vs lambda
 k=1;
N = linspace(2,20,65);
temp= 1./N;
lamda = 1 - (sin(k*pi*0.5*temp)).^2;

plot(N,lamda,'blue','LineWidth',3);
xlabel('ratio, N','FontSize',14,'FontWeight','bold')
ylabel('eigenvalue, \lambda','FontSize',14,'FontWeight','bold')
%axis([0 0.5 0.5 1]);
 set(gca,'fontsize',30)
 set(gca,'fontweight','bold')


hold off;