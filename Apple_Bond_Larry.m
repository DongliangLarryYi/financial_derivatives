m = 101; % total number for stock split
n = 21; % total number of interest rate split
upper_r = 0.1;
upper_stock = 1000;
ss = 0.07;
ks = 0.09; %k_*
rs =0.05; % r_*
w = 0.3*sqrt(430);
dB = 0.024; % coupon rate
dS = 0.023; % divident rate, which is different from DS!!

%Step1

% initial the final state
F(m,n)=0; % m is the number of row, which means different stock price.
r=linspace(0,upper_r,n);
S=linspace(0,upper_stock,m);
F=F+50; % initilize the bondary when default
for i =2:m
    for j=1:n
        F(i,j)=100;
    end
end % initialize the final bond value (percentatge of face value)

Fold=F;
Dr=upper_r/(n-1);
Dt=1/52;
DS=upper_stock/(m-1);
N=floor(10/Dt); %total number of time split

%initialize the transformation matrix for interest rate
A(n,n)=0;
A(1,1)=ss^2*r(1)/(2*Dr^2)-ks*(rs-r(1))/Dr;
A(1,2)=-ss^2*r(1)/(Dr^2)+ks*(rs-r(1))/Dr;
A(1,3)=ss^2*r(1)/(2*Dr^2);
A(n,n-2)=ss^2*r(n)/(2*Dr^2);
A(n,n-1)=-ss^2*r(n)/(Dr^2)-ks*(rs-r(n))/Dr;
A(n,n)=ss^2*r(n)/(2*Dr^2)+ks*(rs-r(n))/Dr;
for i =2:(n-1)
    A(i,i-1)=ss^2*r(i)/(2*Dr^2)-ks*(rs-r(i))/(2*Dr);
    A(i,i)=-ss^2*r(i)/(Dr^2);
    A(i,i+1)=ss^2*r(i)/(2*Dr^2)+ks*(rs-r(i))/(2*Dr);
end
%%%%%%%%%

Fnew=F;

for h=1:N
    W=zeros(n,m);%hold stock price first
    for i=1:m
        W(:,i)=A*(Fold(i,:)')*Dt+Fold(i,:)'; 
    end
    % based calculated W, then update backward
    for j =1:n
        rj=r(j);
        
        % update transformation matrix under different interest rate
        B(m-1,m-1)=0;
        B(1,1)=1-Dt*( -S(2)*w^2/(DS^2)-(rj-dB) );
        B(1,2)=-Dt*( S(2)*w^2/(2*DS^2)+(rj-dS)*S(2)/(2*DS) );
        B(m-1,m-3)=-Dt*(S(m)*w^2/(2*DS^2));
        B(m-1,m-2)=-Dt*( -S(m)*w^2/(DS^2)-(rj-dS)*S(m)/DS );
        B(m-1,m-1)=1-Dt*( S(m)*w^2/(2*DS^2)+(rj-dS)*S(m)/DS-(rj-dB) );
        for i =2:(m-2)
            B(i,i-1)=-Dt*( S(i+1)*w^2/(2*DS^2)-(rj-dS)*S(i+1)/(2*DS) );
            B(i,i)=1-Dt*( -S(i+1)*w^2/(DS^2)-(rj-dB) );
            B(i,i+1)=-Dt*( S(i+1)*w^2/(2*DS^2)+(rj-dS)*S(i+1)/(2*DS) );
        end
        
        C=W(j,2:m);
        
        C(1)=C(1)+Dt*( S(2)*w^2*25/(DS^2)-(rj-dS)*S(2)*25/DS );
        
        Fnew(2:m,j)=inv(B)*(C'); %update along time
    end
    Fold=Fnew;
    disp(sprintf('wait...'))
    disp(h/N)
end

%find the price of bond, linear interpolation
Price = Fold(430/DS,1)*0.8 + Fold(430/DS,2)*0.2

%plot the price surface at t0 for STEP1
surf(r,S(2:m),Fold(2:m, :))



%step2


CC(N,n)=0;
T=10;
for i=1:N
    for j=1:n
        
        % the CIR formula with respect to different r and T-t
        
        t = (i)*(1/52);
        DT=T-t;
        
        
        h=sqrt(ks^2+2*ss^2);
        H1=( 2*h*exp((ks+h)*DT/2)/(2*h+(ks+h)*(exp(DT*h)-1)) )^(2*ks*rs/(ss^2));
        H2=2*(exp(DT*h)-1)/( 2*h+(ks+h)*(exp(DT*h)-1) );
        Bond=100*H1*exp(-H2*r(j));
        if DT==0
            Bond=100;
        end
        
        %add coupon to the call price
        spread=log(100/Bond)/(DT)+0.0015;
        I=2.4*(1-exp(-spread*(DT)))/spread;
        if DT==T
            I=0;
        end
        CC(i,j)=max(Bond+I,100);
    end
end


figure;
%plot the call price for STEP2
surf(r,linspace(0,10,N),CC,gradient(CC))

%Step3 stock price is not related to the call price, we just need to
%compare the result form step1 and step2 to get the updated price of this
%bond


F_c(m,n)=0;

F_c=F_c+50;
for i =2:m
    for j=1:n
        F_c(i,j)=100;
    end
end
Fold_c=F_c;	


%initialize the transformation matrix for interest rate
A(n,n)=0;
A(1,1)=ss^2*r(1)/(2*Dr^2)-ks*(rs-r(1))/Dr;
A(1,2)=-ss^2*r(1)/(Dr^2)+ks*(rs-r(1))/Dr;
A(1,3)=ss^2*r(1)/(2*Dr^2);
A(n,n-2)=ss^2*r(n)/(2*Dr^2);
A(n,n-1)=-ss^2*r(n)/(Dr^2)-ks*(rs-r(n))/Dr;
A(n,n)=ss^2*r(n)/(2*Dr^2)+ks*(rs-r(n))/Dr;
for i =2:(n-1)
    A(i,i-1)=ss^2*r(i)/(2*Dr^2)-ks*(rs-r(i))/(2*Dr);
    A(i,i)=-ss^2*r(i)/(Dr^2);
    A(i,i+1)=ss^2*r(i)/(2*Dr^2)+ks*(rs-r(i))/(2*Dr);
end
%%%%%%%%%


Fnew_c=F_c;

for h=1:N
    W=zeros(n,m);
    for i=1:m
        W(:,i)=A*(Fold_c(i,:)')*Dt+Fold_c(i,:)';
    end
    for j =1:n
        rj=r(j);
        
        
        
        % update transformation matrix under different interest rate
        B(m-1,m-1)=0;
        B(1,1)=1-Dt*( -S(2)*w^2/(DS^2)-(rj-dB) );
        B(1,2)=-Dt*( S(2)*w^2/(2*DS^2)+(rj-dS)*S(2)/(2*DS) );
        B(m-1,m-3)=-Dt*(S(m)*w^2/(2*DS^2));
        B(m-1,m-2)=-Dt*( -S(m)*w^2/(DS^2)-(rj-dS)*S(m)/DS );
        B(m-1,m-1)=1-Dt*( S(m)*w^2/(2*DS^2)+(rj-dS)*S(m)/DS-(rj-dB) );
        for i =2:(m-2)
            B(i,i-1)=-Dt*( S(i+1)*w^2/(2*DS^2)-(rj-dS)*S(i+1)/(2*DS) );
            B(i,i)=1-Dt*( -S(i+1)*w^2/(DS^2)-(rj-dB) );
            B(i,i+1)=-Dt*( S(i+1)*w^2/(2*DS^2)+(rj-dS)*S(i+1)/(2*DS) );
        end
        
        
        
        C=W(j,2:m);
        C(1)=C(1)+Dt*( S(2)*w^2*25/(DS^2)-(rj-dS)*S(2)*25/DS );
        %%%F1=B\(C');
        Fnew_c(2:m,j)=inv(B)*(C');
    end
    Ct=CC(N+1-h,:);
    for i=1:m
        Fnew_c(i,:)=min(Fnew_c(i,:),Ct);
    end
    Fold_c=Fnew_c;
    disp(sprintf('wait...'))
    disp(h/N)
end
F_c=Fnew_c;

Price_c = Fold_c(430/DS,1)*0.8 + Fold_c(430/DS,2)*0.2



figure;
%plot the price for STEP3
surf(r,S(2:m),Fold(2:m, :))

