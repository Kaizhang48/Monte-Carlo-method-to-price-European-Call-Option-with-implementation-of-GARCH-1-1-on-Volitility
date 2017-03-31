    times=30000000;
    year=252;
    epr=0.5;
    T=epr*year;
    step=T;
    dt=epr/step;
    alpha = 0.080; beta = 0.885; gamma = 0.035;
    ryear=0.05;r=ryear*dt;
    S0=100;
    sigma_2=zeros(1,T);
    K=(50:10:150)';
    [K_num,~]=size(K);
    op=zeros(K_num,times);
    V_2=power(0.3,2)*dt;
    d1=gamma*V_2;
    sigma_2(1)=power(0.35,2)*dt;
    for x=1:times
        N=normrnd(0,1,1,T);
        for i=2:T
            sigma_2(i)=d1+(beta+alpha*power(N(i-1),2))*sigma_2(i-1);
        end
        drift=(r-1/2.*sigma_2);
        sigma=power(sigma_2,1/2);
        rand=sigma.*N;
        lnds=cumsum(drift+rand);
        sp=S0*exp(lnds);
        op(:,x)=max([sp(end)-K,zeros(K_num,1)],[],2);
    end
    result=sum(op,2)./times;
    vol=zeros(1,K_num);
    for i=1:K_num
         vol(i)=european_formula_volatility(1,result(i,1),K(i,1),epr,S0,0,ryear,0.1,0.5,0.0000001,10000);
    end

    plot(K',vol,'.','markersize',20)
    title('HOMEWORK1 #3')
    xlabel('Strike Price')
    ylabel('Vol')
