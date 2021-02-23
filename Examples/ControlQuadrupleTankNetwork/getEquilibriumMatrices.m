function [alpha,beta] = getEquilibriumMatrices()
    Cte = getConstantsWL();
    x = sym('x',[1 Cte.n],'real');
    assume(x,'positive');
    u = sym('u',[1,Cte.m],'real');
    assume(u,'positive');
    xdot = sym(zeros(Cte.n,1));
    for i = 1:Cte.n/2
       j = i+Cte.n/2;
       xdot(i) = -(Cte.a(i)/Cte.A(i))*sqrt(2*Cte.g)*x(i)+(Cte.a(j)/Cte.A(i))*sqrt(2*Cte.g*x(j))+...
           Cte.gamma(i)*Cte.k(i)*u(i)/Cte.A(i) == 0;
    end
    i = i+1;
    j = Cte.n/2;
    xdot(Cte.n/2+1) = -(Cte.a(i)/Cte.A(i))*sqrt(2*Cte.g*x(i))+...
           (1-Cte.gamma(j))*Cte.k(j)*u(j)/Cte.A(i) == 0;
    for i = Cte.n/2+2:Cte.n
       j = i-Cte.n/2-1; 
       xdot(i) = -(Cte.a(i)/Cte.A(i))*sqrt(2*Cte.g*x(i))+...
           (1-Cte.gamma(j))*Cte.k(j)*u(j)/Cte.A(i) == 0;
    end
    sol = solve(xdot,[x(Cte.n/2+1:end) u(1:end)]);
    sol = struct2cell(sol);
    alpha = zeros(Cte.n/2,Cte.n/2+nchoosek(Cte.n/2,2));
    for i = Cte.n/2+1:Cte.n
        [alpha(i-Cte.n/2,:),~] = coeffs(sol{i-Cte.n/2,1},x(1:Cte.n/2));
    end
    beta = zeros(Cte.m,Cte.n/2);
    for i = Cte.n+1:3*Cte.n/2
        [beta(i-Cte.n,:),~] = coeffs(sol{i-Cte.n/2,1},x(1:Cte.n/2));
    end
    save('EqMatrices.mat','alpha','beta');
end

