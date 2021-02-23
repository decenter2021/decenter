function [K,uEqDisc] = iLQR(x0,ref,Cte,d,m)

    % Compute number of linearizations
    system= cell(Cte.T(m)+1,7);
    
    % Forward pass variables
    x = cell(1,Cte.T(m)+1);
    uEqDisc = cell(1,Cte.T(m)+1);
    uControlDisc = cell(1,Cte.T(m));
    prevuControlDisc = cell(1,Cte.T(m));
    
    t_disc = 0:Cte.dT:Cte.T(m);
    
    for k = 1:Cte.iLQRIt
       % Forward pass
       if k == 1
           % Initial propagated system
           system(1,:) = getDiscreteDynamicsWL(x0,Cte);
           for i = 2:Cte.T(m)+1
               system(i,:) = system(1,:);
           end
       else
           for i = 1:Cte.T(m)
                if i == 1     
                     % Compute the linearized dynamics for the first instant 
                     system(i,:) = getDiscreteDynamicsWL(x0,Cte);
                     for l = i+1:i+Cte.dTlin-1
                         system(l,:) = system(i,:);
                     end
                     % output
                     x{1,i} = x0;
                     
                     uEqDisc{1,i} = ref{2,min(i,size(ref,2))}+(Cte.h*system{i,2}(1:Cte.n,1:Cte.m))\Cte.h*(ref{1,min(i+1,size(ref,2))}(1:Cte.n,1)-ref{1,min(i,size(ref,2))}(1:Cte.n,1));
                     uControlDisc{1,i} = -transpose(K{i,1})*(x{1,i}-ref{1,min(i,size(ref,2))})+uEqDisc{1,i};
                     uControlDisc{1,i}(uControlDisc{1,i}<0) = 0;
                     uControlDisc{1,i}(uControlDisc{1,i}>Cte.uMax) = Cte.uMax;
                     % Simulate nonlinear dynamics
                     nonLinSol = ode45(@(t,x) xdotContinuous(x,uControlDisc{1,min(floor(t/Cte.dT)+1,round(t_disc(i+1)/Cte.dT))},Cte),[t_disc(i) t_disc(i+1)],x0(1:Cte.n,1));
                     x{1,i+1}(1:Cte.n,1) = deval(nonLinSol, t_disc(i+1));
                     x{1,i+1}(Cte.n+1:3*Cte.n/2) = x{1,i}(Cte.n+1:3*Cte.n/2)+x{1,i+1}(1:Cte.n/2,1)-ref{1,i+1}(1:Cte.n/2,1);
                else
                    
                     if rem(i-1,Cte.dTlin) == 0 
                        system(i,:) = getDiscreteDynamicsWL(x{1,i},Cte);
                        for l = i+1:Cte.T(m)+1
                            system(l,:) = system(i,:);
                        end
                     end
                     uEqDisc{1,i} = ref{2,min(i,size(ref,2))}+(Cte.h*system{i,2}(1:Cte.n,1:Cte.m))\Cte.h*(ref{1,min(i+1,size(ref,2))}(1:Cte.n,1)-ref{1,min(i,size(ref,2))}(1:Cte.n,1));
                     uControlDisc{1,i} = -transpose(K{i,1})*(x{1,i}-ref{1,min(i,size(ref,2))})+uEqDisc{1,i};
                     uControlDisc{1,i}(uControlDisc{1,i}<0) = 0;
                     uControlDisc{1,i}(uControlDisc{1,i}>Cte.uMax) = Cte.uMax;
                     % Simulate nonlinear dynamics
                     nonLinSol = ode45(@(t,x) xdotContinuous(x,uControlDisc{1,min(floor(t/Cte.dT)+1,round(t_disc(i+1)/Cte.dT))},Cte),[t_disc(i) t_disc(i+1)],x{1,i}(1:Cte.n,1));
                     x{1,i+1}(1:Cte.n,1) = deval(nonLinSol, t_disc(i+1));
                     x{1,i+1}(Cte.n+1:3*Cte.n/2) = x{1,i}(Cte.n+1:3*Cte.n/2)+x{1,i+1}(1:Cte.n/2,1)-ref{1,min(i+1,size(ref,2))}(1:Cte.n/2,1);
                     %x{1,i+1} = system{i,1}*(x{1,i}-system{i,5})+system{i,2}*(uControlDisc{1,i}-system{i,6}) + system{i,5} - [zeros(Cte.n,Cte.n/2); eye(Cte.n/2)]*ref{1,min(i+1,size(ref,2))}(1:2,1);
                     
                     % Anti windup
                     for j = 1:Cte.n/2
                        if abs(x{1,i}(Cte.n+j)) > Cte.AntiWU(m)
                            x{1,i}(Cte.n+j) = Cte.AntiWU(m)*abs(x{1,i}(Cte.n+j))/x{1,i}(Cte.n+j); 
                        end
                     end
                end
            end
           
       end
       
       % Stopping criterion
       dif = zeros(1,Cte.T(m));
       if k > 2
           for i = 1:Cte.T(m)
                dif(1,i) = norm(prevuControlDisc{1,i}-uControlDisc{1,i})/norm(uControlDisc{1,i});
           end
           
           if max(dif) < Cte.iLQReps 
               break; 
           end
       end
       prevuControlDisc = uControlDisc;
       
       if m<= 2
           % One-step
            [~,K] = CentralizedLQR(system,Cte.T(m));
       else
           % One-step
            [~,K] = OneStepLQR(system,Cte.E,Cte.T(m));
       end
       
       if k == Cte.iLQRIt
           max(dif)
           fprintf(sprintf('iLQR max iterartions %d \n',m));
       end
       

    end

    K = transpose(K(1:d,1));
    uEqDisc = uEqDisc(1,1:d);
end
