function [x,lambda] = SubP_alpha( y,s,r,alpha )
% This code solves the subproblem of SCP_ls, which can be reformulated as 
% min ||x||_1 + alpha/2*||x - y||^2
% s.t. ||x - s||^2<= r
% for some y, s, alpha, r.

   [n,~] = size(y);
    
    
    % test for interior solution
    tmp = sign(y).*max(abs(y)-1/alpha,0);
    if norm(tmp - s)^2 < r
        x = tmp;
        lambda = 0;
%         fprintf('interior solution\n')
    else        
        %%compute the dividing points
        v = alpha*(y-s);
        PC = -s./(v + ones(n,1));
        PD = -s./(v - ones(n,1));
        
        p1 = sign(v);
        p2 = v + ones(n,1);
        p3 = v - ones(n,1);
        [Case,sg] = CaseSg(p1,p2,p3,PC,PD,v,s,n);%computing each cases and segments of each case
        
      
        
        p = sg(:);
        [pp,I] = sort(p,'ascend');
        [N,~] = size(I);
        count = ones(n,1);
        Fs1 = -0.1*ones(N,1);
        beta = -0.1;
        a = sum(Case(:,1));
        c = sum(Case(:,2));
        ae =  sum(Case(:,5));
        ce =  sum(Case(:,6));
        Fs1(N) = ae*(pp(N))^2 + ce - r;
        if Fs1(N)<= 0
            beta = sqrt((r-ce)/ae);
        else
        for j=1:N
            % computing function values on sigments
            Fs1(j) = a*(pp(j))^2 + c - r;
            if Fs1(j)>=0
                beta = sqrt((r-c)/a);
               % fprintf('beta %6.5e\n', beta)

                break
            end
            k = mod(I(j),n);
            if k == 0
                k=n;
            end
            
            if count(k,1) <= 2
                a = a - Case(k,2*count(k) - 1) + Case(k,2*count(k) + 1);
                c = c - Case(k,2*count(k) ) + Case(k,2*count(k) + 2);
                count(k,1) = count(k,1) +1;
             end
        end
        end
        %plot(1:N,Fs1,'ro');
        
        
        
        %plot(Fs1,'r*');
        
        z = beta*v + s;
        x = max(0,abs(z)-beta).*sign(z);
        lambda = (1/beta - alpha)/2;
        
    end
end

