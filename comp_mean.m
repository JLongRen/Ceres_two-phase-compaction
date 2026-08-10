function [Kd] = comp_mean(K,p,N) 

        K1=K(1:N-1);
        K2=K(2:N);
        mean = zeros(N+1,1); 
        mean(2:N) = (0.5*(K1.^p(2:N)+K2.^p(2:N))).^(1./p(2:N));
        mean(1)=K(1); 
        mean(N+1)=K(N);
        Kd = sparse(diag (mean));

end