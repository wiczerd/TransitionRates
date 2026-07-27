function stat=Ergodic(Tran)
N=size(Tran,1);
A1=Tran-eye(N);
AA=[A1(:,1:N-1),ones(N,1)];
AAinv=AA^(-1);
stat=AAinv(N,:)';  %  stationary distribution across the states