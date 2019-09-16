function Ksir_dir = KSIR(K, y, NumOfSlice)

B = 1;

[n p] = size(K);
[Sorty Index] = sort(y);
K = K(Index,:);
Kmean=mean(K,1);

% extract centered and weigthed slice means
    % for regression 
SizeOfSlice = fix(n/NumOfSlice); % size of each slice
m = mod(n,NumOfSlice);
base = zeros(2,1);
smean_c=zeros(NumOfSlice,p);
for k = 1:NumOfSlice
    count = SizeOfSlice+(k<m+1);
    base(2) = base(2) + count;
    smean_c(k,:) =(mean(K(base(1)+1:base(2),:),1)-Kmean)*sqrt((base(2)-base(1))/n); 
    % k-th slice mean, centered
    base(1) = base(2);
end

% solve the following generalized eigenvalue problem
% Cov(HK)*V = lamda*Cov(K)*V
Cov_K=K'*K/n-Kmean'*Kmean; 
clear K Kmean
Temp = (Cov_K+10^(-10)*eye(p))\smean_c'; % compute inv(Cov_K)W
clear Cov_K

[U D] = eig(smean_c*Temp); % extract U via solving W'inv(Cov_K)WU=UD
[D Index] = sort(diag(D),'descend');
D = D(1:end-1);
U = U(:,Index(1:end-1));
Ksir_dir = B*(Temp*U)*diag(1./sqrt(D)); % normalization 
%Ksir_dir =orth(Ksir_dir);





