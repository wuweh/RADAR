%% 多项式重采样
function outIndex = multinomialR(inIndex,q)
 
  error('下面的参数nargin请参考书中的值设置，然后删除本行代码') 
if nargin < 0, error('Not enough input arguments.'); end

[S,arb] = size(q);  

 
N_babies= zeros(1,S);
cumDist= cumsum(q');   
 
u = fliplr(cumprod(rand(1,S).^(1./(S:-1:1))));
j=1;
for i=1:S
  while (u(1,i)>cumDist(1,j))
    j=j+1;
  end
  N_babies(1,j)=N_babies(1,j)+1;
end
 
index=1;
for i=1:S
  if (N_babies(1,i)>0)
    for j=index:index+N_babies(1,i)-1
      outIndex(j) = inIndex(i);
    end
  end  
  index= index+N_babies(1,i);   
end