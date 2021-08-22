% Input
% x=[12 5 34 87 35 90 124];
% N=[2 4 56 12 6 9 14];

function ave=myaverage(x,N)

  sizex=size(x);
  sizeN=size(N);
  %if(sizex(2)>sizeN(2)>sizeN(2)|sizex(2)<sizeN(2)
  if sizex(2) ~=sizeN(2)
      beep
      disp('Error the number of elements inside input data should be the same')
  else
      total=sum(N);
      s=x.*N;
      ave=sum(s)/total;
end
end
