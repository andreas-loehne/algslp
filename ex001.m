clear all

1;

function [sol_cone,pi_matrix] = create_solvency_cone(pi_matrix)
	n=size(pi_matrix,1);
	gen_matrix = zeros(n, n^2);
	for i=1:n
		for j=1:n
			if i!=j
				gen_matrix(i,n*(j-1)+i)=pi_matrix(i,j);
				gen_matrix(j,n*(j-1)+i)=-1;
			endif
		endfor
	endfor
	sol_cone=polyh(struct('V',zeros(n,1),'D',gen_matrix),'v');
endfunction

n=2;
d=4;

pi_matrix_1=[
      1    223/100      89/50      47/25;
157/100          1      54/25    177/100;
  19/10    211/100          1    199/100;
 101/50      41/20      43/20          1];
 
pi_matrix_2=[
    1      39/20    223/100    207/100;
37/20          1      41/20      54/25;
38/25      97/50          1      37/20;
 11/5      49/25      44/25          1];

 % for i=1:d
 % 	 for j=1:d
 % 		 for k=1:d
 % 			 pi_matrix_1(i,j)<=pi_matrix_1(i,k)*pi_matrix_1(k,j)
 % 		 endfor
 % 	 endfor
 % endfor
 % for i=1:d
 % 	 for j=1:d
 % 		 for k=1:d
 % 			 pi_matrix_2(i,j)<=pi_matrix_2(i,k)*pi_matrix_2(k,j)
 % 		 endfor
 % 	 endfor
 % endfor
 
 
[K1,pi_matrix_1]=create_solvency_cone(pi_matrix_1);
[K2,pi_matrix_2]=create_solvency_cone(pi_matrix_2); 

% set-valued risk measure
x_init=[zeros(n,1);-ones(d-n,1)];
A=K2; % acceptance set
M=[eye(d),[eye(n);zeros(d-n,n)]];
sop.gr=A.inv(M)&((point(x_init)-K1):space(n));
sop.cone=cone(2);

diary('pics/output1.txt');
tic
res=so_solve(sop);
toc
opt.axes=[-3 7 -3 7];
opt.dirlength=3;
opt.name='ex001_';
opt.path='pics/ex001';
show(res,opt);
diary off







