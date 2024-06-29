function S = so_solve(sop)
	q=sop.cone.sdim;
	n=sop.gr.sdim-q;
	S.poi=zeros(n,0);
	S.dir=zeros(n,0);
	S.init.poi=zeros(n,0);
	S.init.dir=zeros(n,0);
	
	S.gr_F_std = sop.gr+(origin(n):sop.cone);
	S.gr_G_std = S.gr_F_std.recc;
	S.C_std    = so_val(S.gr_G_std,zeros(n,1)).eval;
	
	S.PP=(S.gr_F_std.im([zeros(q,n),eye(q)])).eval;		
	if S.PP.isempty
		disp('Set optimization problem is not feasible.');
		return;
	endif
	
	if isempty(so_compute_minimizer(S.gr_G_std,zeros(q,1)))
		disp('Set optimization problem has no solution (unbounded).');
		return;
	endif
	
	P=S.PP.vrep.V;
	D=so_set_diff(S.PP.vrep.D,S.C_std.vrep.D); % directions of PP not belonging to C_std
	
	for i=1:size(P,2)
		[x,x_init]=so_compute_minimizer(S.gr_F_std,P(:,i));
		S.poi(:,end+1)=x;
		S.init.poi(:,end+1)=x_init;
	endfor
	
	for j=1:size(D,2)
		[x,x_init]=so_compute_minimizer(S.gr_G_std,D(:,j));
		S.dir(:,end+1)=x;
		S.init.dir(:,end+1)=x_init;
	endfor
endfunction