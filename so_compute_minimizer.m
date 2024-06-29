function [x,x_init] = so_compute_minimizer(gr,y) 

	% set-valued map F is given by its graph gr	
	% computes a minimizer x of F with y in F(x)

	q=size(y,1);
	n=gr.sdim-q;

	% compute x such that y in F(x)
	x=(so_val_inv(gr,y)).getpoint();
	if ~isempty(x)
		x_init=x;
		Fx=(so_val(gr,x)).eval; % compute the value F(x)
			
		if Fx.dim<Fx.sdim % this can happen if the order cone has empty interior
			aff_Fx=Fx.aff.eval; % affine hull of F(x)
			cnt=0;
			do
				leave=true;
				NN=aff_Fx.hrep.Beq';
				if ~isempty(NN)				
					NN=[NN,-NN*ones(size(NN,2),1)];
				end
				N=so__normalize_cols(NN); % minimal system of outer normals of aff(F(x))
				for i=1:size(N,2)
					disp(sprintf('* expansion %i',cnt));
					cnt=cnt+1;
					w=N(:,i);
					res=so__lp(gr,Fx,w); % solve LP
					if strcmp(res.status,'opt')
						x=res.x;
						y=res.y;
					else
						x=zeros(n,0);
						return;
					endif
					if ~(point(y)<=aff_Fx) % if y not in aff(F(x))
						Fx=(so_val(gr,x)).eval; % update F(x) 
						aff_Fx=Fx.aff.eval;    % update aff(F(x)) 
						leave=false;
						break; 
					endif					
				endfor
			until leave	
			Fx=(Fx+((aff_Fx.lin)')).eval;
			gr=gr+(origin(n):(aff_Fx.lin)');
		endif
		
		N=so__normalize_cols(Fx.hrep.B'); % minimal system of outer normals of F(x)
		W=zeros(q,0); % initialies as empty set
		cnt=0;
		while true
			disp(sprintf('expansion %i',cnt));
			cnt=cnt+1;
			tmp=so_set_diff(N,W);
			if isempty(tmp) % N subseteq W
				break;
			endif
			w=tmp(:,1); % choose w in N \ W
			W=[W,w];
			res=so__lp(gr,Fx,w); % solve LP
			if strcmp(res.status,'opt')
				x=res.x;
				y=res.y;
			else
				x=zeros(n,0);
				break;
			endif
			if ~(point(y)<=Fx) % if y not in F(x)
				Fx=(so_val(gr,x)).eval; % update F(x)  
				N=so__normalize_cols(Fx.hrep.B');  % minimal system of outer normals of F(x)
			endif
			
		endwhile
		disp('');
	else	
		x=zeros(n,0); % there is no x with y in F(x)
		x_init=x;
	endif
endfunction

function res = so__lp(gr, Fx, w)
	q=Fx.sdim;
	n=gr.sdim-q;
	V=Fx.vrep.V;
	ii=size(V,2);
	GI=cell(1,ii);
	for i=1:ii
		v=V(:,i);
		GI{i}=so_val_inv(gr,v);
	endfor
	feas= gr&(intsec(GI):space(q));
	[optval,sol]=lpsolve([zeros(n,1);-w],feas);
	if optval>-Inf && ~isempty(sol)
		res.x=sol(1:n,1);
		res.y=sol(n+1:n+q,1);
		res.status='opt';
	elseif optval==-Inf
		res.x=zeros(1:n,0);
		res.y=zeros(n+1:n+q,0);
		res.status='unbounded';
	else	
		error('unexpected lp status');
	endif	
endfunction

function res = so_val_inv(gr,y)
	q=size(y,1);
	n=gr.sdim-q;
	res = (gr & (space(n):point(y))).im([eye(n),zeros(n,q)]);
endfunction

function matrix_out = so__normalize_cols(matrix_in,eps=1e-06) 
	if size(matrix_in,2)>0
		no=vecnorm(matrix_in);
		if min(no)>eps
		    matrix_out=matrix_in./(ones(size(matrix_in,1),1)*no);
		else
		    error('matrix contains all zero column');
		endif
	else
		matrix_out = matrix_in;
	endif
endfunction
