function res = so_val(gr,x)
	n=size(x,1);
	q=gr.sdim-n;
	res = (gr & (point(x):space(q))).im([zeros(q,n),eye(q)]);
endfunction