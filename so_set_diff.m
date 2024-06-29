function C = so_set_diff(A,B,eps=1e-06)
	if size(A, 1) != size(B, 1)
		error('Matrices must have the same number of rows');
	endif
	C = zeros(size(A,1),size(A,2));
	cnt = 1;
	for i = 1:size(A, 2)
		col_A = A(:, i);
		is_unique = true;
		for j = 1:size(B, 2)
			col_B = B(:, j);
			if vecnorm(col_A-col_B)<eps
				is_unique = false;
				break;
			endif
		endfor
		if is_unique
			C(:,cnt) = col_A;
			cnt=cnt+1;
		endif
	endfor
	C=C(:,1:cnt-1);
endfunction