function show(res,opt=struct())
	% Determine screen size
	
	if ~isfield(opt,'axes')
		opt.axes=[0 2 0 2];
	endif
	if ~isfield(opt,'dirlength')
		opt.dirlength=1;
	endif
	if ~isfield(opt,'name')
		opt.name='plot';
	endif
	
	screenSize = get(0, 'ScreenSize');
	screenWidth = screenSize(3);
	screenHeight = screenSize(4);

	% Set size for each figure window
	figWidth = 400;  % width of each figure window
	figHeight = 300; % height of each figure window
	
	
	
	
	cm=colormap(prism);
	X=[res.poi,res.dir];
	k=size(res.poi,2);
	for i=1:size(X,2)
		figPosX = mod((i-1) * figWidth, screenWidth);  % horizontal position
		figPosY = screenHeight - figHeight - mod(floor((i-1) * figWidth / screenWidth) * figHeight, screenHeight);  % vertical position
		fig=figure(i);
		set(gcf, 'Position', [figPosX, figPosY, figWidth, figHeight]);
		opt.color=[.7 .7 .7];
		if i <= k
			PP=res.PP;
		else
			PP=res.PP.recc;
		endif
		plot(PP,opt);
		hold on
		if i <= k
			opt.color=cm(1,:);
			gr=res.gr_F_std;
		else
			opt.color=cm(5,:);
			gr=res.gr_G_std;
		endif
		P=so_val(gr,X(:,i));
		P.plot(opt);
		hold off
		
		fullname='pics';
		axis(opt.axes);
		sprintf([fullname,"/",opt.name,"%d"],i)
		print (fig, sprintf([fullname,"/",opt.name,"%d"],i), "-dpdflatexstandalone");
		system (sprintf(["pdflatex --output-directory=",fullname," ",fullname,"/",opt.name,"%d"],i));		
		system(["rm ",fullname,"/*.aux"]);
		system(["rm ",fullname,"/*.log"]);		
	endfor
	
	
	disp('Press any key to terminate.');
	pause;
	close all;
	
	system(["inkscape --export-type=""pdf"" ",fullname,"/*.pdf"])
	system(['cd ',fullname,' && for file in *_out*; do mv "$file" "${file/_out/}"; done']);

	
endfunction


