mata:
	step1 = 8
	x = (0..step1)/(step1) * 2*pi()
	v = sin(x)
	step2 = 32
	xq = (0..step2)/(step2) * 2*pi()
	spline3(x,v)
	v2=spline3(x,v)
	///v = v[1..rows(v)-1,.]
	spline3eval(spline3(x,v),xq)

	



end

/*
in matlab:

clear;
x = 0:pi/4:2*pi; 
v = sin(x);
xq = 0:pi/16:2*pi;
vq2 = interp1(x,v,xq,'spline')


*/
