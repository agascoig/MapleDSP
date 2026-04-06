                            "maple in ~/.mapleinit"

# Moon, "Mathematical Methods and Algorithms for Signal Processing" Eigenfilter Method
restart: 
omega__p := `*`(.4,Pi): # passband frequency;
omega__s := `*`(.6,Pi): # stopband frequency;
m := 13: # number of coefficients, must be odd;
alpha := .5: # tradeoff between stopband and passband;
M := iquo(m-1,2): 
p2 := `/`(Pi,2): 
eigmakePQ := " proc(wp,ws,m, $) ::list;"; _local(N,M,p2,mlist,nlist,n,P,Q,p2i,a,b); N := m-1; M := iquo(N,2); p2 := `/`(Pi,2); mlist := Vector(M,i -> i); P := Matrix(DiagonalMatrix(Concatenate(1,1-ws/Pi,((p2-`/`(ws,2)) -~ (sin~(2*mlist . ws) /~ ((4) . mlist)))/Pi))); Q := Matrix(DiagonalMatrix(Concatenate(1,0,((1.5*wp -~ ((2 *~ sin~(mlist . wp)) /~ mlist)) +~ sin~(((2) . mlist) . wp) /~ ((4) . mlist))/Pi))); p2i := `/`(1,2*Pi); " for n from 1 to M do "nlist := [seq(0 .. n-1)]; a := `-`~(n,nlist); b := `+`~(n,nlist); P[[n+1],`+`~(nlist,1)] := Vector[row]([-p2i*ws*(sinc~(a . ws,1) +~ sinc~(b . ws,1))]); P[`+`~(nlist,1),[n+1]] := convert(P[[n+1],`+`~(nlist,1)],Vector[column]); Q[[n+1],`+`~(nlist,1)] := Vector[row]([(wp/Pi) *~ ((((1-sinc(n*wp,1)) -~ sinc~(wp*nlist,1)) +~ .5*sinc~(wp*a,1)) +~ .5*sinc~(wp*b,1))]); Q[`+`~(nlist,1),[n+1]] := convert(Q[[n+1],`+`~(nlist,1)],Vector[column]); " end:": return evalf(P), evalf(Q); " end proc:":
(P, Q) := eigmakePQ(omega__p,omega__s,m): 
R := alpha . P+(1-alpha) . Q: 
(V, DD) := MTM:-eig(R): 
y := min(abs(Diagonal(DD))): 
i := min[index](abs(Diagonal(DD))): 
b := Re(V[() .. (),i]): 
mlist := [seq(1 .. M)]: 
h := Vector(m,i -> 0): 
h[(`-`~(M,mlist)) +~ 1] := `/`~(b[`+`~(mlist,1)],2): 
h[(`+`~(M,mlist)) +~ 1] := `/`~(b[`+`~(mlist,1)],2): 
h[M+1] := b[1]: 
HF := FilterFrequencyResponse(h,size = 1024,output = record,fftnormalization = none): 
HF[magnitudeplot];
HF[phaseplot];
Worksheet:-WorksheetToMapleText("filter_examples.maple");
help("WorksheetToMapleTest");

