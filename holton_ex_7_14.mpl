                            "maple in ~/.mapleinit"

Warning, on line 1, incomplete string; use " to end the string
# Holton DSP, p.467 Example 7.14 : Least-square-error lowpass filter
restart: 
N := 17: 
M := iquo(N-1,2): 
P := 41: 
R := `/`(P-1,2): 
wc := `*`(.5,Pi): 
k := Vector(M+1,i -> i-1)^%T: # row vector;
l := Vector(R+1,i -> i-1)^%T: # row vector;
w := 2*Pi*l^%H/P: # frequency samples, column vector;
D2 := Vector(R+1,i -> `if`(is(w[i] < wc or 2*Pi-w[i] < wc),1,0)): # amplitude samples, column vector;
C := cos~(evalf~(w . k)): # outer product;
a := LeastSquares(C,D2): 
h := Concatenate(1,.5*Reverse(a)[1 .. NumElems(a)-1],a[1],.5*a[2 .. ()]): # impulse response;
HF := FilterFrequencyResponse(h,size = 1024,output = record,fftnormalization = none): # noncausal response;
HF[magnitudeplot];
HF[phaseplot];
freq := HF[frequencies]: 
response := HF[response]: 
plot(freq,abs(response),title = "Absolute response magnitude",color = blue,gridlines = true,labels = ["Normalized Frequency (&times; &pi; rad/sample)", typeset("|H(",omega,")|")],labeldirections = [horizontal, vertical]);

