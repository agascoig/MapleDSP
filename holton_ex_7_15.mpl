                            "maple in ~/.mapleinit"

# Holton DSP, p.469 Example 7.15 : Discrete weighted least-square-error lowpass filter
# The error of points in the stopband is weighted by a factor of K = 10.
restart: 
N := 17: 
R := 41: 
P := iquo(R-1,2): 
K := 10: 
wc := `*`(.5,Pi): 
k := Vector(iquo(N-1,2)+1,i -> i-1)^%T: # row vector;
l := Vector(P+1,i -> i-1)^%T: # row vector;
w := 2*Pi*l^%H/R: # frequency samples, column vector;
D2 := Vector(P+1,i -> `if`(is(w[i] < wc or 2*Pi-w[i] < wc),1,0)): # amplitude samples, column vector;
C := cos~(evalf~(w . k)): # outer product;
w_diag := Vector(P+1,i -> `if`(is(wc < w[i]),K,1)): 
W := DiagonalMatrix(w_diag): 
a := LeastSquares((C^%H . W) . C,(C^%H . W) . D2): 
h := Concatenate(1,.5*Reverse(a)[1 .. NumElems(a)-1],a[1],.5*a[2 .. ()]): # impulse response;
stemplot(h);
HF := FilterFrequencyResponse(h,size = 1024,output = record,fftnormalization = none): 
HF[magnitudeplot];
HF[phaseplot];
freq := HF[frequencies]: 
response := HF[response]: 
plot(freq,abs(response),title = "Absolute response magnitude",color = blue,gridlines = true,labels = ["Normalized Frequency (&times; &pi; rad/sample)", typeset("|H(",omega,")|")],labeldirections = [horizontal, vertical]);

