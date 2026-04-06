                            "maple in ~/.mapleinit"

Warning, on line 1, incomplete string; use " to end the string
# Holton Example 7.20: Integral Least-Square Error Lowpass Filter
# Design an integral least-square-error lowpass filter of length N = 17 with a passband of \317\211
# p
# = 0.45\317\200 and a stopband of \317\211
# s
# = 0.55\317\200, where the error of the passband is weighted by one, and the error in the stopband is weighted by a factor of K = 10.
restart: 
with(LinearAlgebra): with(SignalProcessing): with(ArrayTools):
assume(n,integer): 
omega__p := `*`(.45,Pi): 
omega__s := `*`(.55,Pi): 
omega__p := `*`(.45,Pi): 
omega__s := `*`(.55,Pi): 
N := 17: 
M := iquo(N-1,2): 
K := 10: 
n1 := Vector(M+1,i -> i-1)^%T: 
n2 := Vector(2*M+1,i -> i-1)^%T: 
f := ptmp: 
d := f(omega__p,n1): 
omega__v := f(omega__p,n2)-(K . f(omega__s,n2)): 
omega__v[1] := omega__v[1]+K: 
omega__1 := Matrix(M+1,(i, j) -> omega__v[i+j-1]): # Hankel Matrix;
omega__2 := Matrix(M+1,(i, j) -> `if`(is(i < j),omega__v[j-i+1],omega__v[i-j+1])): # Toeplitz;
a := LeastSquares(`*`(.5,omega__1+omega__2),d^%H): 
h := Concatenate(1,.5*Reverse(a)[1 .. NumElems(a)-1],a[1],.5*a[2 .. ()]): 
HF := FilterFrequencyResponse(h,size = 1024,output = record,fftnormalization = none): 
HF[magnitudeplot];
HF[phaseplot];
freq := HF[frequencies]: 
response := HF[response]: 
plot(freq,abs(response),title = "Absolute response magnitude",color = blue,gridlines = true,labels = ["Normalized Frequency (&times; &pi; rad/sample)", typeset("|H(",omega,")|")],labeldirections = [horizontal, vertical]);

