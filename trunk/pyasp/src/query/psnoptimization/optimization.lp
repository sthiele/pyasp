%Define mismatches
#const q = 4.
#const p = 2.
#const minsize = 0.
#const maxsize = -1.

dfactor(F) :- F = #pow(10, p).
penalty_D(D) :- D=#pow(10, q-2*p).

mismatch(E,V,F-D) :- obs(E,V,D), active(E,V), dfactor(F).
mismatch(E,V,D) :- obs(E,V,D), inactive(E,V).

%Variables to build the objective function
tsize(T) :- T = [subset(_,_,_,N,_) = N].
ne(E) :- E = [obs(_,_,_)].

:- maxsize + 1 [clause(_,N,_) = N], maxsize>=0.
:- T + 1 [clause(_,N,_) = N], maxsize<0, tsize(T).
:- [clause(_,N,_) = N]minsize - 1.

#minimize[ mismatch(_,_,W)   : tsize(T) : penalty_D(D) = #pow(W,2)*T*D,
           clause(_,N,_)     :    ne(E) = E*N].

