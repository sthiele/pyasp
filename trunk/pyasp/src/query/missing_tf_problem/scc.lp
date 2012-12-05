

%gcorrelation(GENA,GENB, X) :- correlation(GENA,GENB, X).
gcorrelation(GENA,GENB,-1) :- correlation(GENA,GENB,-1).
gcorrelation(GENA,GENB, 0) :- correlation(GENA,GENB, 0).
gcorrelation(GENA,GENB, 1) :- correlation(GENA,GENB, 1).

% when have two anti-/corellated genes a tf incommon
one_common_tf(GENA,GENB):- tf(GENA,TF), tf(GENB,TF), gcorrelation(GENA,GENB,1).
one_common_tf(GENA,GENB):- tf(GENA,TF), tf(GENB,TF), gcorrelation(GENA,GENB,-1).

% when have have two uncorrelated genes a distinct tf
one_different_tf(GENA,GENB) :- tf(GENA,TF), not tf(GENB,TF), gcorrelation(GENA,GENB,0).
one_different_tf(GENA,GENB) :- tf(GENB,TF), not tf(GENA,TF), gcorrelation(GENA,GENB,0).

% when have two uncorrelated genes all tf's incommon
all_common_tf(GENA,GENB) :- not one_different_tf(GENA,GENB), gcorrelation(GENA,GENB,0).


gene(X)  :- correlation(X,Y,Z).
gene(Y)  :- correlation(X,Y,Z).

gnFirst(X) :- gene(X), gene(Y), Y < X.
gFirst(X)  :- not gnFirst(X), gene(X).

gnLast(X) :- gene(X), gene(Y), Y > X.
gLast(X)  :- not gnLast(X), gene(X).

gnNext(X,Z) :- gene(X;Y;Z), X < Y, Y < Z, X < Z.
gNext(X,Y)  :- gene(X;Y), X < Y, not gnNext(X,Y).

numG(X):-  X = #sum[ gene(Y): gene(Y) ].

number(X,1) :- gFirst(X).
number(Y,N+1) :- number(X,N), gNext(X,Y),N=1..M-1,numG(M).


conflict(A,B)   :- not one_common_tf(A,B), gcorrelation(A,B,1).
conflict(A,B)   :- not one_common_tf(A,B), gcorrelation(A,B,-1).
conflict(A,B)   :- all_common_tf(A,B),  gcorrelation(A,B,0).
conflict(A,B)   :- conflict(B,A).

genc(X)  :- conflict(X,Y).
genc(Y)  :- conflict(X,Y).

scc(X,X)  :- genc(X).
scc(X,Z)  :- scc(X,Y), conflict(Y,Z), conflict(X,BLA).
scc(X,Z)  :- scc(X,Y), conflict(Z,Y), conflict(X,BLA).
%scc(X,Z)  :- scc(Y,X), conflict(Z,Y), conflict(Y,BLA).
mini2(X,M) :- M = #min[scc(X,Y): number(Y,N)= N], genc(X).
mini(M) :- M = #min[scc(X,Y): number(Y,N)= N], genc(X).
scc2(M,X) :- mini(M), number(Y,M), scc(Y,X).
scconflict(M,A,B) :- scc2(M,A), scc2(M,B), conflict(A,B).
scconflict(M,B,A) :- scc2(M,A), scc2(M,B), conflict(B,A).
numscconflict(M,Z) :-  Z = #sum[scconflict(M,A,B):scconflict(M,A,B) ], mini(M).

%1{ choose(X) : genc(X) }1.
%scc(X)  :- choose(X).
%scc(Y)  :- scc(X), conflict(X,Y).
%scc(X)  :- scc(Y), conflict(X,Y).



%corr(GENA,GENB,C) :- scc(GENA),scc(GENB), gcorrelation(GENA,GENB,C).


%numconflict(Z) :-  Z = #sum[conflict(A,B) : scc(A): genc(A) : scc(B): genc(B) ].
#hide.
%#show corr/3.
#show numconflict/1.
#show mini/1.
#show numscconflict/2.
%#show scc/2.
%#show conflict/2.
%#show scconflict/3.
%#show number/2.