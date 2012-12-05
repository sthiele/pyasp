subset(A,R,end,1,V) :- edge(A,V,R).
subset(A,RA,subset(B,RB,S),N+1,V) :- subset(A,RA,end,1,V), subset(B,RB,S,N,V), vertex(A) < vertex(B).

in(A,R,subset(A,R,S)) :- subset(A,R,S,N,V).
in(B,RB,subset(A,RA,S)) :- in(B,RB,S), subset(A,RA,S,N,V).

