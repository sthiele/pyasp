%Generate boolean functions
{clause(subset(A,R,S),N,V)} :- subset(A,R,S,N,V).
%Avoid redundant functions such as: A or (A and B)
:- clause(S1,N,V), clause(S2,M,V), S1!=S2, in(A,R,S2) : in(A,R,S1).

%Initialization of the stimuli species.
active(E,V) :- exp(E,V,1), stimuli(V).
inactive(E,V) :- exp(E,V,0).

%Fixed species.
fixed(E,V) :- nexp(N), E=1..N, stimuli(V).
fixed(E,V) :- inhibitor(V), exp(E,V,0).

%Boolean logic
active(E,V) :- nexp(N), E=1..N, clause(S,M,V), not fixed(E,V),
               active(E,U)   : in(U,1,S),
               inactive(E,U) : in(U,-1,S).

%Species that are not active, are inactive.
inactive(E,V) :- vertex(V), nexp(N), E=1..N, not fixed(E,V), not active(E,V).

