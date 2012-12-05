% This code projects influence graphs to vertices that are properly constrained.
exp(P) :- obs_vlabel(P,V,S).
exp("empty") :- 0{obs_vlabel(P,V,S)}0.

observed(E,V) :- obs_vlabel(E,V,1).
observed(E,V) :- obs_vlabel(E,V,-1).

gets(E,V, 1) :- obs_elabel(U,V,S), obs_vlabel(E,U,S).
gets(E,V,-1) :- obs_elabel(U,V,S), obs_vlabel(E,U,T), S != T.

% Any variation of a vertex having a positive influence on itself is explained.
input(E,V) :- obs_elabel(V,V,1), exp(E).

% Any variation of a vertex receiving an unspecified influence is explainable.
input(E,V) :- edge(U,V),     obs_elabel(U,V,1),     obs_elabel(U,V,-1), exp(E).
input(E,V) :- edge(U,V), not obs_elabel(U,V,1), not obs_elabel(U,V,-1), exp(E).

% An observed vertex receiving a matching influence can be regarded as input.
input(E,V) :- obs_vlabel(E,V,S), gets(E,V,S).

% A vertex receiving complementary influences can be regarded as input.
input(E,V) :- gets(E,V,S), gets(E,V,T), S != T.

% An unobserved vertex without non-input successors is explainable.
%input4(E,V) :- edge(U,V), input(E,W) : edge(V,W), not observed(E,V), exp(E).

% An unobserved input vertex can explain any single non-input target.
input(E,V) :- edge(U,V), input(E,W) : edge(U,W) : V != W, input(E,U), not observed(E,U), exp(E).
#hide.
#show input/2.
%#show input4/2.
