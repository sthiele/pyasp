% density.lp
% gene(G).
% redge(R1,R2).
% cat(G,R).
% start(R).
% goal(R).
% pmax(U).
#const pmin=1.

% genes that catalyse reactions
cat_gene(X) :- gene(X), cat(X,_).


% first check for possible paths 
% this is completely evaluated while grounding and hence
% can be used to skip the grounding of expensive rules
ppath(S):- start(S), cat_gene(G), cat(G,S).
ppath(Y) :- ppath(X), redge(X,Y), cat_gene(G), cat(G,Y).
pp:- ppath(X), goal(X).
:- not pp.

sgene(G) :- start(S),cat(G,S).
ggene(G) :- goal(S),cat(G,S).

last(L) :- gene(L), not gene(L+1).

% possible start an end genes have max distance pmax 
pse(S,E) :- pp, cat_gene(S),cat_gene(E), S<E, 
            (E-S)+1 <= U, pmax(U), (E-S)+1 >= pmin,
            sgene(B), S > B-U, E < B+U,
             B >=S, B <= E,
            ggene(G), S > G-U, E < G+U,
             G >=S, G <= E.

pse(S,S) :- pp, sgene(S), ggene(S), pmin==1.

% over the cut
pse(S,E) :- pp, cat_gene(S),cat_gene(E), S>E, 
            E < U, pmax(U), S > (L-U)+1, last(L),
            (L-S)+1+E <= U, (L-S)+1+E>= pmin.


% choose start and end genes
1{se(S,E): pse(S,E)}1.


% the length of the chosen dna strand
length((E-S)+1)   :- S<=E, se(S,E).

% over the cut
length((L-S)+1+E) :- S>E, se(S,E), last(L).


% the following genes are the activated
cgene(S) :- se(S,E).
cgene(E) :- se(S,E).

cgene(X) :- se(S,E), S<E, cat_gene(X), X>S, X<E. 

% over the cut
cgene(X) :- se(S,E), S>E, cat_gene(X), X>S. 
cgene(X) :- se(S,E), S>E, cat_gene(X), X<E. 


% active edges
aedge(X,Y) :- redge(X,Y), cat(G1,X), cat(G2,Y), cgene(G1), cgene(G2), G1!=G2.
aedge(X,Y) :- redge(X,Y),cat(G,X), cat(G,Y).


path(S) :- start(S), cgene(G), cat(G,S).
path(Y) :- path(X), aedge(X,Y).
:- not path(X), goal(X).

#minimize [  length(L) = L].

#hide.
#show cgene(X).



