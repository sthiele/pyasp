shared(S) :- reactant(S,R1), reactant(S,R2), R1 < R2.
shared(P) :- product(P,R1), reactant(P,R2), R1 < R2, reversible(R1).
shared(P) :- product(P,R1), product(P,R2), R1 < R2, reversible(R1), reversible(R2).
%numShared(X) :- X = #sum [shared(S) : species(S)].

%quantity of S
vertex(gen(q(S))) :- reactant(S,R).
%flux of S not produced by R
vertex(gen(if(S,R))) :- reactant(S,R),     isinloop(S).
vertex(gen(if(S,R))) :- reactant(S,R),                      shared(S).
vertex(gen(f(S)))    :- reactant(S,R), not isinloop(S), not shared(S).

%quantity of P
vertex(gen(q(P))) :- product(P,R).
%flux of P produced by R
vertex(gen(of(P,R))) :- product(P,R),                      shared(P).
vertex(gen(of(P,R))) :- product(P,R),     isinloop(P). 
vertex(gen(f(P)))    :- product(P,R), not isinloop(P), not shared(P). 

%reaction speed
vertex(gen(rs(R))) :- reactant(S,R).


% fluxes of reactants increases reaction speed
      edge(gen(f(S)),gen(rs(R)))   :- reactant(S,R), vertex(gen(f(S))).
obs_elabel(gen(f(S)),gen(rs(R)),1) :- reactant(S,R), vertex(gen(f(S))).

      edge(gen(f(P)),gen(rs(R)))   :- product(P,R), vertex(gen(f(P))),reversible(R).
obs_elabel(gen(f(P)),gen(rs(R)),1) :- product(P,R), vertex(gen(f(P))),reversible(R).

      edge(gen(if(S,R)),gen(rs(R)))   :- reactant(S,R), vertex(gen(if(S,R))).
obs_elabel(gen(if(S,R)),gen(rs(R)),1) :- reactant(S,R), vertex(gen(if(S,R))).

      edge(gen(of(P,R)),gen(rs(R)))   :- product(P,R), vertex(gen(of(P,R))),reversible(R).
obs_elabel(gen(of(P,R)),gen(rs(R)),1) :- product(P,R), vertex(gen(of(P,R))),reversible(R).

% reaction speed decreases quantity of reactants
      edge(gen(rs(R)),gen(q(S)))    :- reactant(S,R).
obs_elabel(gen(rs(R)),gen(q(S)),-1) :- reactant(S,R).

      edge(gen(rs(R)),gen(q(P)))    :- product(P,R),reversible(R).
obs_elabel(gen(rs(R)),gen(q(P)),-1) :- product(P,R),reversible(R).

% reaction speed increases flux of product
      edge(gen(rs(R)),gen(f(P)))   :- product(P,R), vertex(gen(f(P))).
obs_elabel(gen(rs(R)),gen(f(P)),1) :- product(P,R), vertex(gen(f(P))).

      edge(gen(rs(R)),gen(f(S)))   :- reactant(S,R), vertex(gen(f(S))),reversible(R).
obs_elabel(gen(rs(R)),gen(f(S)),1) :- reactant(S,R), vertex(gen(f(S))),reversible(R).

      edge(gen(rs(R)),gen(of(P,R)))   :- product(P,R), vertex(gen(of(P,R))).
obs_elabel(gen(rs(R)),gen(of(P,R)),1) :- product(P,R), vertex(gen(of(P,R))).

      edge(gen(rs(R)),gen(if(S,R)))   :- reactant(S,R), vertex(gen(if(S,R))),reversible(R).
obs_elabel(gen(rs(R)),gen(if(S,R)),1) :- reactant(S,R), vertex(gen(if(S,R))),reversible(R).

% reaction speed increases quantity of product
      edge(gen(rs(R)),gen(q(P)))   :- product(P,R).
obs_elabel(gen(rs(R)),gen(q(P)),1) :- product(P,R).

      edge(gen(rs(R)),gen(q(S)))   :- reactant(S,R),reversible(R).
obs_elabel(gen(rs(R)),gen(q(S)),1) :- reactant(S,R),reversible(R).


% balance 
% a balance node for a shared substrate S, used by reaction R1 "if(S,R1)" 
% has a negative influence on all the other balance nodes for the same substrate S used by another reactions R2 "if(S,R2)"

vertex(gen(bn(S,R))) :- reactant(S,R), shared(S).
vertex(gen(bn(P,R))) :- product(P,R), shared(P), reversible(R).   
    
      edge(gen(if(S,R)),gen(bn(S,R)))   :- reactant(S,R), vertex(gen(bn(S,R))).
obs_elabel(gen(if(S,R)),gen(bn(S,R)),1) :- reactant(S,R), vertex(gen(bn(S,R))).

      edge(gen(of(P,R)),gen(bn(P,R)))   :- product(P,R), vertex(gen(bn(P,R))),reversible(R).
obs_elabel(gen(of(P,R)),gen(bn(P,R)),1) :- product(P,R), vertex(gen(bn(P,R))),reversible(R).

      edge(gen(bn(S,R1)),gen(if(S,R1)))    :- reactant(S,R),vertex(gen(bn(S,R1))), R1!=R.
obs_elabel(gen(bn(S,R1)),gen(if(S,R1)),-1) :- reactant(S,R),vertex(gen(bn(S,R1))), R1!=R.

      edge(gen(bn(P,R1)),gen(of(P,R1)))    :- product(P,R),vertex(gen(bn(P,R1))), R1!=R.
obs_elabel(gen(bn(P,R1)),gen(of(P,R1)),-1) :- product(P,R),vertex(gen(bn(P,R1))), R1!=R.


%connect output fluxes of different reactions to input flux
      edge(gen(of(P,R2)),gen(if(P,R1)))   :- vertex(gen(of(P,R1))),vertex(gen(if(P,R2))), R1!=R2.
obs_elabel(gen(of(P,R2)),gen(if(P,R1)),1) :- vertex(gen(of(P,R1))),vertex(gen(if(P,R2))), R1!=R2.

      edge(gen(if(P,R2)),gen(of(P,R1)))   :- vertex(gen(of(P,R1))),vertex(gen(if(P,R2))), R1!=R2,reversible(R).
obs_elabel(gen(if(P,R2)),gen(of(P,R1)),1) :- vertex(gen(of(P,R1))),vertex(gen(if(P,R2))), R1!=R2,reversible(R).


#hide.
#show vertex/1.
#show edge/2.
#show obs_elabel/3.



