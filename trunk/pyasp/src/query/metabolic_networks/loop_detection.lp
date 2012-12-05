

reaction(rev(R)) :- reversible(R), reaction(R).
reactant(S,rev(R)) :- reversible(R), reaction(R), product(S,R).
product(S,rev(R)) :- reversible(R), reaction(R), reactant(S,R).

species(S) :- product(S,R).
species(S) :- reactant(S,R).

%numSpecies(X)  :- X=#sum [species(S)].
%#show numSpecies/1.

0 {isinloop(S) : species(S)}1.

loop(S2) :- reactant(S1,R), product(S2,R), isinloop(S1).
loop(S2) :- reactant(S1,R), product(S2,R), loop(S1).

:- isinloop(S), not loop(S).


#hide.
#show isinloop/1.