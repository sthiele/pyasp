% first step compute minimal number of vertices to fix most of the targets:
% gringo-3.0.92 --body-order=unify minimal_fixator_sets.lp crashtest3.lp| clasp-2.0.1-st --restart-on-model --restarts=100
%
% second step compute all minimal fixator sets:
% gringo-3.0.92 --body-order=unify minimal_fixator_sets.lp crashtest3.lp| clasp-2.0.1-st --project --restart-on-model --restarts=100 --opt-all=0,2

% helper for input predicate
input(V) :- input(E,V).

input(P,V) :- input(V), exp(P).
input(V) :- vertex(V), not edge(U,V) : edge(U,V) : U != V.

% infer the experiment name from the observations
exp(E) :- obs_vlabel(E,V,S).
exp("empty") :- 0 { obs_vlabel(E,V,S)} 0.

converse(1,-1).
converse(-1,1).
sign(1).
sign(-1).


% set everything a candidate which is not  forbidden
%cand(X) :- input(X), not targ(X), not forbidden(X).


% every vertex has a unique sign
1 { vlabel(E,V,1), vlabel(E,V,-1) } 1 :- vertex(V), exp(E).

% unless stated otherwise, a regulation has a unique sign
1 { elabel(U,V,1), elabel(U,V,-1) } 1 :- edge(U,V), not obs_elabel(U,V,1),
                                                    not obs_elabel(U,V,-1).

vlabel(E,V,S) :- vertex(V), obs_vlabel(E,V,S).
elabel(U,V,S) :- edge(U,V), obs_elabel(U,V,S).

% influences
infl(E,V, 1)  :- elabel(U,V,S), vlabel(E,U,S), sign (S).
infl(E,V,-1)  :- elabel(U,V,S), vlabel(E,U,T), converse(S,T).

% every vertex labeling must be explained by one influence
:- vlabel(E,V,1), not infl(E,V,1), not input(V), exp(E), vertex(V).
:- vlabel(E,V,-1), not infl(E,V,-1), not input(V), exp(E), vertex(V).
%vlabel(E,V,0) :- not infl(E,V,1),  not infl(E,V,-1), not input(V), exp(E), vertex(V).


% each candidate vertext might be choosen to be in the minimal set
0{ccand(V)}1 :-  cand(V).

% compute everything which is fixed by the candidates
fixed(X) :- ccand(X).
fixed(X) :- fixed(Z): edge(Z,X), vertex(X), not n_fixed(X), 1{edge(Y,X)}.
n_fixed(X) :- infl(E,X, 1), infl(E,X, -1).

% every target must be fixed
%:- targ(V), not fixed(V).


% 2nd optimization priority: minimize the choosen candidates = minimal fixator sets
#minimize [ not false=0, ccand(V) : cand(V) ].
% 1st optimization priority: maximize fixed targets
#maximize [ not false=0, fixed(V) : targ(V) ].

#hide.
fixator(E,V,S) :- vlabel(E,V,S), ccand(V).
#show fixator/3.
fixed_target(X) :- targ(X), fixed(X).
#show fixed_target/1.

