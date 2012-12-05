% first step compute minimal size for an explanation:
% gringo-3.0.92 --body-order=unify minimal_explanations.lp quick-explain.lp | clasp-2.0.1-st --restart-on-model --restarts=100
%
% second step compute all minimal explanations:
% gringo-3.0.92 --body-order=unify minimal_explanations.lp quick-explain.lp | clasp-2.0.1-st --restart-on-model --restarts=100 --project --opt-all=25 

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



% compute all targets which are possibly explained by the candidates
pexplained(V) :- cand(V).
pexplained(V) :- vlabel(E,V, 1), pexplained(U), elabel(U,V,S), vlabel(E,U,S).
pexplained(V) :- vlabel(E,V,-1), pexplained(U), elabel(U,V,S), vlabel(E,U,T), converse(S,T).

% each candidate vertext might be choosen to be in the minimal set
0{ccand(V)}1 :-  cand(V).

% compute all targets which are explained by the choosen candidates
explained(V) :- ccand(V).
explained(V) :- vlabel(E,V, 1), explained(U), elabel(U,V,S), vlabel(E,U,S).
explained(V) :- vlabel(E,V,-1), explained(U), elabel(U,V,S), vlabel(E,U,T), converse(S,T).

% every explainable target must be explained
:- targ(V), not explained(V), pexplained(V).


% 1st optimization priority: minimize the choosen candidates = small explanation
#minimize [ not false=0, ccand(V) : cand(V) ].


#hide.
explanation(E,V,S) :- vlabel(E,V,S), ccand(V).
#show explanation/3.

%explained_targets(X) :- targ(X), explained(X).
%#show explained_targets/1.
