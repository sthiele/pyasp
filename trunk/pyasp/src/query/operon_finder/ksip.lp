% input predicates #shown by preprocess.lp
% pdir_edge/3.
% prev_edge/3.
% sedge/3.
% gedge/3.
% start2/1
% goal2/1
% pmax(U)

pmin(min).

1{ cdir_edge(X,Y,W) : sedge(X,Y,W) }1.
1{ cdir_edge(X,Y,W) : gedge(X,Y,W) }1.

0{ cdir_edge(X,Y,W) }1 :- pdir_edge(X,Y,W), not sedge(X,Y,W), not gedge(X,Y,W).
0{ crev_edge(X,Y,W) }1 :- prev_edge(X,Y,W).


% Can we reach a goal node
at(S) :- start2(S).
at(Y) :- cdir_edge(X,Y,W), at(X).
at(Y) :- crev_edge(X,Y,W), at(X).
at(Y) :- crev_edge(Y,X,W), at(X).
path :- at(G), goal2(G). 
:- not path.

% Everything longer than pmax is no solution
:-  pmax(U), U+1 [ cdir_edge(X,Y,W)=W, crev_edge(X2,Y2,W2)=W2 ].
% Everything smaller than pmin is no solution
:-  pmin(L), [ cdir_edge(X,Y,W)=W, crev_edge(X2,Y2,W2)=W2 ] L-1.

%additional constraints
fedge(X,Y) :- cdir_edge(X,Y,W).
fedge(X,Y) :- crev_edge(X,Y,W).

% no branches
:-fedge(A,X),fedge(A,Z), fedge(A,Y), X!=Y, X!=Z, Y!=Z.
:-fedge(X,A),fedge(A,Z), fedge(A,Y), X!=Y, X!=Z, Y!=Z.
:-fedge(X,A),fedge(Z,A), fedge(A,Y), X!=Y, X!=Z, Y!=Z.
:-fedge(X,A),fedge(Z,A), fedge(Y,A), X!=Y, X!=Z, Y!=Z.

% no unconnected edges
:-fedge(X,Y), not at(X).
:-fedge(X,Y), not at(Y).


% restore direction of found edges
fedge(X,Y,W) :- cdir_edge(X,Y,W).

fedge(X,Y,W) :- fedge(Z,X,_), crev_edge(X,Y,W), Y!=Z.
fedge(X,Y,W) :- fedge(Z,X,_), crev_edge(Y,X,W), Y!=Z.


% Find shortest paths
#minimize [cdir_edge(X,Y,W)=W,crev_edge(X2,Y2,W2)=W2].


#hide.
#show fedge/3.
