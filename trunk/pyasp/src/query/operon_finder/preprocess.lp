% redge/3
% start/1
% goal/1
% pmax(U)


start2(p(G,R)) :-  redge(p(G,R),Y,W), start(R).
goal2(p(G,R)) :-   redge(X,p(G,R),W), goal(R).

fpath(X,0) :- start2(X).
fpath(Y,L+W) :- redge(X,Y,W), fpath(X,L),
                pmax(N), L+W <= N,
                not start2(Y), not goal2(X).

pp :- fpath(X,_), goal2(X).
:- not pp.

bpath(X,0) :- goal2(X),pp.
bpath(X,L+W) :- redge(X,Y,W), bpath(Y,L),
                pmax(N), L+W <= N,
                not start2(Y), not goal2(X).

% Interesting edges lie in the intersection of fpath and bpath

pdir_edge(X,Y,W) :- redge(X,Y,W), not redge(Y,X,W), fpath(X,L1), bpath(Y,L2),
                pmax(N), L1+W <= N, L2+W <= N, L1+L2+W <= N,
                not start2(Y), not goal2(X).



pdir_edge(X,Y,W) :- redge(X,Y,W), redge(Y,X,W), fpath(X,L1), bpath(Y,L2),
                pmax(N), L1+W <= N, L2+W <= N, L1+L2+W <= N,
                not start2(Y), not goal2(X), start2(X).

pdir_edge(X,Y,W) :- redge(X,Y,W), redge(Y,X,W), fpath(X,L1), bpath(Y,L2),
                pmax(N), L1+W <= N, L2+W <= N, L1+L2+W <= N,
                not start2(Y), not goal2(X), goal2(Y).



prev_edge(X,Y,W) :- redge(X,Y,W), redge(Y,X,W), fpath(X,L1), bpath(Y,L2),
                pmax(N), L1+W <= N, L2+W <= N, L1+L2+W <= N,
                not start2(Y), not goal2(X), X<Y, not start2(X), not goal2(Y).

sedge(X,Y,W) :- pdir_edge(X,Y,W), start2(X).

gedge(X,Y,W) :- pdir_edge(X,Y,W), goal2(Y).


#hide.
#show pmax/1.
#show start2/1.
#show goal2/1.
#show pdir_edge/3.
#show prev_edge/3.
#show sedge/3.
#show gedge/3.
