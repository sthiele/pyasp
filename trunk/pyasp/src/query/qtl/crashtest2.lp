vertex(a;b;c;d;e;x).
edge(a,b).
obs_elabel(a,b,1).
edge(b,c).
obs_elabel(b,c,1).
edge(c,d).
obs_elabel(c,d,1).
edge(e,c).
obs_elabel(e,c,1).

obs_vlabel(exp1,d,1).
obs_vlabel(exp1,x,-1).


cand(a;b;e).
targ(c;d).

input(a;x;e).
