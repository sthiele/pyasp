vertex(a;b;ab;e;f).
edge(a,ab).
obs_elabel(a,ab,1).
edge(b,ab).
obs_elabel(b,ab,1).
edge(ab,e).
obs_elabel(ab,e,1).
edge(ab,f).
obs_elabel(ab,f,1).

obs_vlabel(exp1,e,1).
obs_vlabel(exp1,f,1).
exp(exp1).

cand(a;b;ab).
targ(e;f).
input(a;b).
