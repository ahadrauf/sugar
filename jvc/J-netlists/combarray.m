subnet combarray [A B][nfinger=* wbase=* wfinger=* gap=* h=* Lfinger=*  V=*] 
[
    e0 = 8.854e-12
    beam_link_grav  p1 [A a(1)][l=2*(wfinger+gap)  w=wbase  h=h ] %first base
    for j=1:nfinger-1
    [
        beam_link_grav  p1 [a(j) a(j+1)][l=2*(wfinger+gap)  w=wbase  h=h ] %base between fingers
        beam_link_grav  p1 [a(j) b(j)][l=Lfinger  w=wfinger  h=h  oz=-pi/2  L1=wbase/2] %combfinger 
        f3d * [b(j)][F=2*1/2*e0*h/gap*V^2 oz=-pi/2] %combfinger force
    ]
    j=nfinger %last finger
    beam_link_grav  p1 [a(j) B][l=2*(wfinger+gap)  w=wbase  h=h ] %base between fingers
    beam_link_grav  p1 [a(j) b(j)][l=Lfinger  w=wfinger  h=h  oz=-pi/2  L1=wbase/2] %combfinger
    f3d * [b(j)][F=2*1/2*e0*h/gap*V^2 oz=-pi/2] %combfinger force
]

   
 