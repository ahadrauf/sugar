subnet u_contact [contact1 contact2 contact3 contact4][contact=* l=* d=* w=* q1=0 q2=0 q3=0 q4=0 q5=0 q6=0 q7=0 q8=0 q9=0 q10=0 q11=0 q12=0 q13=0 q14=0 q15=0 q16=0 q17=0 q18=0 q16=0 q17=0 q18=0 q19=0 q20=0 q21=0 q22=0 q23=0 q24=0]
[
    if contact==1
    [
        q_beam parent [contact1 contact2][l=l w=w oz=pi   qox1=q4  qoy1=q5  qoz1=q6  qx2=q7  qy2=q8  qz2=q9  qox2=q10 qoy2=q11 qoz2=q12]
        q_beam parent [contact2 contact3][l=d w=w oz=pi/2 qox1=q10 qoy1=q11 qoz1=q12 qx2=q13 qy2=q14 qz2=q15 qox2=q16 qoy2=q17 qoz2=q18]
        q_beam parent [contact3 contact4][l=l w=w oz=0    qox1=q16 qoy1=q17 qoz1=q18 qx2=q19 qy2=q20 qz2=q21 qox2=q22 qoy2=q23 qoz2=q24]
    ]
    if contact==0
    [
        q_beam p1 [contact1 contact2][l=l w=w oz=pi   qox1=0 qoy1=0 qoz1=0   qx2=0 qy2=0 qz2=0   qox2=0 qoy2=0 qoz2=0]
        q_beam p1 [contact2 contact3][l=d w=w oz=pi/2 qox1=0 qoy1=0 qoz1=0   qx2=0 qy2=0 qz2=0   qox2=0 qoy2=0 qoz2=0]
        q_beam p1 [contact3 contact4][l=l w=w oz=0    qox1=0 qoy1=0 qoz1=0   qx2=0 qy2=0 qz2=0   qox2=0 qoy2=0 qoz2=0]
    ]
]


