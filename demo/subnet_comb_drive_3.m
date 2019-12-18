%subnet_comb_drive_3
%This subnet makes a comb drive pair
%jvclark - Feb 2011
%   ==========A==============
%     |   |   |   |   |   |
%   |   |   |   |   |   |   |
%   ==========B==============

subnet subnet_comb_drive_3 [A B][ V=* Nf=* gapf=* Lf=* Wf=* Hf=* overlap=* W_support=* ] 
[
    e0 = 8.854e-12
    F = 2 * 1/2 * e0 * V^2 * Hf / gapf
    n = 1
    L_support = gapf*2 + Wf

    beam3d parent [A tt(0+n)][l=W_support/2 w=W_support h=Hf oz=-pi/2]
    
    %Top center finger with sides
    beam3dlink parent [tt(0+n) tf(0+n)][l=Lf w=Wf oz=-pi/2 L1=W_support/2 h=Hf]         
    beam3d parent [tt(0+n) ttr(1+n)][l=L_support w=W_support h=Hf]
    beam3d parent [tt(0+n) ttl(1+n)][l=L_support w=W_support h=Hf oz=pi]
    %Bottom center u-shape
    beam3d parent [btl(0+n) B][l=L_support/2 w=W_support h=Hf]
    beam3d parent [B btr(0+n)][l=L_support/2 w=W_support h=Hf]
    beam3dlink parent [btr(0+n) bfr(0+n)][l=Lf w=Wf oz=pi/2 L1=W_support/2 h=Hf] 
    beam3dlink parent [btl(0+n) bfl(0+n)][l=Lf w=Wf oz=pi/2 L1=W_support/2 h=Hf]         
    for n = 1 : Nf
    [
        %Top right side
        beam3d parent [ttr(0+n) ttr(1+n)][l=L_support w=W_support h=Hf]
        beam3dlink parent [ttr(1+n) tfr(1+n)][l=Lf w=Wf oz=-pi/2 L1=W_support/2 h=Hf] 
        f3d * [tfr(1+n)] [F=F oz=-pi/2] 

        %Top left side
        beam3d parent [ttl(0+n) ttl(1+n)][l=L_support w=W_support h=Hf oz=pi]
        beam3dlink parent [ttl(1+n) tfl(1+n)][l=Lf w=Wf oz=-pi/2 L1=W_support/2 h=Hf] 
        f3d * [tfl(1+n)] [F=F oz=-pi/2] 

        %Bottom right
        beam3d parent [btr(0+n) btr(1+n)][l=L_support w=W_support h=Hf]
        beam3dlink parent [btr(1+n) bfr(1+n)][l=Lf w=Wf oz=pi/2 L1=W_support/2 h=Hf] 
        f3d * [btr(1+n)] [F=F oz=pi/2] 

        %Top side
        beam3d parent [btl(0+n) btl(1+n)][l=L_support w=W_support h=Hf oz=pi]
        beam3dlink parent [btl(1+n) bfl(1+n)][l=Lf w=Wf oz=pi/2 L1=W_support/2 h=Hf] 
        f3d * [btl(1+n)] [F=F oz=pi/2] 
    ]   

    %Relative positioning
    Lx = (Nf*2)*(gapf*2+Wf) + (gapf*2+Wf)*0.5
    Ly = overlap+W_support/2
    ww = 2n
    beam3dlink_cloak parent [tfr(Nf+1) x1][l=Lx w=ww h=ww oz=pi]    
    beam3dlink_cloak parent [x1 bfl(Nf+1)][l=Ly w=ww h=ww oz=pi/2]    
    beam3dlink_cloak parent [tfl(Nf+1) x2][l=Lx w=ww h=ww ]    
    beam3dlink_cloak parent [x2 bfr(Nf+1)][l=Ly w=ww h=ww oz=pi/2]    

    %Beam caps
    beam3d parent [ttr(0+Nf+1) ttr(1+Nf+1)][l=Wf/2 w=W_support h=Hf]
    beam3d parent [ttl(0+Nf+1) ttl(1+Nf+1)][l=Wf/2 w=W_support h=Hf oz=pi]
    beam3d parent [btr(0+Nf+1) btr(1+Nf+1)][l=Wf/2 w=W_support h=Hf]
    beam3d parent [btl(0+Nf+1) btl(1+Nf+1)][l=Wf/2 w=W_support h=Hf oz=pi]
]

