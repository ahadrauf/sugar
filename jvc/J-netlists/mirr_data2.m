tic;
net = cho_load('mirror25.m');  
cho_display(net);              
[f, egv, dq] = cho_mode(net);  
cho_modeshape(net, f, egv, dq, mode_num);
toc