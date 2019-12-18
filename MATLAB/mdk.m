

[M] = assemble_system(net, 'M', 5);
[D] = assemble_system(net, 'D', 5);
[K] = assemble_system(net, 'K', 5);

[L] = length(K);


[M] = assemble_system(net, 'M', 5);



L = length(K);
k = K;
k11 = k([1:2],[1:2]);
k12 = k([1:2],[3:L]);
k21 = k([3:L],[1:2]);
k22 = k([3:L],[3:L]);
k = k11 - k12*inv(k22)*k21;

L = length(K);
k = M;
k11 = k([1:1],[1:1]);
k12 = k([1:1],[2:L]);
k21 = k([2:L],[1:1]);
k22 = k([2:L],[2:L]);
k = k11 - k12*inv(k22)*k21;
m=k

L = length(K);
k = K;
k11 = k([1:1],[1:1]);
k12 = k([1:1],[2:L]);
k21 = k([2:L],[1:1]);
k22 = k([2:L],[2:L]);
k = k11 - k12*inv(k22)*k21;
k

