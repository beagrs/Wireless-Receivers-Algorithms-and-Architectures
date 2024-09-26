function [start] = detector(p,r, thr)
% Input :   p: preamble (shape = (100, 1)), r: received signal (shape = (Nr, 1)), thr: scalar threshold 
% output:   start: signal start index
Np = size(p,1);
Nr = size(r,1);
c = zeros(Nr-Np, 1);
c_norm = zeros(Nr-Np, 1);
%% TODO
for n = 1:length(c)
    c(n) = sum(((p').').*r(n:(n+Np-1),1));
    c_norm(n) = abs(c(n)).^2/sum(abs((r(n:(n+Np-1),1))).^2);
    if(c_norm(n)) >thr
        start = n+length(p);
        break
    end
end
%after loop no threshold reached
if start == 0
    disp('Frame start not found.')
    start = -1;
end
end

