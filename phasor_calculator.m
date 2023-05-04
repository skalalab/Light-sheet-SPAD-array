function [M, Phi, G, S] = phasor_calculator(f, time, decays, IRF)
w = 2*pi*f;
G_decay = (decays*cos(w*time))./sum(decays,2);
S_decay = (decays*sin(w*time))./sum(decays,2);

G_IRF = sum(IRF.*cos(w*time))/sum(IRF(:));
S_IRF = sum(IRF.*sin(w*time))/sum(IRF(:));

GS = [G_IRF, S_IRF; -S_IRF, G_IRF]*[G_decay'; S_decay']/(G_IRF^2 + S_IRF^2);

G = transpose(GS(1,:));
S = transpose(GS(2,:));
Phi = pi-atan2(S,-G);
M = sqrt(G.^2 + S.^2);
end