%{
function ymodel = reconv(time, a, b, c, d, e, irf)
    L = length(time);
    y = a*exp(b*time)+c*exp(d*time) + e;
    ymodel_full = conv(y, irf, 'full');
    ymodel = ymodel_full(1:L);
end
%}

%
function ymodel = reconv(time, a, b, c, d, irf)
    L = length(time);
    y = a*exp(b*time)+c*exp(d*time);
    ymodel_full = conv(y, irf, 'full');
    ymodel = ymodel_full(1:L);
end
%}