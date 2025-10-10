function zhat = measFcn_attEuler(xk, ~)
zhat = [xk(1:3); xk(4:6)];
end