function r = beta_divergence(A, B, beta)
switch beta
    case 0
        r = sum(sum (A./B - log(A./B+eps) - 1));
    case 1
        r = sum(sum(A.*log(A./B+eps) - A + B));
    case 2
        r = .5 * norm(A - B, 'fro').^2;
end
end