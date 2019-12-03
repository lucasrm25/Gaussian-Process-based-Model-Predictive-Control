function r = ref(tk, xk, t_r, t_l)

%     xk = [2,1]';

    t_c = (t_r + t_l)/2;

    [~,idx] = min( pdist2(xk',t_c,'seuclidean',[1 1].^0.5).^2 );
    
    idx_target = idx +3;
    
    idx_target = mod(idx_target, size(t_c,1));
    
    r = t_c(idx_target,:);
end

