% Build action natrix for variable "actMvar" given coefficient matrix M
% (GBsolver subroutine)
% by Martin Bujnak, mar2008


function [amrows amcols gjcols aidx] = gbs_BuildActionMatrix(cfg, M, algB, amLT, amLTall, algBidx, actMvar)

    % extract action matrix for variable "actMvar"
    
    cols = size(M, 2);
    algBcnt = size(algBidx, 2);

    amrows = zeros(algBcnt, 1);
    
    % GJ
    nonzero = find(sum(M) ~= 0);
    Kk = M(:, nonzero);
    B = gjzpsp(Kk, cfg.prime);
    A = zeros(size(M));
    A(:, nonzero) = B;

    %
    % detect leading '1'
    [val idx] = min(abs(B' - 1));
    idx = idx(find(val == 0));
    
    % add algB columns

    % columns that are really required in GJ
    gjcols = union(nonzero(idx), (cols+1)-algBidx);
        
    % create monomials lookup table
    aidx = (-1)*ones(1,cols);
    aidx(gjcols) = 1:size(gjcols, 2);
    
    if 0
        % debug...
        [r_ones c_ones]=find(A(:,cols-amLT+1) == 1);
        rowssum = find(sum(A(r_ones, zero_el),2) == 1);
        if size(r_ones,1) < amLTcnt || size(rowssum, 1) < amLTcnt
            
            fprintf('error !! wrong ellimination\n');
        end
    end

    fprintf('extracting the action matrix for variable ''%s'' (used coef. matrix size %dx%d) \n', char(actMvar), size(Kk, 1), length(gjcols));
    
    amcols = aidx((cols+1)-algBidx ) ;
    
    am = zeros(algBcnt);
    for i=1:algBcnt

        pos = find(algBidx == amLTall(i));
        if isempty(pos)
            
            % take remainder
            [r_one]=find(A(:,cols-amLTall(i)+1) == 1);
            amrows(i) = r_one;

        else
            % in the B
            amrows(i) = -pos;
        end
    end
    
    if 0
        
        Kk = M(:, gjcols);
        B = rrefZp(Kk, cfg.prime);
        A2 = zeros(size(M));
        A2(:, gjcols) = B;


        am2 = zeros(algBcnt);
        for i=1:algBcnt

            pos = find(algBidx == amLTall(i));
            if isempty(pos)

                % take remainder
                [r_one]=find(A2(:,cols-amLTall(i)+1) == 1);
                am2(i, :) = -A2(r_one, (cols+1)-algBidx);
            else

                % in the B
                am2(i, pos) = 1;
            end
        end
    end

end