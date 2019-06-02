% Generate all polynomials required to build an action matrix for
% variable "actMvar" (detect actMvar if variable was not specified)
% (GBsolver subroutine)
% by Martin Bujnak, mar2008


function [M Mcoefs symcoefs amVar amLT amLTall algBidx algB] = gbs_PreparePolySystemCoefsMatrix(cfg, eq, known, unknown, algB, amVar)

    bdoReduce = cfg.bdoReduce;
    ordering = cfg.ordering;
    prime = cfg.prime;

    % create random instance
    if ~isfield(cfg, 'eqinstance')

        cfg.eqinstance = cfg.InstanceGenerator(cfg, eq, known, unknown);
    end
    
    if nargin < 5 || isempty(algB)
        
        % find quotion ring basis ( monomial basis of quotient ring C[x1,...xn]\I )
        [algB] = cfg.GBSolver(cfg, eq, known, unknown);
        
        if isempty(algB)
            
            M=[];
            Mcoefs=[];
            symcoefs=[];
            amVar=[];
            amLT=[];
            amLTall=[];
            algBidx=[];
            return;
        end
    end
    
    % prepare algebra B and find which monomials we need to build the
    % action matrix
    fprintf('analyzing quotient ring basis\n');
    
    if nargin < 9
        
        amVars = unknown;
    else
        amVars = amVar;
    end
    
    amStats = {};
    for i=1:length(amVars)
        
        actMvar = amVars{i};
        
        [amLT amLTcnt amLTall algB algBidx] = gbs_ParseAlgebraB(algB, actMvar, unknown);
        
        amStats{i}.actMvar = actMvar;
        amStats{i}.amLT = amLT;
        amStats{i}.amLTcnt = amLTcnt;
        amStats{i}.amLTall = amLTall;
        amStats{i}.algB = algB;
        amStats{i}.algBidx = algBidx;
    end
    
    algBcnt = size(algB, 2);

    %
    % extract equations into coeffients + monomial form
    fprintf('extracting coefficient & monomials\n');
    [monomials p p_currdeg symcoefs maxdeg] = gbs_ParseEquations(cfg, eq, known, unknown);
    
    
    %
    % just statistics (can be removed from this code)
    [monomials] = ReorderMonomials(monomials, unknown, ordering);
    
    fprintf('...used %d monomials : \n   ', length(monomials));
    for mons = monomials
        fprintf('%s ', char(mons));
    end
    fprintf('\n');

    fprintf('...max. poly. degree %d\n', maxdeg);

    
    if ~isempty(cfg.crashlog)

        fprintf('...crash log file detected - recovering from previous state using: %s\n', cfg.crashlog);
        [used lstused max_deg] = gbs_ResetFromLog(cfg.crashlog);

        maxdeg = max_deg;
        fprintf('...need to generate polynomials to degree %d\n', maxdeg);
    end
    

    %
    % generate monomials up-to maxdeg (matrix cols)
    allmons = {'1'};
    alldegs = zeros(1, length(unknown));
    allmonsdeg = 0;
    for i=1:maxdeg
        [mons degs] = GenerateMonomials(i, unknown, ordering);
        allmons = [mons allmons];
        alldegs = [degs; alldegs];
        allmonsdeg = [i*ones(1, length(mons)) allmonsdeg];
    end
    
    % create action var. masks
    cols = size(allmonsdeg, 2);
    for a=1:length(amStats)
        amStats{a}.zero_el = setdiff(1:cols, (cols+1)-amStats{a}.algBidx);
    end
    
    if cfg.exportEqs
        
        fprintf('...exporting parsed equations');
        
        % debug - create coefficient matrix with all initial polynomials
        cols = size(allmonsdeg, 2);
        Minit = zeros(length(eq), cols);
        for i=1:length(eq)
            for j=1:p{i}.monscnt
                pmon2 = p{i}.deg(j,:);
                order = GetMonomialOrder(pmon2, unknown); 
                Minit(i, cols-(order-1)) =  p{i}.coefsIDX(j);
            end        
        end
        
        nonzero = find(sum(Minit) ~= 0);
        Mmons = Minit(:, nonzero);
        
        gbs_DBGExport(known, unknown, eq, symcoefs, allmons, Minit, Mmons);
    end

    fprintf('adding polynomials\n');
    fprintf('...initializing matrix size %d,%d\n', length(eq), length(allmons));

    % initial size
    cols = size(allmonsdeg, 2);
    rowsused = 0;
    rowsalloc = length(allmons);
    M = zeros(rowsalloc, cols);
    Mcoefs  = zeros(rowsalloc, cols);

    % generate polynomials
    allOk = false;
    testeddeg = allmonsdeg(1);

    bAllEqReady = false;

    while ~allOk

        % add monomial multiples up to 'testeddeg' degree
        while ~bAllEqReady

            [v eqi] = min(p_currdeg);
            if v >= testeddeg 
                break;
            end

            % add monomial multiples 
            todeg = v + 1;
            if (todeg < p{eqi}.maxdeg)
                todeg = p{eqi}.maxdeg;
            end

            if (todeg > alldegs(1))

                % reallocate cols
                [mons degs] = GenerateMonomials(todeg, unknown, ordering);
                allmons = [mons allmons];
                alldegs = [degs; alldegs];
                allmonsdeg = [todeg*ones(1, length(mons)) allmonsdeg];            

                [M Mcoefs] = gbs_ReallocMatrix(M, Mcoefs, rowsalloc, cols+length(mons));

                cols = size(allmonsdeg, 2);

                for a=1:length(amStats)
                    amStats{a}.zero_el = setdiff(1:cols, (cols+1)-amStats{a}.algBidx);
                end
                
                fprintf('...coefficient matrix reallocaction, adding %d(%d degree) monomials\n', length(mons), todeg);
            end

            starmon = p{eqi}.monsidx;
            while (p{eqi}.maxdeg+allmonsdeg(cols-starmon)) <= todeg

                mon = alldegs(cols-starmon, :);

                % reallocate matrix on demand
                if rowsused == rowsalloc
                    
                    rowsalloc = 2*rowsalloc;
                    [M Mcoefs] = gbs_ReallocMatrix(M, Mcoefs, rowsalloc, cols);
                end

                % multiply eqi-th equation by monomial
                rowsused = rowsused + 1;
                for j=1:p{eqi}.monscnt

                    pmon2 = p{eqi}.deg(j,:) + mon;
                    order = GetMonomialOrder(pmon2, unknown); 
                    M(rowsused, cols-(order-1)) =  p{eqi}.coefs(j);
                    Mcoefs(rowsused, cols-(order-1)) =  p{eqi}.coefsIDX(j);
                end
                
                starmon = starmon + 1;
            end
            p_currdeg(eqi) = p{eqi}.maxdeg + allmonsdeg(cols-(starmon-1));
            p{eqi}.monsidx = starmon;

        end

        fprintf('...reducing matrix with size %dx%d\n', rowsused, cols);

        % reduce not required polynomials (if we generated too much)
        allOk = true;
        
        reduced_tail = rowsused;
        reduce_step = 1;
        
        if ~isempty(cfg.crashlog) && ~isempty(used)
            
            reduced_tail = lstused;
            selrows=[1:lstused fliplr(used')];
            prevState = true;
        else
            selrows = 1:rowsused;
            prevState = false;
        end

        reduced_final = selrows;
        
        while allOk

            [res] = gbs_CheckActionMatrixConditions(M, selrows, amStats, false, prime);
            allOk = res > 0;
            
            if allOk

                amRes = res;
                amVar = amStats{res}.actMvar;
                
                if ~prevState              
                    
                    if ~bdoReduce
                        fprintf('\n');
                        break;
                    end
                    
                    fprintf('all required polynomials for action variable ''%s'' generated, removing unnecessary\n', char(amVar));
                    reduce_step = floor(rowsused / 32);
                    if (reduce_step < 4) 
                        reduce_step = 4;
                    end;
                else
                    reduce_step = 2*reduce_step;
                end

                if allOk && prevState
                    fprintf('succeeded\n');
                end

                % reduction phase starts here...
                prevState = true;                   %set start state
                reduced_final = selrows;            %current OK polynomials

            elseif prevState

                fprintf('failed *\n');
                
                % undo step
                if reduce_step > 1
                    
                    reduced_tail = reduced_tail + reduce_step;
                    reduce_step = floor(reduce_step / 4);
                    
                    if (reduce_step < 1) 
                        reduce_step = 1; 
                    end;
                end
 
                % need more polynomials
                selrows = reduced_final;
                allOk = prevState;
            end

            if reduced_tail <= 0 || reduce_step == 0

                break;
            end
            
            % remove another equation
            if (reduced_tail < reduce_step) 
                reduce_step = reduced_tail - 1;
            end
            
            if 0
                % remove starting from 1
                filter = (size(selrows, 2) + 1) - ((reduced_tail-reduce_step+1):reduced_tail);
            elseif 0
                % random remove
                filter = ceil(rand() * reduced_tail);
                reduce_step = 1;
                fprintf('...removing equation %d ', filter);
            else
                filter = (reduced_tail-reduce_step+1):reduced_tail;
            end
            
            selrows = setdiff(selrows, filter);

            if reduced_tail > 0

                if reduce_step < 2
                    fprintf('...removing equation %d ', reduced_tail);
                else
                    fprintf('...removing equation %d - %d ', (reduced_tail-reduce_step+1), reduced_tail);
                end
            end
            
            reduced_tail = reduced_tail - reduce_step;
        end
        
        if (~allOk)

            fprintf('failed +\n');
        end

        rowsused = length(reduced_final);
        M(1:rowsused, :) = M(reduced_final, :);
        Mcoefs(1:rowsused, :) = Mcoefs(reduced_final, :);
        testeddeg = testeddeg + 1;
    end

    if ~allOk
        
        M = [];
        Mcoefs = [];
        
        amVar = [];
        amLT = [];
        amLTall = [];
        algBidx = [];

    else
        
        amVar = amStats{amRes}.actMvar;
        amLT = amStats{amRes}.amLT;
        amLTall = amStats{amRes}.amLTall;
        algBidx = amStats{amRes}.algBidx;

        M = M(1:rowsused, :);
        Mcoefs = Mcoefs(1:rowsused, :);

        nonzero = find(sum(M) ~= 0);

        fprintf('generated %d polynomials in %d monomials (%d nonzero)\n', rowsused, size(M, 2), length(nonzero));
    end
end

function [res] = gbs_CheckActionMatrixConditions(M, selrows, amStats, isElim, prime)

    if ~isElim

        nonzero = find(sum(M) ~= 0);
        Kk = M(selrows, nonzero);

        gjtime = cputime;

        %B2 = gjzp(Kk, prime);
        %B = rrefZp(Kk, prime);
        B = gjzpsp(Kk, prime);
        
        %if sum(sum(B-B2))
        %    keyboard;
        %end
        
        gjtime = cputime - gjtime;
        if (gjtime > 0.5)
            fprintf('\t[t:%.2fsec]\t', gjtime);
        end

        M = zeros(size(selrows, 2), size(M, 2));
        M(:, nonzero) = B;
    end

    cols = size(M, 2);
        
    % check conditions for all action matrices
     for i=1:length(amStats)
    
        % test required leading terms
        if max(amStats{i}.amLT) > cols
            continue;
        end
        [r_ones c_ones]=find(M(:, cols-amStats{i}.amLT+1) == 1);
        if size(r_ones,1) < amStats{i}.amLTcnt
            continue;
        end
        rowssum = find(sum(M(r_ones, amStats{i}.zero_el),2) == 1);
        if size(rowssum, 1) < amStats{i}.amLTcnt
            continue;
        end
        
        % all tests passed OK
        res = i;
        return;
    end
    
    res = 0;
end

function [M Mcoefs] = gbs_ReallocMatrix(M, Mcoefs, torows, tocols)

    crows = size(M, 1);
    ccols = size(M, 2);
    
    if crows ~= torows
        
        M = [M; zeros(torows - crows, ccols)];
        Mcoefs = [Mcoefs; zeros(torows - crows, ccols)];
    end
    
    if ccols ~= tocols
               
        % shift coefs
        newcols = zeros(torows, tocols - ccols);
        M = [newcols M];
        Mcoefs = [newcols Mcoefs];
    end
end

function [c] = InvZp(x, p)

    [g,c] = gcd(x,-p);
    c = mod(c,p);
end
