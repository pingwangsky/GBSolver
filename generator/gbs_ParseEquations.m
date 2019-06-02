% Parse input equations
% (GBsolver subroutine)
% by Martin Bujnak, mar2008


function [monomials p p_currdeg symcoefs maxdeg] = gbs_ParseEquations(cfg, eq, known, unknown)

    if ~isfield(cfg, 'eqinstance')

        eqs = cfg.InstanceGenerator(cfg, eq, known, unknown);
    else
        eqs = cfg.eqinstance;
    end
    
    prime = cfg.prime;
    
    monomials=[];
    eqcnt = length(eq);
    usedcoefidx = 1;
    clear symcoefs;
    clear p;
    for i=1:eqcnt

        eq(i) = expand(eq(i));
        
        if eq(i) == 0
            
            warning('empty/zero equation (%d), remove it and run solver again', i);
        end

        p{i}.eq = eq(i);

        % extract monomials
        [p{i}.mons p{i}.deg] = ExtractMonomials(eq(i), unknown);
        p{i}.maxdeg = max( sum(p{i}.deg, 2) );

        % extract coefs
        [subst smonomials] = CreateSubstBlock(p{i}.mons);
        p{i}.coefs = ExtractCoeffs(eq(i), smonomials, subst);
        p{i}.coefs_rc = ExtractCoeffs(eqs(i), smonomials, subst);

        p_currdeg(i) = p{i}.maxdeg;
        p{i}.mons = [p{i}.mons {'1'}];
        p{i}.deg = [p{i}.deg; zeros(1, size(p{i}.deg, 2))];
        p{i}.monsidx = 0;
        p{i}.monscnt = length(p{i}.mons);
        monomials = union(monomials, p{i}.mons);

        % convert (fractions) coeffs to Zp
        zpcoefs = [];
        coefidx = [];
        for j=1:length(p{i}.coefs)

           cmpl = p{i}.coefs_rc(j);
           [a]=sscanf(char(cmpl),'%d/%d');
           if length(a) > 1
               zpcoefs(j) = double(Zp(a(1), prime) / Zp(a(2), prime));
           else
               zpcoefs(j) = double(Zp(a(1), prime));
           end
           
           symcoefs(usedcoefidx) = p{i}.coefs(j);
           coefidx(j) = usedcoefidx;
           usedcoefidx = usedcoefidx+1;
        end

        p{i}.coefs = zpcoefs;
        p{i}.coefsIDX = coefidx;

    end 
    maxdeg = max(p_currdeg);
    p_currdeg = p_currdeg - min(p_currdeg);
end
