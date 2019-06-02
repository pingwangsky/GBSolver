% coeffs extraction alg...
% (c) Martin Bujnak, martinb@dataexpert.sk, aug2006
%
% [coefs] = ExtractCoeffs(eq, terms, subst)
%
% input 
%   eq - symbolic equation
%   terms - vector of terms to extract (example {'x6' ...})
%   subst - substitution terms.. (example {'x^6 = x6' ...})
%
% return
%   coefs - coefficients in order of terms + remaining terms

function [coefs] = ExtractCoeffs(eq, terms, subst)

    s_s_s = 'syms';
    for i_i_i = 1:length(terms)
        s_s_s = [s_s_s ' ' terms{i_i_i}];
    end
    eval(s_s_s);

    exp = sort(expand(eq));
    
    % apply substititions...
    for i_i_i = 1:length(subst)
        exp = expand(maple('algsubs', subst{i_i_i}, exp));
    end
    
    r_s_s = exp;
    for i_i_t = 1:length(terms)

        s_s_s = ['coef = maple(''coeff'',r_s_s,''' terms{i_i_t} ''');'];
        eval(s_s_s)
        coefs(i_i_t) = coef;
        
        % test
        s1_s_s = ['r_s_s = r_s_s - expand(coef*' terms{i_i_t} '); maple(''coeff'',r_s_s,''' terms{i_i_t} ''');'];
        eval(s1_s_s)
        
        if isempty(coef)
            keyboard;
        end
    end
   
    if isempty(r_s_s)
        keyboard;
    end
    
    coefs(i_i_t+1) = r_s_s;
end