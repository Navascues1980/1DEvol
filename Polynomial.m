classdef Polynomial
    %class for polynomials of the Pauli matrices and time
    
    properties
        %list of monomials
        list_mon
        %list of coefficients
        list_coeffs
    end
    
    methods
        function obj = Polynomial(list_coeffs,list_mon)
            %class constructor
            if nargin==0
                %empty constructor to handle arrays
                obj=Polynomial.zero(1);
            else
                if length(list_coeffs)==length(list_mon)
                    obj.list_coeffs=reshape(list_coeffs, 1,length(list_coeffs));
                    %obj.list_coeffs=list_coeffs;
                    obj.list_mon=reshape(list_mon, 1, length(list_mon));
                    %obj.list_mon=list_mon;
                else
                    error('The lengths of the list of monomials and coefficients must coincide!');
                end
            end
        end
        
        
        function supp=sparsity(obj)
            %returns the number of non-zero terms
            supp=length(obj.list_mon);
        end
        
        function supp=support(poly)
            %returns the qubits supporting the polynomial
            supp=[];
            for k=1:poly.sparsity
                supp=union(supp, poly.list_mon(k).support);
                
            end
            supp=sort(supp);
        end
        
        function num=n(obj)
           %number of qubits
           num=obj.list_mon(1).n;
            
        end
        
        function new_obj=trim(obj)
           %it eliminates monomials with zero coefficient
           mon=[];
           coeffs=[];
           for k=1:obj.sparsity
               if not(obj.list_coeffs(k)==0)
                  mon=[mon,obj.list_mon(k)];
                  coeffs=[coeffs,obj.list_coeffs(k)];
                  
               end
               
           end
           
           if isempty(mon)
              %if there are no coefficients left, return the zero
              %polynomial
              new_obj=Polynomial.zero(obj.n);
           else
              new_obj=Polynomial(coeffs,mon);
           end
            
        end
        
        function res=iszero(obj)
            %returns true if the polynomial is zero
            if isequal(obj.trim.list_coeffs,0)
                res=true;
            else
                res=false;
            end
            
        end
        
        function poly_sorted=sort(poly)
            %orders the monomials in the polynomial
            list_mon=Monomial.sort(poly.list_mon);
            list_coeffs=[];
            for k=1:poly.sparsity
                [ind, ~]=list_mon(k).search(list_mon);
                
                list_coeffs=[list_coeffs,poly.list_coeffs(ind)];
                
            end
            poly_sorted=Polynomial(list_coeffs,list_mon);
            
        end
        
        
        function obj3=plus_poly(obj1,obj2)
            %returns the sum of two actual polynomials
            mon3=obj2.list_mon;
            coeffs3=obj2.list_coeffs;
            for k=1:obj1.sparsity
                %check if the monomial is repeated                
                [ind,sugg_pos]=obj1.list_mon(k).search(mon3);
                if ind==-1
                    %if the monomial is not there, add it in the right
                    %position
                    mon3=[mon3(1:(sugg_pos-1)),obj1.list_mon(k),mon3(sugg_pos:end)];
                    coeffs3=[coeffs3(1:(sugg_pos-1)),obj1.list_coeffs(k),coeffs3(sugg_pos:end)];
                   
                else
                    %if it is there, update the polynomial
                    coeffs3(ind)=coeffs3(ind)+obj1.list_coeffs(k);
                    
                    
                end
            end
            obj3=Polynomial(coeffs3,mon3);
            obj3=obj3.trim;
            
        end
        
        function obj3=plus(obj1,obj2)
            %returns the sum of two polynomials, or a polynomial and a
            %number
            if isnumeric(obj1)
                %if the first argument is a number, make it a polynomial
                %and add them both
                obj3=obj1*Polynomial.unit(obj2.n)+obj2;
                obj3=obj3.trim;                
            elseif isnumeric(obj2)
                %the same, in the opposite case
                obj3=obj1+obj2*Polynomial.unit(obj1.n);
                obj3=obj3.trim;                
            else
                %if not, call the function for polynomial products
                obj3=plus_poly(obj1,obj2);
            end
            
            %trim the final result
            obj3=obj3.trim;
            
        end
        
        
        function obj3 = mtimes(obj1,obj2)
            %multiplication of two polynomials
            if isnumeric(obj1)
                %if the first argument is a number, just multiply the
                %coefficients of the second argument by obj1 and trim
                obj3=Polynomial(obj1*obj2.list_coeffs,obj2.list_mon);
                obj3=obj3.trim;
            elseif isnumeric(obj2)
                %the same, in the opposite case
                obj3=Polynomial(obj2*obj1.list_coeffs,obj1.list_mon);
                obj3=obj3.trim;
            elseif isa(obj2,'Monomial')
                %if the second argument is a monomial, turn it into a
                %polynomial
                obj3=obj1*obj2.mon2poly;
            elseif isa(obj1,'Monomial')
                %if the first argument is a monomial, turn it into a
                %polynomial
                obj3=obj1.mon2poly*obj2;
            else
                obj3=Polynomial.zero(obj1.n);
                for k=1:obj1.sparsity
                    for j=1:obj2.sparsity
                        coeff1=obj1.list_coeffs(k);
                        coeff2=obj2.list_coeffs(j);
                        mon=obj1.list_mon(k)*obj2.list_mon(j);
                        obj3=obj3+(coeff1*coeff2)*mon;
                    end
                end
            end
            
                
        end
        
        function obj3=minus(obj1,obj2)
        %difference between two polynomials
            obj3=obj1 + (-1)*obj2;
        end
        
        function poly2=mpower(poly, m)
           %computes poly^m
           if (mod(m,1)==0) && (m>=0)
               poly2=Polynomial.unit(poly.n);
               for k=1:m
                   poly2=poly2*poly;
               end
           else
              error('The exponent of a polynomial must be a non-negative integer!');
           end
        end
        
        
        function obj_trans=translate(obj, delta)
            %translates a polynomial by a given amount delta (if the size of the
            %chain allows)
            list_mon=[];
            for k=1:obj.sparsity
                list_mon=[list_mon, obj.list_mon(k).translate(delta)];
            end
            %the final polynomial is already ordered if obj already was
            obj_trans=Polynomial(obj.list_coeffs,list_mon);
            
        end
        
        function cadena=spell(obj)
           %prints the polynomial in LATEX form
           if obj.iszero==true
               cadena="0";
           else
               cadena="";
           
                for k=1:obj.sparsity
               
                    coeff=obj.list_coeffs(k);
                    %add a plus sign if either the real part of the number is
                    %positive or the real part is equal to zero and the
                    %imaginary part is positive
                            
                    if real(coeff)*imag(coeff)==0
                        if imag(coeff)==0
                            if real(coeff)>0
                                cadena=cadena+"+" + num2str_m1(real(coeff));
                            else
                                cadena=cadena + num2str_m1(real(coeff));
                            end
                        end
                        if real(coeff)==0
                                if imag(coeff)>0
                                    cadena=cadena+"+" + num2str_m1(imag(coeff))+ "i";
                                else
                                    cadena=cadena + num2str_m1(imag(coeff))+ "i";
                                end
                        end
                    else
                        cadena=cadena + "(" + num2str(coeff) + ")";
                   
                    end
                    %add the monomial
                    if not((sum(obj.list_mon(k).seq)==0))||(coeff==1)
                        cadena=cadena+obj.list_mon(k).spell;
                    end
                end
           end
        
           
           function cade=num2str_m1(num)
                %given a number, it returns its value as a strong as long
                %as it is not 1
                if num == 1
                    cade="";                    
                elseif num ==-1
                    cade="-";
                else
                    cade=num2str(num);
                end
           end
            
        end
        
        function [sol, mistake]=evaluate(obj,y_vector)
            %given a polynomial obj, it evaluates it by replacing the
            %monomials in obj by values in list_values (which is coded as a
            %Polynomial) It returns "mistake", if y_vector cannot express
            %the polynomial
            sol=0;
            mistake=false;
            for k=1:obj.sparsity
                [ind_mon,~]=obj.list_mon(k).search(y_vector.list_mon);
                if ind_mon==-1
                    mistake=true;
                else
                    sol=sol+obj.list_coeffs(k)*y_vector.list_coeffs(ind_mon);
                end
                
            end
        end
        
            
        function disp(obj)
            %display matrices of polynomials
            
            for j=1:size(obj,1)
                cadena="";
                for k=1:size(obj,2)
                    palabra=obj(j,k).spell;
                    cadena=cadena+ palabra+ blanks(50-strlength(palabra));
                    
                end
                disp(cadena)
            end
            
        end
        
            
    end
    
    
    methods(Static)
        
        function obj=zero(n)
            %creates the zero polynomial
            obj=Polynomial([0],[Monomial(zeros(n+1,1))]); 
            
        end
        
        function obj=unit(n)
            %creates the polynomial 1
            obj=Polynomial([1],[Monomial(zeros(n+1,1))]);
            
        end
        
        function obj=time(n)
            %creates the polynomial 1
            obj=Polynomial([1],[Monomial( [1;zeros(n,1)] )]);
            
        end
        
        function obj=sigma(j,m,n)
            %it creates the polynomial corresponding to \sigma_j^{(m)} for
            %a system with n particles
            %first, define the monomial
            
            seq=zeros(n+1,1);
            seq(m+1)=j;
            mon=Monomial(seq);
            %construct the polynomial
            obj=Polynomial([1],mon);
            
            
        end
        
        function [M_eval, mistake]=evaluate_matrix(M, poly)
            %evaluates a matrix of polynomials: returns mistake=true if
            %poly does not allow to evaluate all polynomials
            %compatible with YALMIP coefficients of poly
            L=size(M,1);
            M_eval=sdpvar(size(M,1), size(M,2),'full');
            mistake=false;
            for j=1:size(M,1)
                for k=1:size(M,2)
                    [aux, mistake2]=M(j,k).evaluate(poly);
                    M_eval(j,k)=aux;
                    mistake=or(mistake, mistake2);
                end
            end
            
            
            
        end
        
        
    end
        
end

