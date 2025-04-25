classdef Monomial
    %class to handle monomials of Pauli matrices and time
    
    properties
        %sequence of operators: the first number represents the considered power of
        %time; the rest are numbers from 0 to 3 that represent either the
        %identity (0) or any of the Pauli matrices
        seq
        %number of Pauli matrices
        n
    end
    
    methods
        function obj = Monomial(seq)
            %class constructor
            obj.n = length(seq)-1;
            %make the sequence a column vector
            obj.seq=reshape(seq,length(seq),1);
        end
        
        function deg = degree(obj)
            %computes the degree of the monomial
            deg = length(find(obj.seq>0));
        end
        
        function supp=support(mon)
            %returns the qubits supporting the monomial
            supp=find(mon.seq(2:(mon.n+1))>0);
            supp=sort(supp);
        end
        
        function poly=mon2poly(mon)
            %turns a monomial into a polynomial
            poly=Polynomial([1], [mon]);
            
        end
        
        
        function num=seq2num(obj)
            %associates a number to the monomial sequence
            template=obj.n:(-1):0;
            template=4.^template;
            num=template*obj.seq;
            
        end
        
        
        
        function sol=lt(obj1, obj2)
            %returns "yes", if obj1 is smaller than obj2 in graded
            %lexicographic order
            if not(obj1.n==obj2.n)
               error("I cannot compare monomials for systems with a different number of qubits.");                 
            end
            
            if obj1.degree>obj2.degree
                sol=false;
            elseif obj1.degree<obj2.degree
                sol=true;
            else
                if obj1.seq2num<obj2.seq2num
                    sol=true;
                else
                    sol=false;
                end
                
            end
            
        end
        
        function sol=gt(obj1,obj2)
            %returns "yes", if obj1 is greater than obj2 in graded
            %lexicographic order
            sol=lt(obj2,obj1);
            
        end
        
        function sol=eq(obj1,obj2)
            %returns "yes", if obj1 equals obj2
            sol=isequal(obj1.seq, obj2.seq);            
            
        end
        
        function sol=le(obj1,obj2)
            %returns "yes", if obj1 is smaller than or equal to obj2 in graded
            %lexicographic order
            if (obj1<obj2) || (obj1==obj2)
                sol=true;
            else
                sol=false;
            end
        end    
        
        function sol=ge(obj1,obj2)
            %returns "yes", if obj1 is greater than or equal to obj2 in graded
            %lexicographic order
            if (obj1>obj2) || (obj1==obj2)
                sol=true;
            else
                sol=false;
            end
        end    
        
        function sol=mtimes(obj1, obj2)
           %computes the product of two monomials, resulting in a
           %polynomial
           %first, compute the time contribution
           seq3=zeros(obj1.n+1,1);
           seq3(1)=obj1.seq(1)+obj2.seq(1);
           total_coeff=1;
           for k=2:obj1.n+1
               [coeff,seq3(k)]=combine_paulis(obj1.seq(k),obj2.seq(k));
               total_coeff=total_coeff*coeff;
               
           end
           
           %generate polynomial
           sol=Polynomial(total_coeff,Monomial(seq3));
           
            function [coeff,c]=combine_paulis(a,b)
                %given the identifiers of two Pauli matrices, it returns
                %the result of their product as a number c, denoting the
                %resulting Pauli matrix and a coefficient coeff
                if a==0
                    c=b;
                    coeff=1;
                elseif b==0
                    c=a;
                    coeff=1;
                elseif a==b
                    c=0;
                    coeff=1;
                else
                    c=setdiff([1,2,3],[a,b]);
                    coeff=1i*(-1)^(mod(b-a,3)+1);
                
                end
            end
            
        end
        
        
        function obj_trans=translate(obj,delta)
            %given a monomial obj, it displaces it in the 1D chain by an
            %amount delta
            if delta>0
                if not(isequal(obj.seq((obj.n-delta+2):(obj.n+1)), zeros(delta,1)))
                    error('Impossible translation!');
                end
                seq=[obj.seq(1); zeros(delta,1);obj.seq(2:obj.n+1-delta)];
            elseif delta<0
                if not(isequal(obj.seq(2:(-delta+1)), zeros(-delta,1)))
                    error('Impossible translation!');
                end
                seq=[obj.seq(1); obj.seq((-delta+2):obj.n+1);zeros(-delta,1)];
                
            else
                %if there is no translation, return the same object
                seq=obj.seq;
            end
            obj_trans=Monomial(seq);
            
        end
        
        function [ind,sugg_pos]=search(mon, list_mon)
           %checks if a monomial mon is in an ordered (from low to high) list
           % of monomials list_mon. If that's the case, it returns the 
           %index ind such that list_mon(ind)=mon; otherwise, it returns
           %-1 It also returns a suggested position: if the new monomial is
           %there, then the list will keep being ordered
           if length(list_mon)>0
           %interval of indices where one might find mon
           high_ind=length(list_mon);
           low_ind=1;
           while high_ind-low_ind>1
               %divide the interval [low_ind, high_ind] in two parts,
               %separated by cut_off
               
               cut_off=low_ind+floor((high_ind-low_ind)/2);
               if mon>list_mon(cut_off)
                   low_ind=cut_off;
                   
               elseif mon<list_mon(cut_off)
                   high_ind=cut_off;
               else
                   high_ind=cut_off;
                   low_ind=cut_off;
                   
               end
               
           end
           
           
               
           if list_mon(high_ind)==mon
               ind=high_ind;
               sugg_pos=ind;
           elseif list_mon(low_ind)==mon
               ind=low_ind;
               sugg_pos=ind;
           else
               ind=-1;
               if list_mon(high_ind)<mon
                   sugg_pos=high_ind+1;
               elseif list_mon(low_ind)>mon
                   sugg_pos=low_ind;
               else
                   sugg_pos=low_ind+1;
               end
                   
           end
           
           else
               ind=-1;
               sugg_pos=1;
           end
           
        end
        
        function cadena=spell(obj)
            %returns a text chain in LATEX expressing the monomial
            if obj.seq(1)==0
                cadena='';
            else
                cadena="t";
                if obj.seq(1)>1
                    cadena=cadena + "^" + num2str(obj.seq(1));
                end
            end
            for k=1:obj.n
                %Pauli variable
                var_p=obj.seq(k+1);
                if not(obj.seq(k+1)==0)
                    cadena=cadena+"\sigma^{("+num2str(k)+")}_{"+num2str(var_p)+"}";                    
                end
            end
            %if the result is empty, the monomial must be the identity
            if isempty(cadena)
                cadena='1';
            end
            
        end
        
        function disp(obj)
            %display lists of monomials
            for k=1:length(obj)
                disp(obj(k).spell);
                
            end
            
            
        end
        
        
    end
    
    methods (Static)
       function new_list=sort(list_mon)
            %sorts a list of monomials from lower to higher order
            %length of the list
            L=length(list_mon);
            if (L==1)
                new_list=list_mon;
            else
                list_aux1=Monomial.sort(list_mon(1:ceil(L/2)));
                L_1=length(list_aux1);
                list_aux2=Monomial.sort(list_mon(ceil(L/2)+1:L));
                L_2=length(list_aux2);
                new_list=[];
                ind1=1;
                ind2=1;
                for k=1:L
                    if (ind1<=L_1)&&(ind2<=L_2)
                        if (list_aux1(ind1)<list_aux2(ind2))
                            new_list=[new_list,list_aux1(ind1)];
                            ind1=ind1+1;
                        else
                            new_list=[new_list,list_aux2(ind2)];
                            ind2=ind2+1;
                        end
                    else
                        if (ind1>L_1)
                            new_list=[new_list, list_aux2(ind2:L_2)];
                            break;
                            
                        else
                            new_list=[new_list, list_aux1(ind1:L_1)];
                            break;
                            
                        end
                            
                            
                    end
                end
            
            end
        end
        
        
        function res=isordered(list_mon)
           %checks if a list of monomials is ordered 
           res=true;
           for k=1:length(list_mon)-1
               res=res && (list_mon(k)<list_mon(k+1));
               
           end
            
        end 
        
        function lista_pauli=list_local_pauli_vectors(m, n)
            %lists all m-local Pauli operators in vector form (without the
            %time coordinate)
                        
            
            %first, list all possible local operators of qubits 1,...,m
            first_qubits=cartesian(4,m)-ones(4^m,m);
            L_first=4^m;
            %next, traslate them to fill the whole chain
            lista_pauli=zeros((n-m+1)*L_first,n);
            for k=0:n-m
                lista_pauli(k*L_first+1:(k+1)*L_first,:)=[zeros(L_first,k),first_qubits,zeros(L_first,n-m-k)];
            end
            %eliminate duplicates
            lista_pauli=unique(lista_pauli, 'rows');
            
        end
        
        function lista_fin=add_time(lista_pauli,t_power)
            %given a matrix whose rows correspond to local Pauli monomials (in vector form),
            %it returns a matrix of full monomials in vector form, with the
            %time coordinate going from 0 to t_power
            
            L_pauli=size(lista_pauli,1);
            lista_fin=zeros((t_power+1)*L_pauli,size(lista_pauli,2)+1);
            for k=0:t_power
                lista_fin(k*L_pauli+1:(k+1)*L_pauli,:)=[k*ones(L_pauli,1),lista_pauli];
            end
            
            
            
        end
        
        
        function lista_mon=local_pauli(m,t_power, n)
            %it lists the set of m-local monomials of the Pauli matrices of
            %a system of n qubits (multiplied by t^k, k=0,...,t_power)
            
            if m>n
               error('The neighboring distance exceeds the size of the chain.')
            end
            
            %create a list lista_pauli of all possible m-local Pauli operators
            %in vector form (without the time coordinate)
            lista_pauli=Monomial.list_local_pauli_vectors(m, n);
            L_pauli=size(lista_pauli,1);
            
            %next, add the time variable
            lista_fin=Monomial.add_time(lista_pauli,t_power);
            
            %create a list of monomials
            lista_mon=[];
            for k=1:size(lista_fin,1)
               lista_mon=[lista_mon,Monomial(lista_fin(k,:))];                
            end
            
            %important: don't forget to sort the list!
            lista_mon=Monomial.sort(lista_mon);
        end
        
        function lista_mon=local_pauli_squared(m,t_power, n)
            %it lists the monomials corresponding to local_pauli(m,t:power,
            %n)^2
            
            if m>n
               error('The neighboring distance exceeds the size of the chain.')
            end
            
            %create a list lista_pauli of all possible products of m-local Pauli operators
            %in vector form (without the time coordinate)
            
            %first, we deal with the case of no collisions between the two
            %local operators
                       
            %list all possible local operators of qubits 1,...,2m
            indep_qubits=cartesian(4,2*m)-ones(4^(2*m),2*m);
            L_indep=size(indep_qubits,1);
            
            pauli_nc=[];
            for sep=1:n-2*m
                %sep= separation between the two local operators
                %construct products on qbits 1,...,2m+sep
                local_op=[indep_qubits(:,1:m),zeros(L_indep, sep),indep_qubits(:,(m+1):2*m)];
                
                for trans=0:n-2*m-sep
                    %trans=translation of the operator
                    pauli_nc=[pauli_nc;[zeros(L_indep, trans), local_op,zeros(L_indep, n-trans-2*m-sep)]];
                    
                end
                
            end
            %retain just unique operators
            pauli_nc=unique(pauli_nc, 'rows');
            
            %list all Pauli operators that result from a collision
            tam=min(2*m,n);
            lista_col=Monomial.list_local_pauli_vectors(tam, n);
            
            %join the two lists and eliminate duplicates
            lista_pauli=unique([pauli_nc;lista_col],'rows');
            
            
            %next, add the time variable
            lista_fin=Monomial.add_time(lista_pauli,2*t_power);
            
            %create a list of monomials
            lista_mon=[];
            for k=1:size(lista_fin,1)
               lista_mon=[lista_mon,Monomial(lista_fin(k,:))];                
            end
                        
            %important: don't forget to sort the list!
            lista_mon=Monomial.sort(lista_mon);
            
            
        end
        
        
        function new_list=unique(list_mon)
            %given a list of monomials, it eliminates those that are
            %repeated
            
            new_list=list_mon(1);
            for k=2:length(list_mon)                
                [ind, suggested]=list_mon(k).search(new_list);
                if ind==-1
                    new_list=[new_list(1:suggested-1), list_mon(k), new_list(suggested:end)];
                end
                
            end
            
        end
        
        function all_prods=product_mons(list1,list2)
            %given two lists of monomials, it returns a list containing the
            %support of all the products
            all_prods=[];
            for k=1:length(list1)
                for j=1:length(list2)
                    new_pol=list1(k)*list2(j);
                    all_prods=[all_prods,new_pol.list_mon(1)];                    
                end                
            end
            
            %eliminate repeated elements
            all_prods=Monomial.unique(all_prods);
            
            
        end
        
        function mon=sigma(j,m,n)
            %it creates the monomial corresponding to \sigma_j^{(m)} for
            %a system with n particles
            %first, define the monomial
            
            seq=zeros(n+1,1);
            seq(m+1)=j;
            mon=Monomial(seq);
            
            
        end
        
        function mon=unit(n)
            %it creates the monomial corresponding to 1
            
            seq=zeros(n+1,1);
            mon=Monomial(seq);
            
            
        end
        
        
    end
    
    
    
    
end

