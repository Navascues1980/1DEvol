classdef Evolution
    %class for Hamiltonians of 1D many-body systems
    
    properties
        %Hamiltonian: Pauli polynomial, involving qubits 1 and 2 (yes, just
        %2-local, TI Hamiltonians)
        h
        %evolution time
        tau
    end
    
    methods
        function obj = Evolution(h, tau)
            %class constructor
            obj.h =h;
            obj.tau=tau;
        end
        
        function num = n(obj)
            %number of qubits
            num=obj.h.n;
        end
        
        
        function [poly_diff,border]=diff_H(obj, poly)
            %differentiates poly with respect to time implicitly, through the Hamiltonian
            %it returns a zero and border = true if the polynomial cannot
            %be differentiated, due to border effects
            
            border=false;
            %find support of polynomial
            qubits=poly.support;
            if isempty(qubits)
                %if the polynomial only contains terms in t, return 0
                poly_diff=Polynomial.zero(poly.n);
            elseif min(qubits)==1 || max(qubits)==obj.n
                %error('I cannot compute the implicit time derivative of the input polynomial, due to border effects.');
                border=true;
                poly_diff=Polynomial.zero(obj.n);
            else
                small=min(qubits);
                large=max(qubits);
                H=Polynomial.zero(obj.n);
                %construct a Hamiltonian to cover the full support of the
                %polynomial
                for k=0:(large-small+1)
                    H=H+obj.h.translate(small-2+k);
                end
                
                %compute the corresponding commutator
                poly_diff=1i*obj.tau*(H*poly-poly*H);
                poly_diff=poly_diff.trim();
                
            end
            
        end
        
        function poly_diff=diff_H_open(obj, poly)
            %differentiates poly with respect to time implicitly, assuming
            %that the system has obj.n particles and open boundary conditions
            
            %construct the full Hamiltonian
            H=0;
            for k=0:obj.n-2
                
                H=H+obj.h.translate(k);
            end
            %compute the corresponding commutator
            poly_diff=1i*obj.tau*(H*poly-poly*H);
            poly_diff=poly_diff.trim();
        end
        
        function [poly_diff, border]=diff(obj, poly)
        %differentiates poly with respect to t, both implicitly and
        %explicitly
        %it returns border=true if the polynomial cannot be differentiated,
        %due to border effects
        [poly_diff_H, border]=obj.diff_H(poly);
        poly_diff=Evolution.diff_time(poly)+poly_diff_H;
            
            
        end
        
        function poly_diff=diff_open(obj, poly)
        %differentiates poly with respect to t, both implicitly and
        %explicitly, assuming that the system is finite (with obj.n)
        %particles and has open boundary conditions
        poly_diff_H=obj.diff_H_open(poly);
        poly_diff=Evolution.diff_time(poly)+poly_diff_H;
            
            
        end
        
    end
    
    methods (Static)
        function poly_diff=diff_time(poly)
            %differentiates poly explicitly with respect to time
            list_mon=[];
            list_coeffs=[];
            for k=1:poly.sparsity
                seq=poly.list_mon(k).seq;
                if seq(1)>0
                    list_mon=[list_mon, Monomial([seq(1)-1; seq(2:length(seq))])];
                    list_coeffs=[list_coeffs, poly.list_coeffs(k)*seq(1)];
                end
            end
            %if there are no terms left, return 0
            if length(list_mon)==0
                poly_diff=Polynomial.zero(poly.n);
            else
                %order the resulting polynomial
                poly_diff=Polynomial(list_coeffs, list_mon);
                poly_diff=poly_diff.sort;
            end
            
        end
        
        
    end
end

