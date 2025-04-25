classdef Optimization
    %class to model Hamiltonian evolution in the thermodynamic limit
    
    properties
        %Hamiltonian term
        h
        %objective function (Pauli polynomial corresponding to the desired observable)
        objective
    end
    
    methods
        function obj = Optimization(h,objective)
            %class constructor
            obj.h = h;
            obj.objective = objective;
        end
        
        function num=n(obj)
            %number of qubits
            num=obj.h.n;
        end
        
        function M = localizing_matrix(obj, poly, k, k_time)
            %matrix of polynomials of the form \gamma_{ij}=r_i*poly*r_j,
            %where (r_i)_i is a basis of k-local Pauli polynomials with
            %products of t up to degree k_time
            
            %list of monomials labeling the columns of the moment matrix
            lista_mon=Monomial.local_pauli(k,k_time, obj.n);
            
            %compute moment matrix
            %M=zeros(length(lista_mon), length(lista_mon));
            M(length(lista_mon), length(lista_mon))=Polynomial();
            for j=1:length(lista_mon)
                for k=1:length(lista_mon)
                    M(j,k)=(lista_mon(k).mon2poly)*poly*(lista_mon(j).mon2poly);
                    
                end
            end
        end
        
        function var=yalmip_variables(obj,k, k_time)
            %it generates yalmip variables for moment matrices of order k,
            %k_time
            %compute the list of all monomials that can be expressed as
            %products of k-local monomials
            all_prods=Monomial.local_pauli_squared(k,k_time, obj.n);
            %create a list of yalmip variables of the same size and turn
            %the whole thing into a polynomial
            %normalize the variables: y_1=1
            [ind, cosa]=Monomial.unit(obj.n).search(all_prods);
            if (ind==-1)
                error('I did not find the identity polynomial.');
            else
                %vector of variables
                var_vec=sdpvar(length(all_prods),1);
                var_vec(ind)=1;
                var=Polynomial(var_vec,all_prods);
            end
        end
        
        function cons=LTI_constraints(obj,poly_var)
            %given a set of moments in polynomial form, it enforces 
            %LTI constraints when possible
            cons=[];
            mons=poly_var.list_mon;
            vars=poly_var.list_coeffs;
            while length(mons)>1
                mon=mons(1);
                var=vars(1);
                %eliminate monomial and variable
                mons=mons(2:end);
                vars=vars(2:end);
                qubits=mon.support;
                if not(isempty(qubits))
                    %if (min(qubits)>1) && (max(qubits)<obj.n)
                    for k=-min(qubits)+1:(obj.n-max(qubits))
                        [ind, suggested]=mon.translate(k).search(mons);                             
                        if not(ind==-1)
                            %if the translated monomial appears
                            %elsewhere, equal the two variables
                            cons=cons+[cons,var==vars(ind)];
                            %eliminate variable and polynomial from
                            %the list
                            
                            vars(ind)=[];
                            mons(ind)=[];                                 
                        end
                    end
                    
                end
            end
        end
        
        function cons=diff_constraints(obj,tau,poly_var_0,poly_var_t, poly_var_1)
            %generates differential constraints for all polynomials of the
            %free variables
            
            
            evol= Evolution(obj.h, tau);            
            cons=[];
            
            for k=1:poly_var_t.sparsity
                mon=poly_var_t.list_mon(k);
                qubits=mon.support;
                %differentiate polynomial
                [diff_poly, border]=evol.diff(mon.mon2poly);
                %if there are no border effects (and the monomial was not
                %the identity, continue
                if (border==false) && not(mon==Monomial.unit(obj.n))
                    %try to express the new polynomial as a linear
                    %combination of moments in poly_var_t
                    [val_t, mistake]=diff_poly.evaluate(poly_var_t);
                    if (mistake==false)
                        %reconstruct the monomial mon, but eliminate the
                        %time
                        seq_no_time=mon.seq;
                        seq_no_time(1)=0;
                        mon_no_time=Monomial(seq_no_time);
                        %if there is no time component, evaluate mon at
                        %t=0; otherwise, make it zero
                        if mon.seq(1)==0
                            val_0=mon_no_time.mon2poly.evaluate(poly_var_0);
                            poly_0=mon_no_time.mon2poly;
                        else
                            val_0=0;
                            poly_0=Polynomial.zero(obj.n);
                        end
                        val_1=mon_no_time.mon2poly.evaluate(poly_var_1);
                        cons=[cons, val_t==val_1-val_0];
                        
                        
                    end
                end
            end
            
            disp('Differential constraints, generated.');
        end
        
        function cons=diff_constraints_open(obj,tau,poly_var_0,poly_var_t, poly_var_1)
            %generates differential constraints for all polynomials of the
            %free variables under open boundary conditions
            
            
            evol= Evolution(obj.h, tau);            
            cons=[];
            
            for k=1:poly_var_t.sparsity
                mon=poly_var_t.list_mon(k);
                qubits=mon.support;
                %differentiate polynomial
                diff_poly=evol.diff_open(mon.mon2poly);
                %if the monomial was not the identity, continue
                if not(mon==Monomial.unit(obj.n))
                    %try to express the new polynomial as a linear
                    %combination of moments in poly_var_t
                    [val_t, mistake]=diff_poly.evaluate(poly_var_t);
                    if (mistake==false)
                        %reconstruct the monomial mon, but eliminate the
                        %time
                        seq_no_time=mon.seq;
                        seq_no_time(1)=0;
                        mon_no_time=Monomial(seq_no_time);
                        %if there is no time component, evaluate mon at
                        %t=0; otherwise, make it zero
                        if mon.seq(1)==0
                            val_0=mon_no_time.mon2poly.evaluate(poly_var_0);
                            poly_0=mon_no_time.mon2poly;
                        else
                            val_0=0;
                            poly_0=Polynomial.zero(obj.n);
                        end
                        val_1=mon_no_time.mon2poly.evaluate(poly_var_1);
                        cons=[cons, val_t==val_1-val_0];
                        
                        
                    end
                end
            end
            
            disp('Differential constraints, generated.');
        end
        
        function [cons, energy,var]=min_Hamil_cons(obj,k)
            %creates constraints to minimize the Hamiltonian term (with LTI
            %states)
            
            %first, define the (normalized) free variables for the main
            %functional
            var=yalmip_variables(obj,k, 0);
            disp('SDP Variables generated.');
            %generate moment matrix
            M= localizing_matrix(obj, Polynomial.unit(obj.n), k, 0);
            
            %evaluate moment matrix and enforce positive semidefiniteness
            [M_eval, mistake]=Polynomial.evaluate_matrix(M,var);
            cons=[M_eval>=0];
            
            %enforce translation invariance
            cons_LTI=obj.LTI_constraints(var);
            cons=cons + cons_LTI;
            
            %generate objective function
            energy= obj.h.evaluate(var);
            
        end
        
        
        
        function [cons, objective_func, var_0, var_t, var_1]=time_indep_cons_except_TI(obj, k, k_time, initial_state)
            %prepares the time-independent constraints of the optimization
            %problem (except TI). It returns the constraints, the objective
            %function and the real SDP variables related to each of the
            %three linear functionals
            
            %initial_state denotes the initial state of the chain, an
            %element of the computational basis, expressed as a vector of
            %0s and 1s
            
            %verify that the size of the initial state is correct
            if not(length(initial_state)==obj.n)
               error('The size of the initial state does not correspond to the number of qubits of the Optimization instance.');
                
            end
            
            %first, define the (normalized) free variables for the three
            %linear functionals
            var_t=yalmip_variables(obj,k, k_time);
            disp('Variables for \omega_t generated.');
            var_0=yalmip_variables(obj,k, 0);
            disp('Variables for \omega_0 generated.');
            var_1=yalmip_variables(obj,k, 0);
            disp('Variables for \omega_1 generated.');
            %generate moment matrices
            M_t = localizing_matrix(obj, Polynomial.unit(obj.n), k, k_time);
            disp('Polynomial moment matrix 1 generated.');
            M_0= localizing_matrix(obj, Polynomial.unit(obj.n), k, 0);
            disp('Polynomial moment matrix 2 generated.');
            M_1= localizing_matrix(obj, Polynomial.unit(obj.n), k, 0);
            disp('Polynomial moment matrix 3 generated.');
            
            
            %evaluate them and enforce positive semidefiniteness
            [M_t_eval, mistake_t]=Polynomial.evaluate_matrix(M_t,var_t);
            [M_0_eval, mistake_0]=Polynomial.evaluate_matrix(M_0,var_0);
            [M_1_eval, mistake_1]=Polynomial.evaluate_matrix(M_1,var_1);
            if mistake_t || mistake_0 || mistake_1
                error("I could not evaluate the moment matrices.");
            end
            
            cons=[M_t_eval>=0]+[M_0_eval>=0]+[M_1_eval>=0];
            
            disp('Positivity constraints deriving from moment matrices enforced.');
            
            %localizing matrix for T-T^2
            t=Polynomial.time(obj.n);
            %localizing matrix for t>=0
            %Loc_t = localizing_matrix(obj, t, k, k_time-1);
            %localizing matrix for 1-t>=0
            %Loc_t2 = localizing_matrix(obj, 1-t, k, k_time-1);
            Loc_t_t2 = localizing_matrix(obj, t-t^2, k, k_time-1);
            disp('Polynomial localizing matrices generated.');
            
            %evaluate them and impose positive semidefiniteness
            %[Loc_t_eval, mistake_t]=Polynomial.evaluate_matrix(Loc_t,var_t);
            %[Loc_t2_eval, mistake_t2]=Polynomial.evaluate_matrix(Loc_t2,var_t);
            [Loc_t_t2_eval, mistake_t_t2]=Polynomial.evaluate_matrix(Loc_t_t2,var_t);
            %if mistake_t || mistake_t2
            if mistake_t_t2
                error("I could not evaluate the localizing matrices.");
            end
            
            %cons=cons +[Loc_t_eval>=0,Loc_t2_eval>=0];
            cons=cons +[Loc_t_t2_eval>=0];
            
            disp('Positivity constraints deriving from localizing matrices enforced.');
            
            
            %constraints to fix the initial state: |0>^{\otimes n}
            for k=1:obj.n
                %generate sigma_3^k and find the corresponding index
                [ind, suggest]=Monomial.sigma(3,k, obj.n).search(var_0.list_mon);
                
                cons=cons+ [var_0.list_coeffs(ind)==(-1)^(initial_state(k))];
            end
            
            disp('Initial state constraints, fixed.');
            
            %compute objective function
            objective_func=real(obj.objective.evaluate(var_1));
            
            disp('Objective function, computed.');
            
            
        end
        
        function [cons, objective_func, var_0, var_t, var_1]=time_indep_cons(obj, k, k_time)
            %prepares the time-independent constraints of the optimization
            %problem. It returns the constraints, the objective
            %function and the real SDP variables related to each of the
            %three linear functionals
            
            %generate all constraints except TI, assuming that the initial
            %state is |0>^{\otimes \infty}
            [cons, objective_func, var_0, var_t, var_1]=time_indep_cons_except_TI(obj, k, k_time, zeros(1,obj.n));
            %add translation invariance
            cons=cons + [obj.LTI_constraints(var_t),obj.LTI_constraints(var_0),obj.LTI_constraints(var_1)];
        end
        
        function [cons, objective_func, var_0, var_t, var_1]=time_indep_cons_open(obj, k, k_time, initial_state)
            %prepares the time-independent constraints of the optimization
            %problem, assuming a finite system with open boundary conditions. 
            %It returns the constraints, the objective
            %function and the real SDP variables related to each of the
            %three linear functionals
            
            %generate all constraints except TI
            [cons, objective_func, var_0, var_t, var_1]=time_indep_cons_except_TI(obj, k, k_time,initial_state);
            
        end
        
        
        function [cons, objective_func]=pre_optimize(obj, tau,k, k_time)
            %given a target time, it returns all necessary constraints to
            %carry out the optimization
            
            %first, generate the tau-independent constraints
            [cons, objective_func, var_0, var_t, var_1]=time_indep_cons(obj, k, k_time);
            %next, generate the tau-dependent (differential) constraints
            cons_diff=obj.diff_constraints(tau,var_0,var_t, var_1);
            cons=cons+cons_diff;
        
        end
        
        
        
        function [upp, low]=bounds(obj,tau,k, k_time)
            %generates upper and lower bounds for the objective function
            %after a time tau
            
            time_indep_cons_func=@(x,y)obj.time_indep_cons(x,y);
            diff_constraints_func=@(x,y,z,t)obj.diff_constraints(x,y,z,t);
            [upp,low]=Optimization.bounds_funcs(tau,k, k_time,time_indep_cons_func,diff_constraints_func);
            
            
        end
        
        function [upp, low]=bounds_open(obj,tau,k, k_time,initial_state)
            %generates upper and lower bounds for the objective function
            %after a time tau
            
            time_indep_cons_func=@(x,y)obj.time_indep_cons_open(x,y,initial_state);
            diff_constraints_func=@(x,y,z,t)obj.diff_constraints_open(x,y,z,t);
            [upp,low]=Optimization.bounds_funcs(tau,k, k_time,time_indep_cons_func,diff_constraints_func);
            
            
        end
        
        
    end
    
    methods (Static)
        function optim=example(n)
            %generates the local term of a TI Hamiltonian, with n
            %particles, and also an objective function to compute
            X1=Polynomial.sigma(1,1,n);
            Y1=Polynomial.sigma(2,1,n);
            Z1=Polynomial.sigma(3,1,n);
            X2=Polynomial.sigma(1,2,n);
            Y2=Polynomial.sigma(2,2,n);
            Z2=Polynomial.sigma(3,2,n);
            
            %Hamiltonian term
            %h=X1*X2 + Y1*Y2 + Z1*Z2;
            h=X1*Y2;
            
            %magnetization
            objective=Z1;
            
            optim=Optimization(h, objective);
            
        end
        
        function optim=example2(n)
            %generates the local term of a TI Hamiltonian, with n
            %particles, and also an objective function to compute
            X1=Polynomial.sigma(1,1,n);
            Y1=Polynomial.sigma(2,1,n);
            Z1=Polynomial.sigma(3,1,n);
            X2=Polynomial.sigma(1,2,n);
            Y2=Polynomial.sigma(2,2,n);
            Z2=Polynomial.sigma(3,2,n);
            
            %Hamiltonian term
            h=X1*Y2;
            
            %the objective function is the average value of the Heisenberg
            %Hamiltonian
            h2=X1*X2 + Y1*Y2 + Z1*Z2;
            objective=0;
            for k=0:n-2
                objective=objective+h2.translate(k)*(1/(n-1));
            end
            
            
            
            optim=Optimization(h, objective);
            
        end
        
        function optim=example3(n)
            %generates the local term of a TI Hamiltonian, with n
            %particles, and also an objective function to compute
            X1=Polynomial.sigma(1,1,n);
            Y1=Polynomial.sigma(2,1,n);
            Z1=Polynomial.sigma(3,1,n);
            X2=Polynomial.sigma(1,2,n);
            Y2=Polynomial.sigma(2,2,n);
            Z2=Polynomial.sigma(3,2,n);
            
            %Hamiltonian term
            h=(X1*X2 + Y1*Y2 + Z1*Z2)*(1/4);
            %h=X1*Y2;
            
            %the objective function is the magnetization of the last qubit
            objective=Polynomial.sigma(3,n,n);
            
            
            
            optim=Optimization(h, objective);
            
        end
        
        
        function [upp, low]=bounds_funcs(tau,k, k_time,time_indep_cons_func,diff_constrainst_func)
            %generates upper and lower bounds for the objective function
            %after a time tau
            
            %we allow tau to have more than one value
            upp=zeros(length(tau),1);
            low=zeros(length(tau),1);
            
            %first, generate the tau-independent constraints
            tic;
            [cons, objective_func, var_0, var_t, var_1]=time_indep_cons_func(k, k_time);
            tiempo=toc;
            fprintf('Preliminary numerics required %4.1f hours.\n',tiempo/3600)
            for k=1:length(tau)
                
                tic;
                %generate the tau-dependent (differential) constraints
                cons_diff=diff_constrainst_func(tau(k),var_0,var_t, var_1);
            
                %minimize objective function
                optimize(cons+cons_diff,objective_func);
                low(k)=double(objective_func);
            
                %maximize objective function
                optimize(cons+cons_diff,-objective_func);
                upp(k)=double(objective_func);
                tiempo=toc;
                %print the results
                fprintf('This data point required %4.1f hours.\n',tiempo/3600)
                fprintf('For tau=%4.2f, the bounds are [%5.4f, %5.4f].\n',tau(k),low(k),upp(k))
                
                
            end
            
        end
        
        
    end
    
end

