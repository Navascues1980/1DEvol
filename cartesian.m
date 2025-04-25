function lista=cartesian(d,n)
%returns a matrix whose rows are elements of the cartesian product
%\{1,...,d\}^n
lista=zeros(d^n,n);
for num=1:d^n
    lista(num,:)=allNumbers(num, d*ones(1,n));
end


end