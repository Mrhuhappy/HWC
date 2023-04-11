mrna_exp=importdata('./mrna_exp.csv');
G=importdata('./gene_exp1.txt');
means=importdata('./means.txt');
Matrix=zeros(18067,1204); 
for i=1:18067
    for j=1:1204
        if  mrna_exp.data(i,j)>G.data(i,5)
            Matrix(i,j)=mrna_exp.data(i,j);
        end
        if  mrna_exp.data(i,j)<=G.data(i,5)
            Matrix(i,j)=means(i,1);
        end
    end
end
writematrix(Matrix, './remrna_exp.txt');