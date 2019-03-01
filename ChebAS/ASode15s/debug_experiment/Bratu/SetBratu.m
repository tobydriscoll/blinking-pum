function [Operator,M,y0] = SetBratu(U_tree,lambda,pack)
%This function sets up the blink objects for each leaf. Here 'blink' is set
%to the NonlinOp property.

for i=1:length(U_tree.leafArray)
    
    %Set blink object and blink motion
    Operator{i} = BratuOp({U_tree.leafArray{i}},lambda,pack);
    
    %Get initial value
    [U] = Operator{i}.initial();
    
    U_tree.leafArray{i}.Setvalues(U(:));
    U_tree.leafArray{i}.sample(U(:));
   
     
end

    if pack
        U_tree.pack();
    end
    
    for i=1:length(U_tree.leafArray)      
        
        %Set mass matrix
        M{i} = Operator{i}.massMatrix();
    end

    y0 = U_tree.Getvalues();
    
end



