function [Burger,M,y0] = SetBurgers(U_tree,V_tree,R,is_pack)
%This function sets up the blink objects for each leaf. Here 'blink' is set
%to the NonlinOp property.

for i=1:length(U_tree.leafArray)
    
    %Set blink object and blink motion
    Burger{i} = BurgersOp({U_tree.leafArray{i},V_tree.leafArray{i}},R,is_pack);
    
    %Get initial value
    [U,V] = Burger{i}.initial();
    
    U_tree.leafArray{i}.Setvalues(U(:));
    U_tree.leafArray{i}.sample(U(:));
    
    V_tree.leafArray{i}.Setvalues(V(:));
    V_tree.leafArray{i}.sample(V(:));
     
end

if is_pack
    U_tree.pack();
    V_tree.pack();
end
    
    for i=1:length(U_tree.leafArray)      
        
        %Set mass matrix
        M{i} = Burger{i}.massMatrix();
    end

    y0 = [U_tree.Getvalues();V_tree.Getvalues()];
    
end

