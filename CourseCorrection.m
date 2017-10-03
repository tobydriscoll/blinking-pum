function [ output ] = CourseCorrection(rhs,tree,domain)

cs = ASCoursePreconditioner(tree,domain,rhs);

z = rhs - LaplacianForward(tree,domain,cs);

output = cs + ASPreconditioner(tree,domain,z);
end

