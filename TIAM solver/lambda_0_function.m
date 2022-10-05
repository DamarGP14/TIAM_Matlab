function lambda_0=lambda_0_function(mesh,D,Re_E)

for i=1:16 %Definition of lambda_0

lambda_0_vector(i)=D(i,i)+Re_E{i}(length(mesh)/2);

end

lambda_0=min(lambda_0_vector);

end