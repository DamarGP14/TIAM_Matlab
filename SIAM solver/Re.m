function Re=Re(omega,mesh,Im_b)

%Omega must be part of a mesh with the same characteristics
for i=1:length(mesh)
    
    if (mesh(i) == omega)
       
        A=i; %point where epsilon equals omega
        break
    end
end

Re_1_cum=0;

for j=2:length(mesh)-1
    
    if (mesh(j) ~= omega) && (mesh(j+1)~=mesh(A))
        
        Re_1=((Im_b(j+1)-Im_b(A))/(mesh(A)-mesh(j+1))+(Im_b(j)-Im_b(A))/(mesh(A)-mesh(j)))*0.5*(mesh(j+1)-mesh(j));
        
        Re_1_cum=Re_1_cum+Re_1;
        continue
    end
end

if A==1 || A==length(mesh)
    
    Re_2=0;
else

Re_2=0.25*(mesh(A+1)-omega)*((Im_b(A+1)-Im_b(A))/(mesh(A+1)-omega)+(Im_b(A)-Im_b(A-1))/(omega-mesh(A-1)));

end
%No me coge Re_2 porque al mezclar mallas, los valores de las sumas de
%mallas no son iguales que la malla original


%Re=(-Re_1_cum-Re_2+Im_b(A)*log((mesh(length(mesh))-omega)/(omega-mesh(1))))/pi;
%Re=(-Re_1_cum-Re_2+Im_b(A)*log((mesh(length(mesh))-omega)/(omega-mesh(1))))/pi;

if omega==mesh(length(mesh))
Re=(-Re_1_cum-Re_2+Im_b(A)*log((mesh(length(mesh))+(mesh(length(mesh))-mesh((length(mesh)-1)))-omega)/(omega-mesh(1))))/pi;
elseif omega==mesh(1)
Re=(-Re_1_cum-Re_2+Im_b(A)*log((mesh(length(mesh))-omega)/(omega-mesh(1)+mesh(2)-mesh(1))))/pi;
else
Re=(-Re_1_cum-Re_2+Im_b(A)*log((mesh(length(mesh))-omega)/(omega-mesh(1))))/pi;
end
Re=-Re;
end