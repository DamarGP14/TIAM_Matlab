%Code for SIAM. We compute the DOS for spin up and down. They are plotted
%at the end.
tic
clear all %clear all previous variables

Delta=1;%-Imaginary part of hybridization

E_0=-4; %Single occupied state energy
E_2=0; %Double occupied occupied state energy

U=E_2-2*E_0;%Coulomb repulsion

T=0.008;%Temperature
beta=1/T; %beta

%We define a  mesh

points=2000;%Number of points of the mesh

limit=40;%Limit omega, the mesh goes from -limit to +limit
w_i=0.5;%the value where we pass from a linear distribution of points 
      %to an exponetial one
cero=1e-50; %Our value for 0, as we use logs in the mesh, we can't just put 0

%Left: We need to fix here the "zero" problem and the mesh 

mesh_1=linspace(-limit,-w_i,points/4+1);
mesh_2=-exp(linspace(log(w_i),log(cero),points/4));
mesh_3=exp(linspace(log(cero),log(w_i),points/4));
mesh_4=linspace(w_i,limit,points/4+1);

mesh_1(length(mesh_1))=[];
mesh_4(1)=[];

mesh=[mesh_1 mesh_2 mesh_3 mesh_4];

mesh=mesh';

conv_parameter=10; %Set up of convergence parameter >>0
tolerance=1e-5; %Tolerance. Loop will end when conv_parameter<tolerance


%Set up of Spectral functions as Lorentians centered at 0, width 1
spec_b=0.5*lorentzian(mesh,0,1);
spec_up=0.2*lorentzian(mesh,0,1);
spec_down=0.2*lorentzian(mesh,0,1);
spec_d=0.1*lorentzian(mesh,0,1);


%Set up of arrays. Just to save some computational time
Im_b=zeros(length(mesh),1); 
Im_up=zeros(length(mesh),1);
Im_down=zeros(length(mesh),1);
Im_d=zeros(length(mesh),1);

Re_b=zeros(length(mesh),1);
Re_up=zeros(length(mesh),1);
Re_down=zeros(length(mesh),1);
Re_d=zeros(length(mesh),1);

spec_bn=zeros(length(mesh),1);
spec_upn=zeros(length(mesh),1);
spec_downn=zeros(length(mesh),1);
spec_dn=zeros(length(mesh),1);


iter=0; %number of iterations

lambda_0=0;

fprintf('Computing Spectral functions \n')

fancy_bar=waitbar(0,'Computing Spectral functions. Will it converge? Only time will say. (Check command window) ');

%We start the while loop

while conv_parameter > tolerance 
    
    %Calculation of imaginary parts. These are 'tilda' quantities
    
    for i=1:length(mesh)
        
        omega=mesh(i);
        
        Im_b(i)=trapz(mesh,Im_bf_tilda(omega,mesh,Delta,spec_up,spec_down,beta,lambda_0));  
        Im_up(i)=trapz(mesh,Im_s_upf_tilda(omega,mesh,Delta,spec_b,spec_d,beta,lambda_0));  
        Im_down(i)=trapz(mesh,Im_s_downf_tilda(omega,mesh,Delta,spec_b,spec_d,beta,lambda_0));  
        Im_d(i)=trapz(mesh,Im_df_tilda(omega,mesh,Delta,spec_up,spec_down,beta,lambda_0));  
    
    end
    
    %We compute the Real part of each eigenenergy. Note that we have to
    %multiply by fermi(-w) to get the 'non-tilda' quantity, which is the
    %one that we need to use to get the 'non-tilda' real part.
   
    for i=1:length(mesh)
        
        omega=mesh(i);
        
        Re_b(i)=Re(omega,mesh,Im_b.*fermi(-mesh,beta,lambda_0));
        Re_up(i)=Re(omega,mesh,Im_up.*fermi(-mesh,beta,lambda_0));
        Re_down(i)=Re(omega,mesh,Im_down.*fermi(-mesh,beta,lambda_0));
        Re_d(i)=Re(omega,mesh,Im_d.*fermi(-mesh,beta,lambda_0));
        
    end
    
    %We compute lambda_0 using the method described in Sposetti
    
    lambda_0_vector=[Re_b(length(mesh)/2), E_0+Re_up(length(mesh)/2), E_0+Re_down(length(mesh)/2), E_0+E_0+U+Re_d(length(mesh)/2)];
    
    lambda_0=min(lambda_0_vector); %As Re is fixed (rightnow), lambda_0 wont change and we check that our method to determine it is valid
    
    %We store the value of the current spectral functions to make a
    %comparison
    
    spec_bn=spec_b;
    spec_upn=spec_up;
    spec_downn=spec_down;
    spec_dn=spec_d;
    
    %We compute the spectral functions
        
   spec_b=spec_bf(mesh,lambda_0,Re_b,Im_b.*fermi(-mesh,beta,lambda_0),Im_b);
   spec_up=spec_upf(mesh,lambda_0,Re_up,Im_up.*fermi(-mesh,beta,lambda_0),Im_up,E_0);
   spec_down=spec_downf(mesh,lambda_0,Re_down,Im_down.*fermi(-mesh,beta,lambda_0),Im_down,E_0);
   spec_d=spec_df(mesh,lambda_0,Re_d,Im_d.*fermi(-mesh,beta,lambda_0),Im_d,(E_0+E_0)/2,U);
   
   %We compare the maximum of each spectral function with the maximum of
   %itself in the previous iteration
   
   max_comparison=[abs(max(spec_b)-max(spec_bn)), abs(max(spec_up)-max(spec_upn)), abs(max(spec_down)-max(spec_downn)), abs(max(spec_d)-max(spec_dn))];
   
   %We take the maximum in the previous array and the while loop will check
   %if it fits the convergence criteria
   
   conv_parameter=max(max_comparison);
   
   iter=iter+1; %Count the number of iterations
   
   %Print the conv_parameter to live-check convergence
   fprintf('Convergence parameter = %e \n', conv_parameter)

   
end
    
%Print the number of iterations needed

pause(1)
waitbar(1,fancy_bar,'Converged!!!!')

fprintf('Job done \n')

fprintf('Number of iterations needed %d \n\n',iter)

%We compute the proyection factor Z(1)/Z(0)
    
Z=trapz(mesh,fermi(mesh,beta,lambda_0).*(spec_b+spec_up+spec_down+spec_d));

fprintf('Computing DoS \n')

pause(1)
waitbar(0,fancy_bar,'Computing DoS');
    
%Now we compute the convolution for the DoS
 
diff=zeros(1,length(mesh));
index=zeros(1,length(mesh));
diff_2=zeros(1,length(mesh));
index_2=zeros(1,length(mesh));
product_a=zeros(1,length(mesh));
product_b=zeros(1,length(mesh));
integrand=zeros(1,length(mesh));

diff_d=zeros(1,length(mesh));
index_d=zeros(1,length(mesh));
diff_2_d=zeros(1,length(mesh));
index_2_d=zeros(1,length(mesh));
product_a_d=zeros(1,length(mesh));
product_b_d=zeros(1,length(mesh));
integrand_2=zeros(1,length(mesh));

a=mesh;

for j=1:length(a)
    
b=a-a(j); %We traslate the mesh in each iteration
c=a-a(j);


for i=1:length(a) %a(i)*b(index(i))
    
    [diff(i), index(i)]=min(abs(b-a(i)));
    [diff_d(i), index_d(i)]=min(abs(c-a(i)));
    
end

for i=1:length(c)%b(i)*a(index(i))
    
    [diff_2(i), index_2(i)]=min(abs(a-b(i)));
    [diff_2_d(i), index_2_d(i)]=min(abs(a-c(i)));
    
end



%First multiplication

for i=1:length(a)%product with mesh a
    
    product_a(i)=spec_b(i)*spec_up(index(i));
    product_a_d(i)=spec_up(i)*spec_d(index_d(i));
    
end

for i=1:length(c)%product with mesh b
    
    product_b(i)=spec_up(i)*spec_b(index_2(i));
    product_b_d(i)=spec_d(i)*spec_up(index_2_d(i));

end

%we now mix the two meshes, sort them in ascendent order and apply the same
%ordering to the products

big_mesh=[a' b'];
big_mesh_d=[a' c'];

big_product=[product_a product_b];
big_product_d=[product_a_d product_b_d];

[sorted_mesh, order]=sort(big_mesh);
[sorted_mesh_d, order_d]=sort(big_mesh_d);

sorted_product=big_product(order);
sorted_product_d=big_product_d(order_d);

integrand(j)=trapz(sorted_mesh,fermi_factor1(a(j),sorted_mesh,beta,lambda_0).*sorted_product)/(pi*Z);
integrand_2(j)=trapz(sorted_mesh_d,fermi_factor1(a(j),sorted_mesh,beta,lambda_0).*sorted_product_d)/(pi*Z);

waitbar(j/length(mesh))
end

%Plot the DOS
DoS=integrand+integrand_2;

f13=figure;
plot(mesh,DoS)
[titl,s]=title(['DoS'],['E_0=',num2str(E_0),  ', U=',num2str(U) ', T=', num2str(T)]);
titl.FontSize=16;
grid on
xlim([-20 20])
xl=xlabel('\bf \omega');
xl.FontSize=14;
yl=ylabel('\bf \rho_{\uparrow} (\omega)');
yl.FontSize=14;

close(fancy_bar)
toc
