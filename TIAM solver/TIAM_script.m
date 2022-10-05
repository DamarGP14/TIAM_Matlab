tic %Command to measure the time of ejecution
clear all %Delete all previous variables 

syms E1 E2 U1 U2 U12 t J %Symbolic characters for the parameters

%Form of the symbolic Hamiltonian as a matrix. The order of states is the given in the Master Thesis
%Symbolic diagonalization of this Hamiltonian takes computing time so it is
%already given with the code

%Note that if the Hamiltonian is modified, it must be diagonalized

H_sym=[
    0   0	0	0       0	0       0           0               0   0           0           0               0       0               0                   0
    0	E2  0	0       -t	0       0           0               0   0           0           0               0       0               0                   0
    0	0	E2 	0       0	0       0           0               -t   0           0           0               0       0               0                   0
    0   0	0	2*E2+U2	0	0       t           0               0   -t           0           0               0       0               0                   0
    0	-t	0	0       E1	0       0           0	            0   0           0           0               0       0               0                   0
    0	0	0	0       0	E1+E2+U12-J/4 0         0               0   0           0           0               0       0               0                   0
    0	0	0	t       0	0       E1+E2+U12+J/4   0               0   -J/2           0           0               t       0               0                   0
    0	0	0	0       0	0       0           E1+2*E2+U2+2*U12 0  0           0           0               0       t               0                   0
    0	0	-t	0       0	0       0           0               E1  0           0           0               0       0               0                   0
    0   0	0	-t       0	0       -J/2           0               0	E1+E2+U12+J/4   0           0               -t       0               0                   0
    0	0	0	0       0	0       0           0               0   0           E1+E2+U12-J/4   0               0       0               0                   0
    0	0	0	0       0	0       0           0               0   0           0           E1+2*E2+U2+2*U12  0      0               t                   0
    0	0	0	0       0	0       t           0               0   -t           0           0               2*E1+U1  0               0                   0
    0	0	0	0       0	0       0           t               0   0           0           0               0       2*E1+U1+E2+2*U12	0               0
    0	0	0	0       0	0       0           0               0   0           0           t               0       0               2*E1+U1+E2+2*U12    0
    0   0   0   0       0   0       0           0               0   0           0           0               0       0               0               2*E1+U1+2*E2+U2+4*U12
];

%Definition of the parameters of the problem

%The hybridization functions can be modified

D=0.05; %Standard energy unit, usually the bandwidth of the conduction band
hyb_1=@(mesh) 0.3*D; %-Imaginary part of the hybridization of bath 1
hyb_2=@(mesh) 0.3*D; %Bath 2

E1=-4*D; %Energy of the 1st impurity
E2=-4*D; %energy of the 2nd impurity

U1=8*D; %Coulomb repulsion 1st impurity
U2=8*D; %Coulomb repulsion 2nd impurity

t=1e-10*D;%Parameter of spin conserving transfer of electrons between impurities 0.018
J=0*D; %Direct exchange parameter J<0 anti-ferromagnetic J>0 ferromagnetic

U12=0*D;%Coulomb repulsion between impurities

T=0.008*D; %Temperature
beta=1/T; %Inverse temperature

%Mesh definition

for i=1
points=2000;%Number of points of the mesh

limit=10;%Limit omega, the mesh goes from -limit to +limit
w_i=2;%the value where we pass from a linear distribution of points 
      %to an exponetial one
cero=1e-10; %Our value for 0, as we use logs in the mesh, we can't just put 0

%Left: We need to fix here the "zero" problem and the mesh 

mesh_1=linspace(-limit,-w_i,points*0.2+1);
mesh_2=-2*exp(linspace(log(0.5*w_i),log(0.5*cero),points*0.3));
mesh_3=2*exp(linspace(log(0.5*cero),log(0.5*w_i),points*0.3));
mesh_4=linspace(w_i,limit,points*0.2+1);

mesh_1(length(mesh_1))=[];
mesh_4(1)=[];

mesh=[mesh_1 mesh_2 mesh_3 mesh_4];

mesh=mesh';
end

%Initial conditions and definition of Spectral variables

for i=1:16
    spectral{i}=(1/16)*lorentzian(mesh,0,1);
end


%Diagonalization, we diagonalize the symbolic Hamiltonian because it yields
%better results

load('good_diagonalized_t_JSS_minus_var.mat')

%[V_sym,D_sym] = eig(H_sym); %This command does the diagonalization. It is
%suspended because it requires too much time. The previous line loads the
%already diagonalized Hamiltonian

%Substitution of values
V=subs(V_sym);
D=subs(D_sym);

V=double(V);
D=double(D);

V=real(V);
D=real(D);

%We need to normalize the columns of V

for i=1:16
    
    normalization=norm(V(:,i));
    
    V_norm(:,i)=V(:,i)/normalization;
end

%V has the new basis in terms of the old basis (canonical)
%V_inv contains the old basis in terms of the new basis

V_inv=inv(V_norm);

%We delete spureous low values

for i=1:16
    
    normalization=norm(V_inv(:,i));
    
    V_inv_norm(:,i)=V_inv(:,i)/normalization;
end

for i=1:16
    for j=1:16
        if abs(V_norm(i,j))<1e-9
            V_norm(i,j)=0;
        end
    end
end

for i=1:16
    for j=1:16
        if abs(V_inv_norm(i,j))<1e-9
            V_inv_norm(i,j)=0;
        end
    end
end




%MAtrix defining the impurity-bath relations in the original basis
%1 destruction with bath 1; 2 destruction with bath 2
%-1 creation with bath 1; -2 creation with bath 2
%

M_hyb=[
    0 -2 -2 0   -1 0 0 0   -1 0 0 0   0 0 0 0
    2 0 0 -2    0 -1 0 0   0 -1 0 0   0 0 0 0
    2 0 0 -2    0 0 -1 0   0 0 -1 0   0 0 0 0
    0 2 2  0    0 0 0 -1   0 0 0 -1   0 0 0 0
    
    1 0 0 0     0 -2 -2 0  0 0 0 0    -1 0 0 0
    0 1 0 0     2 0 0 -2   0 0 0 0    0 -1 0 0
    0 0 1 0     2 0 0 -2   0 0 0 0    0 0 -1 0
    0 0 0 1     0 2 2  0   0 0 0 0    0 0 0 -1
    
    1 0 0 0     0 0 0 0    0 -2 -2 0   -1 0 0 0
    0 1 0 0     0 0 0 0    2 0 0 -2    0 -1 0 0
    0 0 1 0     0 0 0 0    2 0 0 -2    0 0 -1 0
    0 0 0 1     0 0 0 0    0 2 2 0     0 0 0 -1
    
    0 0 0 0     1 0 0 0    1 0 0 0     0 -2 -2 0
    0 0 0 0     0 1 0 0    0 1 0 0     2 0 0 -2
    0 0 0 0     0 0 1 0    0 0 1 0     2 0 0 -2
    0 0 0 0     0 0 0 1    0 0 0 1     0 2 2 0
    ];

%Now we convert our matrix to the new basis so that we have the terms that
%we need to include in the eigen energy calculation

syms hyb_1_p hyb_1_m hyb_2_p hyb_2_m

for i=1:16
    
    for j=1:16
        
        if M_hyb(i,j) == 1
            
            M_hyb_indexed{i,j}=hyb_1_p;
            
        elseif M_hyb(i,j) == -1
            
            M_hyb_indexed{i,j}=hyb_1_m;
             
        elseif M_hyb(i,j) == 2
            
            M_hyb_indexed{i,j}=hyb_2_p;
            
        elseif M_hyb(i,j) == -2
            
            M_hyb_indexed{i,j}=hyb_2_m;
        else
            M_hyb_indexed{i,j}=0;
        end
    end
end


%M_hyb_new=V_inv*M_hyb_indexed*V;
%Note that we are working with sym operators that we are going to
%substitute with vectors. The operations of the type hyb_1_+*hyb_2_+ will
%not be well defined, if encountered we need to make the substitution 
% * --> .*

for i = 1 : 16
  
  for j = 1 : 16
    cumsum = zeros(length(hyb_1(mesh)),1);
    for k = 1 : 16
      cumsum = cumsum + M_hyb_indexed{i,k}.*V_norm(k,j);
    end
    Product_1{i,j} = cumsum;
  end
end

for i = 1 : 16
  
  for j = 1 : 16
    cumsum = zeros(length(hyb_1(mesh)),1);
    for k = 1 : 16
      cumsum = cumsum + V_inv_norm(i,k).*Product_1{k,j};
    end
    M_hyb_new{i,j} = cumsum;
  end
end

%Now the terms with a - sign are the ones to be hyb(-mesh). So we check
%which are them.
%So to check hyb_1_p we turn every other to zero, substitute with a
%positive number and check the sign, we do this in a new M_hyb_new to not
%erase it
%Couldnt find an easy way to vectorize this process: need correction.

for k=1

hyb_1_p=1;
hyb_1_m=0;
hyb_2_p=0;
hyb_2_m=0;

dummy_M=M_hyb_new;

dummy_M=subs(dummy_M);
dummy_M=double(dummy_M);

for j=1:16 %columns
    
    for i=1:16 %rows
        
        if dummy_M(i,j) > 0
            
            M_hyb_def{i,j}=dummy_M(i,j)*dummy_M(i,j)*hyb_1(mesh);
            M_normalization(i,j)=dummy_M(i,j);
            
        elseif dummy_M(i,j) < 0
            
            M_hyb_def{i,j}=  dummy_M(i,j)*dummy_M(i,j)*hyb_1(-mesh);
            M_normalization(i,j)=-dummy_M(i,j);
            
        else 
            
            M_hyb_def{i,j}=0;
            M_normalization(i,j)=0;
            
        end
    end
end

hyb_1_p=0;
hyb_1_m=1;
hyb_2_p=0;
hyb_2_m=0;

dummy_M=M_hyb_new;

dummy_M=subs(dummy_M);
dummy_M=double(dummy_M);

for j=1:16 %columns
    
    for i=1:16 %rows
        
        if dummy_M(i,j) > 0
            
            M_hyb_def{i,j}=M_hyb_def{i,j}+dummy_M(i,j)*dummy_M(i,j)*hyb_1(-mesh);
            M_normalization(i,j)=M_normalization(i,j)+ dummy_M(i,j);
            
        elseif dummy_M(i,j) < 0
            
            M_hyb_def{i,j}=M_hyb_def{i,j} + dummy_M(i,j)*dummy_M(i,j)*hyb_1(mesh);
            M_normalization(i,j)=M_normalization(i,j)-dummy_M(i,j);
            
        else 
            
            M_hyb_def{i,j}=M_hyb_def{i,j} + 0;
            M_normalization(i,j)=M_normalization(i,j)+0;
            
        end
    end
end

hyb_1_p=0;
hyb_1_m=0;
hyb_2_p=1;
hyb_2_m=0;

dummy_M=M_hyb_new;

dummy_M=subs(dummy_M);
dummy_M=double(dummy_M);

for j=1:16 %columns
    
    for i=1:16 %rows
        
        if dummy_M(i,j) > 0
            
            M_hyb_def{i,j}=M_hyb_def{i,j}+dummy_M(i,j)*dummy_M(i,j)*hyb_1(mesh);
            M_normalization(i,j)=M_normalization(i,j)+ dummy_M(i,j);
            
        elseif dummy_M(i,j) < 0
            
            M_hyb_def{i,j}=M_hyb_def{i,j} + dummy_M(i,j)*dummy_M(i,j)*hyb_1(-mesh);
            M_normalization(i,j)=M_normalization(i,j) - dummy_M(i,j);
            
        else 
            
            M_hyb_def{i,j}=M_hyb_def{i,j} + 0;
            M_normalization(i,j)=M_normalization(i,j)+ 0;
            
        end
    end
end

hyb_1_p=0;
hyb_1_m=0;
hyb_2_p=0;
hyb_2_m=1;

dummy_M=M_hyb_new;

dummy_M=subs(dummy_M);
dummy_M=double(dummy_M);

for j=1:16 %columns
    
    for i=1:16 %rows
        
        if dummy_M(i,j) > 0
            
            M_hyb_def{i,j}=M_hyb_def{i,j}+dummy_M(i,j)*dummy_M(i,j)*hyb_1(-mesh);
            M_normalization(i,j)=M_normalization(i,j)+ dummy_M(i,j);
            
        elseif dummy_M(i,j) < 0
            
            M_hyb_def{i,j}=M_hyb_def{i,j} + dummy_M(i,j)*dummy_M(i,j)*hyb_1(mesh);
            M_normalization(i,j)=M_normalization(i,j)- dummy_M(i,j);
            
        else 
            
            M_hyb_def{i,j}=M_hyb_def{i,j} + 0;
            M_normalization(i,j)=M_normalization(i,j)+0;
            
        end
    end
end

%Now we normalize each site

% for i=1:16 %columns
%     
% %     for i=1:16 %rows
% %         
% %         if M_normalization(i,j) == 0
% %             
% %             M_normalization(i,j)=1;
% %         end
% %         M_hyb_def{i,j}=M_hyb_def{i,j}/M_normalization(i,j);
% %         M_hyb_def{i,j}=abs(M_hyb_def{i,j});
% %     end
% for j=1:16 %rows
%     
% 
%     M_hyb_def{j,i}=(M_hyb_def{j,i}*4)/norm(M_normalization(:,i));
% end
%end

end



%Now we define the functions that we will need to use
%We define an array of functions that take as variables the participating
%pseudo particles for each process

%We check which are the matrix elements that are non zero and add those
%contributions to the functions

for i=1:16 %Columns
    
    for j=1:16 %rows
    M_hyb_new_array{j}=M_hyb_def{j,i};
    end
    
    Im_E_function{i} = @(omega,spectral) Im_eig_call(omega,mesh,beta,spectral,M_hyb_new_array);
    
    Spectral_function{i}=@(D,Re_E,Im_E,Im_E_tilda,lambda_0) spectral_call(mesh,lambda_0,D,Re_E,Im_E,Im_E_tilda);
    
end

conv_parameter=10; %Set up of convergence parameter >>0
tolerance=1e-4; %Tolerance. Loop will end when conv_parameter<tolerance


%Initial Lambda_0 parameter
lambda_0=0;

%We start the while loop

iter=0;

fprintf('Computing Spectral functions \n')

fancy_bar=waitbar(0,'Computing Spectral functions. Will it converge? Only time will say. (Check command window) ');

 while conv_parameter > tolerance
    
    for i=1:16
        
        Im_E_tilda{i}=zeros(length(mesh),1);
        
        for j=1:length(mesh)
            
            omega=mesh(j);
            
            Im_E_tilda{i}(j)=Im_E_function{i}(omega,spectral);
            
        end
    end

    for i=1:16
        
        Im_E{i}=Im_E_tilda{i}.*fermi(-mesh,beta);
    end
    
    for i=1:16
        
        Re_E{i}=zeros(length(mesh),1);
        
        for j=1:length(mesh)
            
            omega=mesh(j);
            
            Re_E{i}(j)=Re(omega,mesh,Im_E{i});
        end
        
    end
    
    lambda_0=lambda_0_function(mesh,D,Re_E);
    
    spectral_old=spectral;
    
    for i=1:16
        
        spectral{i}=Spectral_function{i}(D(i,i),Re_E{i},Im_E{i},Im_E_tilda{i},lambda_0);
        
    end
    
    conv_parameter=conv_check(spectral,spectral_old);
    
    fprintf('Convergence parameter = %e \n', conv_parameter)

    
    iter=iter+1;
    
 end

pause(1)
waitbar(1,fancy_bar,'Converged!!!!')

fprintf('Job done \n')

fprintf('Number of iterations needed %d \n\n',iter)
 
 
%We change the basis of the spectral

for i=1:16
    spectral_new{i}=0;
    for j=1:16
        
    spectral_new{i}=spectral_new{i}+spectral{j}*abs(V_inv(j,i))*abs(V_inv(j,i));
    end
end


%MAtrix to stablish which are the convolutions to me made in the original
%basis
%+ --> A(e+w) --> A(e)
%Spin up 1; Spin down 1; Spin up 2; Spin down 2

M_conv=[
    -1 -1 -1 -1
    -2 -2 -2  1
    -3 -3  1 -2
    -4 -4  2  2
    
    -5  1 -3 -3
    -6  2 -4  3
    -7  3  3 -4
    -8  4  4  4
    
     1 -5 -5 -5
     2 -6 -6  5
     3 -7  5 -6
     4 -8  6  6
     
     5  5 -7 -7
     6  6 -8  7
     7  7  7 -8
     8  8  8  8
     ];
 
 %Computation of the Z_proyection factor
 
 Z=0;
 for i=1:16
     Z=Z+trapz(mesh,fermi(mesh,beta).*spectral_new{i});
 end
 
 %Computation of the convolution
 
 fprintf('Computing DoS \n')

pause(1)
waitbar(0,fancy_bar,'Computing DoS');
 
 diff=zeros(length(mesh),1);
 index=zeros(length(mesh),1);
 
 diff_2=zeros(length(mesh),1);
 index_2=zeros(length(mesh),1);
 
 for i=1:4
     DoS{i}=zeros(length(mesh),1);
 end
 
 
 for i=1:length(mesh)
     
     traslated_mesh=mesh-mesh(i);%We traslate the mesh en each iteration
     
     for j=1:length(mesh)
         
         [diff(j), index(j)]=min(abs(traslated_mesh-mesh(j)));
     end
 
     for j=1:length(mesh)
         
         [diff_2(j), index_2(j)]=min(abs(mesh-traslated_mesh(j)));
     end
     
         big_mesh=[mesh' traslated_mesh'];
        [sorted_mesh,order]=sort(big_mesh);
     
    for k=1:1 %4 columns, 4 desired real DoS. When the spins and impurities are symmetric, it is recommended to do k=1:1 to save time
        
        for h=1:16
            product_1{h,k}=0;
            product_2{h,k}=0;
            for m=1:16
                
                if M_conv(h,k) == -M_conv(m,k)
                    
                    if M_conv(h,k) > 0
         
                        for j=1:length(mesh)
            
                            product_1{h,k}(j)=spectral_new{h}(j)*spectral_new{m}(index(j));
                            product_2{h,k}(j)=spectral_new{h}(index_2(j))*spectral_new{m}(j);
                        end
                    
                    end
                else 
                    product_1{h,k}=product_1{h,k};
                    product_2{h,k}=product_2{h,k};
                end
            end
            if product_1{h,k} ~= 0
            big_product{h,k}=[product_1{h,k} product_2{h,k}];
            sorted_product{h,k}=big_product{h,k}(order);
            
            DoS{k}(i)=DoS{k}(i)+trapz(sorted_mesh,fermi_factor1(mesh(i),sorted_mesh,beta).*sorted_product{h,k})/(pi*Z);
            end
        end
    end
    
    waitbar(i/length(mesh))
 
 end

%plotting of the results, we plot the DOS for the spin up in impurity 1 
 
f1=figure;
grid on
xlim([-1 1])
xl=xlabel('\bf \omega');
xl.FontSize=14;
yl=ylabel('\bf \rho_{1,\uparrow} (\omega)');
yl.FontSize=14;
hold on
plot(mesh,DoS{1},'LineWidth',1.5)

sum_rule=trapz(mesh,DoS{1}); %Computation of the total integrated density, should be 1
 
close(fancy_bar)

toc %Time of ejecution check
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
