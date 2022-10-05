function Im_eig_call=Im_eig_call(omega,mesh,beta,spectral,hybridization)

%Spectral are the indexed spectral functions, total of 16

%Hybridization is a cell array 1D containing which are the baths that a certain
%spectral function is connected and its coefficients.

    position=0;

    for j=1:16 %Rows
    
    if abs(max(hybridization{j})) > 1e-9
        
        position(end+1)=j;
        
    end
    
    end
    
    position(1)=[];
        
    for k=1:length(position)
        
    indexed_spectral{k} = spectral{position(k)};
    indexed_hybridization{k} = hybridization{position(k)};
       
    
    end

Im_eig_call_intern=0;

for i=1:length(indexed_spectral)

Im_eig_call_intern=Im_eig_call_intern + indexed_spectral{i}.*indexed_hybridization{i};

end

Im_eig_call=trapz(mesh,fermi_factor1(omega,mesh-omega,beta).*Im_eig_call_intern)/pi;

end


