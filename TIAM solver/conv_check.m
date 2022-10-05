function conv_check=conv_check(spectral,spectral_old)

for i=1:length(spectral)
    
    vector(i)=abs(max(spectral{i})-max(spectral_old{i}));
end

conv_check=max(vector);

end