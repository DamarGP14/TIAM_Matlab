function L=lorentzian(x,center,width)

L=width./(pi*((x-center).*(x-center)+width*width));

end