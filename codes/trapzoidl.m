function [ trapint ] = trapzoidl( trapi,Nz,Nx,Ny,deldiff )

%  Detailed explanation goes here

point=(trapi(1,1,1)+trapi(1,1,Ny)+trapi(1,Nx,1)+trapi(1,Nx,Ny)+trapi(Nz,1,1)+trapi(Nz,1,Ny)+trapi(Nz,Nx,1)+trapi(Nz,Nx,Ny));

line=0;
for y=2:Ny-1
  line=line+trapi(1,1,y)+trapi(1,Nx,y)+trapi(Nz,1,y)+trapi(Nz,Nx,y);
    
end

for x=2:Nx-1
   line=line+trapi(1,x,1)+trapi(1,x,Ny)+trapi(Nz,x,1)+trapi(Nz,x,Ny);
   
end

for z=2:Nz-1
   line=line+trapi(z,1,1)+trapi(z,1,Ny)+trapi(z,Nx,1)+trapi(z,Nx,Ny);
    
end

surface=0;

for x=2:Nx-1
    for y=2:Ny-1
     
     surface=surface+trapi(1,x,y)+trapi(Nz,x,y);    
    end
    
end

for z=2:Nz-1
    for y=2:Ny-1
        
        surface=surface+trapi(z,1,y)+trapi(z,Nx,y);
    end
end

for z=2:Nz-1
    for x=2:Nx-1
        
        surface=surface+trapi(z,x,1)+trapi(z,x,Ny);
    end
end

volume=0;
for z=2:Nz-1
    for x=2:Nx-1
        for y=2:Ny-1
            
            volume=volume+trapi(z,x,y);
        end
    end
end
trapint=deldiff.^3.*(point+2*line+4*surface+8*volume);

end

