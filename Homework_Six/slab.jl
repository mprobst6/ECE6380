using Plots

function FDTDslab()
    L::Int16 = 2000 # guide length
    A::Int16 = 200 # guide width
    Nz::Int16 = 800 # cells along length use 200->10, 150->13.3, 400-> 5m
    Nx::Int16 = 80 # cells along width
    delz::Float16 = L/Nz # length of cell along z
    delx::Float16 = A/Nx # length of cell along x
    Nsamp::Int16 = round((Nx+1)/2)
    Nzm1::Int16 = Nz-1 # last term

    Sz::Int16 = 1 # cell index of source along length
    Sx::Int16  = 1 # cell index of source along width

    lambda = 300
    BigT = 1.0e-6             # transient excitation of 1 microsecond duration
    mu = pi*4.0e-7            # permeability
    epsilon = 8.854e-12       # permittivity
    c = 2.998e8               # speed of light
    eta=sqrt(mu/epsilon)      # eta
    Hy = zeros(Nx,Nz)      # Hy, size is Nx x Nz
    Ex = zeros(Nx,Nz)      # Ex, size is Nx x Nz
    Ez = zeros(Nx,Nz)      # Ez, size is Nx x Nz
    Jx = zeros(Nx,Nz)      # Jx, size is Nx x Nz
    Hymem = zeros(Nx)     # vector of length Hx

    epsx = ones(Nx,Nz)
    epsz = ones(Nx,Nz)

# going to be 1000 m long, with 100 cells this is 10 m per cell, start at 1/4 the way through for 25 cells to 
# needs to be 300*0.4=120m. Cell width is 2500/Nz
    zstart = convert(UInt16,Nz/2) # starting index for the slab
    zend = convert(UInt16,zstart + 0.4*300/(delz))   # ending index for the slab
    epslb = 5  # epsilon for the slab
    
    # assign the epsilon values for the slab
    for i in zstart:zend
        if i == zstart
            for j in 1:Nx
                epsx[j,i]=(epslb+1)/2
                epsz[j,i]=epslb
            end
        elseif i == zend
            for j in 1:Nx
                epsx[j,i]=(epslb+1)/2
            end
        else
            for j in 1:Nx
                epsx[j,i]=epslb
                epsz[j,i]=epslb
            end
        end
    end
    plot(epsx[convert(UInt16,Nx/2),:]); savefig("epsilon.png")

    # determine time step based on background permittivity
    deltmax = 1/(c*sqrt((1/delz)^2+(1/delx)^2))

    
    delt = 0.9*deltmax

    dtomudz = delt/(mu*delz)      # compute necessary constants
    dtomudx = delt/(mu*delx)      # compute necessary constants
    dtoepdz = delt/(epsilon*delz) # compute necessary constants
    dtoepdx = delt/(epsilon*delx) # compute necessary constants
    beta=1/(2*delz/c/delt + 1)
    alpha=beta*(2*delz/c/delt - 1)

    anim = @animate for kk in 1:1200
        t = kk*delt

        for j in 1:Nx
            Ex[j,1]=sin(2*pi*c*t/lambda)
        end

        # update equations
        for j in 1:Nx
            Hymem[j] = Hy[j,Nzm1] # fill the back row (z=Nz-1, all x values)
        end

        # update equations for Hy in the middle of the mesh
        for i in 1:Nzm1
            for j in 1:Nx
                if j == 1 # top PEC
                    Hy[j,i]=Hy[j,i]-dtomudz*(Ex[j,i+1]-Ex[j,i])+dtomudx*Ez[j+1,i] # no Ez(j-1,i)
                elseif j == Nx
                    Hy[j,i]=Hy[j,i]-dtomudz*(Ex[j,i+1]-Ex[j,i])-dtomudx*Ez[j,i]
                else
                    Hy[j,i]=Hy[j,i]-dtomudz*(Ex[j,i+1]-Ex[j,i])+dtomudx*(Ez[j+1,i]-Ez[j,i])
                end
            end
        end

        # at the absorbing boundary condition
        for j in 1:Nx
            if j == 1
                Hy[j,Nz]=alpha*Hy[j,Nz]+beta*(Hy[j,Nzm1]+Hymem[j])+beta*(Ez[j+1,Nz])/eta
            elseif j == Nx
                Hy[j,Nz]=alpha*Hy[j,Nz]+beta*(Hy[j,Nzm1]+Hymem[j])+beta*(-Ez[j,Nz])/eta
            else
                Hy[j,Nz]=alpha*Hy[j,Nz]+beta*(Hy[j,Nzm1]+Hymem[j])+beta*(Ez[j+1,Nz]-Ez[j,Nz])/eta
            end
        end

        # update Ez
        for i in 2:Nz
            for j in 1:Nx
                Ex[j,i]=Ex[j,i]-dtoepdz*(Hy[j,i]-Hy[j,i-1])/epsx[j,i]-(delt/epsilon)*Jx[j,i]    
            end
        end

        for i in 1:Nz
            for j in 2:Nx
                Ez[j,i]=Ez[j,i]+dtoepdx*(Hy[j,i]-Hy[j-1,i])/epsz[j,i]
            end
        end

        # print("result at time ",string(t),"\n")
        # for i in 1:Nz
        #     z = delz*(i-1)
        #     print("location: ",z,",  field",Ex[Nsamp,i],"\n")
        # end
    
        plot(1:Nz,Ex[2,1:Nz],ylim=(-2,2))
    end
    print("maximum: ",maximum(Ex),"\n")
    print("incident: ",maximum(Ex[1,1:convert(UInt16,Nz/4)]),"\n")
    print("transmission: ",maximum(Ex[1,convert(UInt16,3*Nz/4):Nz]),"\n")
    gif(anim,"anim_pfs50.gif",fps=200)
    


    return
end

FDTDslab()