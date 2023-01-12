#!/usr/bin/env julia

using GRUtils
using DifferentialEquations
using DiffEqSensitivity, ForwardDiff
using LaTeXStrings

using Gtk
using CSV
using DataFrames
using DelimitedFiles

# What is your Temperature?
print("What's your catalyst?");
catalyst = "Co2Ni" #readline() #!!!
print("Your Catalyst is ", catalyst,"\n");

#----Constants-----
A= 7.46720*7.46720*1e-20;
m_NH3 = 17.03 * 1.66054e-27
m_N2 = 14 * 1.66054e-27
m_H2 = 2 * 1.66054e-27
m_H = 1.66054e-27
R = 8.3144598       # gas constant
N0 = 1e-3
Pa = 1e5
P = 1e5
P0 = 1e5

#---functions------
function nu_gas(Temp, P, A, m)
    """
    Reaction rate constant for adsorption
    
    T           - Temperature in K
    P           - Pressure in Pa
    A           - Surface area in m^2
    m           - Mass of reactant in kg
    """
    kb = 1.38064852E-23 # boltzmann constant
    R = 8.3144598       # gas constant
    P*A /sqrt(2 * pi * m * kb * Temp)
end

function nu_surf(Temp)
    """
    Calculate reaction rate constant for a surface reaction
    
    T       - Temperature in K
    nu      - Pre-exponential factor in s^-1
    Eact    - Activation energy in J/mol
    """
    kb = 1.38064852E-23 # boltzmann constant
    R = 8.3144598 # gas constant
    h = 6.62607004e-34  # planck constant
    (kb/h*Temp) 
end 
#---------------

k_rate= "/home/dholi/MKM/CODE/Co2Ni/dG0_Co2Ni_100000.0Pa_01.dat" #open_dialog("Pick a file", GtkNullContainer(), String[])
print(k_rate);
dG= DataFrame(CSV.File(k_rate; header=1, delim=" "))
T = dG.Temp

Gf1=[]; Gb1=[]; G1= [];
Gf2=[]; Gb2=[]; G2= [];
Gf3=[]; Gb3=[]; G3= [];
Gf4=[]; Gb4=[]; G4= [];
Gf5=[]; Gb5=[]; G5= [];
Gf6=[]; Gb6=[]; G6= [];
Gf7=[]; Gb7=[]; G7= [];
Gf8=[]; Gb8=[]; G8= [];
Gf9=[]; Gb9=[]; G9= [];
Gf10=[]; Gb10=[]; G10= [];

G = []; Gf0 =[]

press = Array([0.0, 0.0,1.0]) .* Pa #
#press = press .* Pa# 1.01325

print("\nYour catalyst is :", catalyst, "\n\n")
print("Temp amm_ads fst_deH sst_deH tst_deH NN NNH N2 N2_H N2_des H2_des\n")

for i in 1:length(T)
    #print("\nHere 1\n\n")
    push!(Gf1,float(dG."dG_ads_1")[i])   ; push!(Gb1, float(dG."dG_des_1")[i])  ; push!(G1, Gf1[i]-Gb1[i]) ;
    push!(Gf2,float(dG."dG_f_2")[i]) ; push!(Gb2, float(dG."dG_b_2")[i])    ; push!(G2, Gf2[i]-Gb2[i]) ;
    push!(Gf3,float(dG."dG_f_3")[i]) ; push!(Gb3, float(dG."dG_b_3")[i])    ; push!(G3, Gf3[i]-Gb3[i]) ;
    push!(Gf4,float(dG."dG_f_4")[i]) ; push!(Gb4, float(dG."dG_b_4")[i])    ; push!(G4, Gf4[i]-Gb4[i]) ;
    push!(Gf5,float(dG."dG_f_5")[i]) ; push!(Gb5, float(dG."dG_b_5")[i])    ; push!(G5, Gf5[i]-Gb5[i]) ;
    push!(Gf6,float(dG."dG_f5")[i])  ; push!(Gb6, float(dG."dG_b5")[i])     ; push!(G6, Gf6[i]-Gb6[i]) ;
    push!(Gf7,float(dG."dG_f_6")[i]) ; push!(Gb7, float(dG."dG_b_6")[i])    ; push!(G7, Gf7[i]-Gb7[i]) ;
    push!(Gf8,float(dG."dG_f6")[i])  ; push!(Gb8, float(dG."dG_b6")[i])     ; push!(G8, Gf8[i]-Gb8[i]) ;
    push!(Gf9,float(dG."dG_f_7")[i]) ; push!(Gb9, float(dG."dG_b_7")[i])    ; push!(G9, Gf9[i]-Gb9[i]) ;
    push!(Gf10,float(dG."dG_f_8")[i]); push!(Gb10, float(dG."dG_b_8")[i])   ; push!(G10, Gf10[i]-Gb10[i]);

    push!(G, Array([G1[i],G2[i], G3[i], G4[i], G5[i], G6[i], G7[i], G8[i], G9[i], G10[i]]));
    push!(Gf0, Array([Gf1[i],Gf2[i], Gf3[i], Gf4[i], Gf5[i], Gf6[i], Gf7[i], Gf8[i], Gf9[i], Gf10[i]]));
    
    #print("\nHere 2\n\n")
    
    function calc_rate2(y, Gf)
        global G
        Gr = Gf .- G[i]
        
        P_N2 = press[1]
        P_H2 = press[2]
        P_NH3= press[3]

        Th_NH3 = y[1]
        Th_NH2 = y[2]
        Th_NH  = y[3]
        Th_N   = y[4]
        Th_NN  = y[5]
        Th_NNH = y[6]
        Th_N2  = y[7]
        Th_H   = y[8]
        Th     = y[9]

        #P_N2 = y[10]
        #P_H2 = y[11]
        #P_NH3 = y[12]

        coef=Array([3.0, 3.0, 3.0, 2.0, 0.5, 1.0, 0.5, 1.0, 1.5, 4.5])
        
        nu_frw = Array([nu_gas(T[i], P0 , A, m_NH3),
                    nu_surf(T[i]) ,
                    nu_surf(T[i]),
                    nu_surf(T[i]),
                    nu_surf(T[i]),
                    nu_surf(T[i]),
                    nu_surf(T[i]),
                    nu_surf(T[i]),
                    nu_gas(T[i], P0, A, m_N2),
                    nu_gas(T[i], P0, A, m_H2)])
        
        nu_rev = Array([nu_gas(T[i], P0, A, m_NH3),
                    nu_surf(T[i]) ,
                    nu_surf(T[i]),
                    nu_surf(T[i]),
                    nu_surf(T[i]),
                    nu_surf(T[i]),
                    nu_surf(T[i]),
                    nu_surf(T[i]),
                    nu_gas(T[i], P0, A, m_N2),
                    nu_surf(T[i])])
        
        fwd_factor = Array([P_NH3/P0* Th,
                    Th_NH3 * Th* N0,
                    Th_NH2 * Th* N0,
                    Th_NH * Th* N0,
                    Th_N^2* N0,
                    Th_NH *Th_N* N0,
                    Th_NN* N0,
                    Th_NNH* N0,
                    Th_N2,
                    Th_H^2* N0])
        rev_factor = Array([Th_NH3,
                    Th_NH2 * Th_H* N0,
                    Th_NH * Th_H* N0,
                    Th_N * Th_H* N0,
                    Th_NN* N0,
                    Th_NNH* N0,
                    Th_N2 * Th* N0,
                    Th_N2 * Th_H* N0,
                    P_N2/P0* Th,
                    Th^2 * P_H2/P0* N0])
        r = (coef .* nu_frw .* exp.(-Gf./(R.*T[i])) .* fwd_factor .- coef .* nu_rev .* exp.(-Gr./(R.*T[i])).* rev_factor) #!!!
        #print(nu_frw .* exp.(-Gf./(R.*T[i])))
        return r
    end
    
    #print("\nHere 3\n\n")
    
    function calc_overall_rate2(y, Gf)
        global G
        Gr = Gf .- G[i]  #!!!
        #press = Array([0.4,1.0,0.16]) #0.4,1.2,0.16 #0.0, 0.0, 1.0
        #press = press .* Pa# 1.01325
        #P0 = sum(press)
        # N2 formation
        r1 = 1.5*nu_gas(T[i], P0, A, m_N2)*(exp(-Gf[9]/(R*T[i])) * y[7] - exp(-Gr[9]/(R*T[i])) * press[1]/P0 * y[9]) # y[10]=press[1]
        return r1
    end
    #print("\nHere 4\n\n")
    
    function state_eqn2(y, Gf, t)
        r = calc_rate2(y, Gf)
        eq1 = (r[1]- r[2])              # dθ(NH3])/dt
        eq2 = (r[2]- r[3])              # dθ(NH2])/dt
        eq3 = (r[3]- r[4]- r[6])         # dθ(NH])/dt
        eq4 = (r[4]- 2*r[5]- r[6])       # dθ(N])/dt 2*r[5
        eq5 = (r[5]- r[7])              # dθ(N]-N])/dt
        eq6 = (r[6]- r[8])              # dθ(N]-NH])/dt
        eq7 = (r[7]+ r[8]- r[9])         # dθ(N2])/dt
        eq8 = (r[2]+ r[3]+ r[4]+ r[8]- 2*r[10])   # dθ(H])/dt
        eq9 = - (eq1+ eq2+ eq3+ eq4+ eq5+ eq6+ eq7+ eq8) # dθ/dt
        #eq10= r[7]
        #eq11= r[8]
        #eq12= -r[1]
        Array([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9]) #, eq10, eq11, eq12
    end
    #print("\nHere 5\n\n")
    #time_end = 1e6 #!!!
    
    if T[i] < 250
        time_end = 1e16 #!!!
        time_in = 1e12
        nt = Int(1e6) #!!!
        
    elseif T[i] in 260:320    
        time_end = 1e14 #!!!
        time_in = 1e10
        nt = Int(1e6) #!!!
        
    elseif T[i] == 340    
        time_end = 1e12 #!!!
        time_in = 1e8
        nt = Int(1e6) #!!!

    elseif T[i] in 360:580    
        time_end = 1e11 #!!!
        time_in = 1e6
        nt = Int(1e6) #!!!
        
    else #if T > 550
        time_end = 1e6 #!!!
        time_in = 10
        nt = Int(1e6) #!!!

    end 

    tspan = (0., time_end)
    th_NH3_0 =  0.0
    theta = 1-th_NH3_0  # Initial free sites
    y0 = [th_NH3_0; 0.0; 0.00; 0.00; 0.0; 0.0; 0.0; 0.0; theta ]#   ; 0; 0; P]
    prob2 = ODEProblem(state_eqn2, y0, tspan, Gf0[i]) ; #, dtmin =1e-16
    
    function rate_wrapper2(p)
        p = (p)
        #nt = Int(5e6) #!!!
        _prob = remake(prob2, p=p)
        theta_sol = Array(solve(_prob, Rodas4P(), saveat=LinRange(time_in, time_end,nt),
                sensealg=SensitivityADPassThrough(), reltol=1e-6,abstol=1e-10, dtmin = 1e-20))' #
                #ForwardDiffSensitivity     #Rodas4P  
                #QuadratureAdjoint          #KenCarp4
                #SensitivityADPassThrough   #ImplicitHairerWannerExtrapolation
        p_matrix = reshape(p, 1, size(p, 1)) # [1, np]
        p_repeat = repeat(p_matrix, nt, 1) # [nt, np]
        r = Array{Real, 1}(undef, nt)
        for i in 1:nt
            r[i] = calc_overall_rate2(theta_sol[i, :], p_repeat[i, :])
        end
        return (log.(r))
    end
    
    drdp_N2 = ForwardDiff.jacobian(rate_wrapper2, Gf0[i] ); 
    
    #print("\nHere 8\n")
    
    DrC = round.((drdp_N2[end, :])*(-R*T[i]); digits=3)
    
    #print("\nHere 9\n")

    addT = pushfirst!((DrC), T[i]); #add Temperature to the the 1st column
    nextrow = push!(string.(addT),"\n"); #add next row to line
    print(join(nextrow," "))
end
