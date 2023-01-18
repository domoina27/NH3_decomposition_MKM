#!/usr/bin/env julia

"""
@author: Domoina Holiharimanana
date: 11/09/2022
Project name: NH3 decomposition over abundant earth metals
Southern University of Illinois Carbondale
"""
# =============================================================================
# Import modules/packages
# =============================================================================
using GRUtils
using DifferentialEquations
using DiffEqSensitivity
using ForwardDiff
using OrdinaryDiffEq
using LaTeXStrings

#using Gtk
using CSV
using DataFrames
using DelimitedFiles
using ParameterizedFunctions

# =============================================================================
# Secify your catalyst
# =============================================================================
#println("WHAT is your catalyst??? \nChoose between: Co, Co2Ni, CoNi2 and Ni")
catalyst = "Co2Ni"   #readline()
println("Your catalyst is: ",catalyst)
simulation = "Batch"
# =============================================================================
# Naming each elementary steps
# =============================================================================
Beta = 1 # K/s
P = 1e5 # pressure in Pa
P0= 1e5 
A= 7.46720*7.46720*1e-20; #Surface in m^2
Pascal = 1e5  #Unit conversion
Joule= 96.485e3  #1eV = 96.485e3 Joule/mol
N0 = 1#0.001 #available active sites
t_end = 5e5 #s
t_endstr= "5e5"

# Constants
kb = 1.38064852E-23 # boltzmann constant
R = 8.3144598 # gas constant
h = 6.62607004e-34  # planck constant

# =============================================================================
# Rate constant
# =============================================================================
function k_surf(T, Eact)
    """
    Calculate reaction rate constant for a surface reaction
    
    T       - Temperature in K
    nu      - Pre-exponential factor in s^-1
    Eact    - Activation energy in J/mol
    """
    kb = 1.38064852E-23 # boltzmann constant
    R = 8.3144598 # gas constant
    h = 6.62607004e-34  # planck constant
    (kb/h*T)* exp((-Eact) / (R * T)) 
end 

function k_ads(T, P, A, m)
    """
    Reaction rate constant for adsorption
    
    T           - Temperature in K
    P           - Pressure in Pa
    A           - Surface area in m^2
    m           - Mass of reactant in kg
    """
    kb = 1.38064852E-23 # boltzmann constant
    R = 8.3144598       # gas constant
    P*A / sqrt(2 * pi * m * kb * T)*exp(-0.1*Joule/(R*T))
end

function k_des(T, A, m, sigma, theta_rot, Edes)
    """
    Reaction rate constant for desorption
    
    T           - Temperature in K
    A           - Surface area in m^2
    m           - Mass of reactant in kg
    sigma       - Symmetry number
    theta_rot   - Rotational temperature in K
    Edes        - Desorption energy in J/mol
    """
    kb = 1.38064852e-23 # boltzmann constant
    h = 6.62607004e-34  # planck constant
    R = 8.3144598       # gas constant
    I= 1856

    kb * T / h^3 * A * (2 * pi * m * kb) /(sigma * theta_rot) * exp(-Edes / (R*T))
end
# Vibrational partition functtion and derivative of the partition function

function Qvib(T, frequency)  
    R = 8.3144598       # kg⋅m2⋅s−2⋅K−1⋅mol−1   Ideal gas constant
    kb = 1.38064852e-23 #J⋅K−1                  Boltzmann Constant
    h = 6.62607004e-34  # m2 kg / s             planck constant
    c = 299792458       # m / s                 speed of light
    qvib= Float64[]
    for i in 1:(length(frequency))
        function q_vib(frequency)  
            den = (kb)/(h*c)*0.01 #!!!!convert m to cm #h*c/(kb*T)
            1/(1-exp(-frequency/(den*T)))
        end 
        push!(qvib, q_vib(frequency[i]))
        # println((qvib))
    end 
    prod(qvib)
end 

function Qtran(T,m) #!!!
    R = 8.3144598       # kg⋅m2⋅s−2⋅K−1⋅mol−1   Ideal gas constant
    kb = 1.38064852e-23 #J⋅K−1                  Boltzmann Constant
    h = 6.62607004e-34  # m2 kg / s             planck constant
    c = 299792458       # m / s  
    P = 1e5; # Pressure in Pa
    #V = 6.9032426e-28
    V = kb*T/P;    # Pressure in Pa
    #V= 15*7*7*1e-30 #kb*T/(Pascal)   
    #V = 0.001 #m^3 1L
    (sqrt(2*pi*m*kb*T)/h)^(3)*V#;*1000*1e-6
end 

function Qtran2D(T,m)
    R = 8.3144598       # kg⋅m2⋅s−2⋅K−1⋅mol−1   Ideal gas constant
    kb = 1.38064852e-23 #J⋅K−1                  Boltzmann Constant
    h = 6.62607004e-34  # m2 kg / s             planck constant
    c = 299792458       # m / s                 speed of light
    Acat= 7.46720*7.46720*1e-20 #
    (2*pi*m*kb*T)/h^2 *A #*1000*1e-6
end 

function Qrot(T,sigma,B)
    """
    sigma is the symetry number
    B is the rotational constnt
    """
    kb = 1.38064852e-23 #J⋅K−1                  Boltzmann Constant
    h = 6.62607004e-34  # m2 kg / s             planck constant
    c = 299792458       # m / s                 speed of light
    cm_to_m = 1e-2      #conversion unit from cm-1 to m-1
    q_rotation = 1/sigma/c *kb*T/(h*B)*cm_to_m
end 

function QrotNL(T,sigma,ABC)
    """
    sigma is the symetry number
    B is the rotational constnt
    """
    kb = 1.38064852e-23 #J⋅K−1                  Boltzmann Constant
    h = 6.62607004e-34  # m2 kg / s             planck constant
    c = 299792458       # m / s                 speed of light
    cm_to_m = 1e-2      #conversion unit from cm-1 to m-1
    q_rotation = 1/sigma *(kb*T/h/c*cm_to_m)^(3/2)*sqrt(pi/ABC)
end 

function ZPE(frequency)
    # zpe = zeros(len(frequency))
    h = 6.62607004e-34  # m2 kg / s             planck constant
    c = 299792458       # m / s 
    NA = 6.0221409e23   #Avogadro"s number
    sum(frequency)*100*299792458*6.62607004e-34*6.0221409e23/2
end 

function ZPE_eV(frequency)
    # zpe = zeros(len(frequency))
    h = 6.62607004e-34  # m2 kg / s             planck constant
    c = 299792458       # m / s 
    NA = 6.0221409e23   #Avogadro"s number
    sum(frequency)*100*299792458*6.62607004e-34*6.0221409e23/2 * 0.00001036427230133138
end 

function EXP(Eact,T)
    R = 8.3144598 
    exp(-Eact / (R * T)) 
end
# ---------------------------------------------------------------------------------------------------

# =============================================================================
# Details for each specific catalysts
# =============================================================================
if (catalyst=="Ni")
    # frequency
    frequency = DataFrame(CSV.File("FRQ/Ni_batch_corrected.csv";  header=1, delim=","))
    #Hydrogen coverage
    alpha_H= -0.03*96.48530749926e3; #intercept_H = 0; #.6622*2*96.48530749926e3 #convert eV to KJ/mol
    #Nitrogen coverage
    alpha_N= -2.13*96.48530749926e3; intercept_N = 0; #1.1521*96.48530749926e3*2 #convert eV to KJ/mol
    #activation barriers
    Eng = (59218.73549, 107248.73, 108809.8592, 52317.63, 87250.3205, 120091.32, 93430.27295, 161641.78, 139063.8154, 164229.31, 36574.94024, 91.20, 61273.80983, 29321.46, 103389.8058, 32915.81, 97068.27); #with ZPE
    #Eng = (72512.78, 125540.90, 113172.2678, 71470.35, 96038.39527, 136039.66, 100315.8504, 167340.88, 139748.1857, 169334.25 , 38673.20 , 15584.5, 67225.0236, 50192.97   , 113877.37 , 38836.92, 105539.25)          
    #        Ef1   ,   Ef2    ,    Eb2    ,     Ef3  ,      Eb3   ,    Ef4    ,   Eb4     ,    Ef5   ,    Eb5     ,     Ef_5  ,   Eb_5   ,   Ef6  ,     Eb6   ,   Ef_6     ,   Eb_b   ,  Ef7   ,     Ef8  
    #           0       1          2             3           4         5          6            7          8              9          10       11         12          13          14         15        16
        # coverage correction
    function Exp(x)
        y2 =  1 - 1/exp(900*(x-0.02)^6); #y2 =  1 - 1/exp(1100*(x-0.09)^6) #
    end 

elseif  (catalyst=="Co")
    frequency = DataFrame(CSV.File("FRQ/Co_batch_new.csv";  header=1, delim=","))
    #Hydrogen coverage
    alpha_H= -0.04*96.48530749926e3; intercept_H = 0#.6622*2*96.48530749926e3 #convert eV to KJ/mol
    #Nitrogen coverage
    alpha_N= -1.67*96.48530749926e3; intercept_N = 0#1.1521*96.48530749926e3*2 #convert eV to KJ/mol
    #activation barriers in J/mol
    Eng = (59090.61, 96759.64, 105556.77, 53770.92, 97319.09, 102829.23, 96307.11, 172937.14, 114111.42, 173028.61, 33342.57, 50257.98, 88011.05, 52917.84, 114377.27, 35448.65, 87974.69); # with ZPE
    #Eng = (64775.49, 114327.55, 109517.39 , 64652.05, 97319.09107, 117204.77, 100692.3631, 180015.69, 117885.9204 ,  177932.09, 36497.35 ,51114.96, 92935.27148   ,  70985.64  , 120280.40 , 41466.81, 95063.37)     #Co2Ni     
    #        Ef1   ,   Ef2    ,    Eb2    ,     Ef3  ,      Eb3   ,    Ef4    ,   Eb4     ,    Ef5   ,    Eb5     ,     Ef_5  ,   Eb_5   ,     Ef6    ,     Eb6   ,   Ef_6     ,   Eb_b    ,  Ef7     ,     Ef8  
    function Exp(x)
        y2 =  1 - 1/exp(1100*(x-0.09)^6)
    end 

elseif  (catalyst=="Co2Ni")
    # frequency
    frequency = DataFrame(CSV.File("FRQ/Co2Ni_batch.csv";  header=1, delim=","))
    #Hydrogen coverage
    alpha_H= -0.07*96.48530749926e3; 
    intercept_H = 0#.6622*2*96.48530749926e3 #convert eV to KJ/mol    
    #Nitrogen coverage
    alpha_N= -1.95*96.48530749926e3; 
    intercept_N = 0#1.1521*96.48530749926e3*2 #convert eV to KJ/mol    
    #activation barriers in J/mol
    Eng = (62309.61, 92947.21, 112026.5211, 48069.46, 90932.91768,107626.23, 89265.30712, 136522.78,	76208.17, 159849.23, 34550.57235, 16767.52, 69033.37246, 68159.10, 124683.266, 43509.63, 92439.05); #with ZPE
    #Eng = (59305.89839, 92947.21, 112026.5211, 48069.46, 90932.91768, 107626.23, 89265.30712, 136522.78, 76157.64205, 159849.23, 34550.57235, 16767.52, 69033.37246, 68159.10, 124683.266, 43509.63, 92439.05)
    # Eng = (71180.58, 110886.82, 115082.9827, 35554.10, 67174.02725, 122610.57, 96010.14438, 146357.62, 78859.90605, 164903.23 , 34550.57 ,  16949.01  ,74596.32934,  86662.28  , 133130.46 , 48593.48 , 100539.09)     #Co2Ni     
    #        Ef1   ,   Ef2    ,    Eb2    ,     Ef3  ,      Eb3   ,    Ef4    ,   Eb4     ,    Ef5   ,    Eb5     ,     Ef_5  ,   Eb_5   ,     Ef6    ,     Eb6   ,   Ef_6     ,   Eb_b    ,  Ef7     ,     Ef8  
    function Exp(x)
        y2 =  1 - 1/exp(90*(x+0.08)^6)
    end 

elseif  (catalyst=="CoNi2")
    # frequency
    frequency = DataFrame(CSV.File("FRQ/CoNi2_batch.csv";  header=1, delim=","))
    #Hydrogen coverage
    alpha_H= -0.07*96.48530749926e3; intercept_H = 0#.6622*2*96.48530749926e3 #convert eV to J/mol    
    #Nitrogen coverage
    alpha_N= -2.1*96.48530749926e3; intercept_N = 0#1.1521*96.48530749926e3*2 #convert eV to J/mol    
    #activation barriers
    Eng =(62309.61, 96798.68, 105556.4972, 49846.33, 86830.51582, 109104.50, 87020.72599, 148476.60, 120647.9918, 159309.26, 36228.08713, 5016.20, 62072.86448, 41516.34, 130724.8527, 41590.04, 94489.15); # with ZPE
    #Eng =(57838.96279, 96798.68, 105556.4972, 49846.33, 86830.51582, 109104.50, 87020.72599, 148476.60, 120647.9918, 159309.26, 36228.08713, 5016.20, 62072.86448, 41516.34, 130724.8527, 41590.04, 94489.15)
    # Eng = (71230.68 , 112201.60, 108240.9114, 69100.26, 94594.785  , 124616.15 , 93565.42  , 155932.40, 122701.49  , 162956.50 , 39377.17 ,  5249.21  , 66516.49   ,  64427.06  ,  140693.62, 46879.27 , 102358.68)          
    #    Ef1 71230.68   ,   Ef2    ,    Eb2    ,     Ef3  ,      Eb3   ,    Ef4    ,   Eb4     ,    Ef5   ,    Eb5     ,     Ef_5  ,   Eb_5   ,     Ef6    ,     Eb6   ,   Ef_6     ,   Eb_b    ,  Ef7     ,     Ef8  
    function Exp(x)
        y2 =  1 - 1/exp(900*(x-0.02)^6)
    end 
else
    println("ERROR! Your catalyst is not known")
end
println("---Start collecting Frequencies----")
# ============================================================================= 
# # 0/ Constants
# =============================================================================
# step 1 - NH3 + * <=> NH3*        Ammonia adsorption
frq_IS1 =  collect(skipmissing(frequency.IS1)) ; 
frq_FS1 =  collect(skipmissing(frequency.IS2)) ; 

# step 2 - NH3* + * <=> NH2* + H*   Dehydrogenation
frq_IS2 =  collect(skipmissing(frequency.IS2)) ; 
frq_TS2 =  collect(skipmissing(frequency.TS2)) ; 
frq_FS2 =  collect(skipmissing(frequency.FS2)) ; 

#  step 3 -NH2* + * <=> NH* + H*    Dehydrogenation
frq_IS3 =  collect(skipmissing(frequency.IS3)) ;  
frq_TS3 =  collect(skipmissing(frequency.TS3)) ; 
frq_FS3 =  collect(skipmissing(frequency.FS3)) ; 

#  step 4 - NH* + * <=> NH* + H*    Dehydrogenation
frq_IS4 =  collect(skipmissing(frequency.IS4)) ;  
frq_TS4 =  collect(skipmissing(frequency.TS4)) ; 
frq_FS4 =  collect(skipmissing(frequency.FS4)) ; 

#  step 5 - N* + N* <=> *N-N*       N-N recombination
frq_IS5 =  collect(skipmissing(frequency.IS5)) ;  
frq_TS5 =  collect(skipmissing(frequency.TS5)) ; 
frq_FS5 =  collect(skipmissing(frequency.FS5)) ; 


#  step 5_ - NH* + N* <=> *NH-N*    NH-N recombination
frq_IS_5 =  collect(skipmissing(frequency.IS_5)) ;  
frq_TS_5 =  collect(skipmissing(frequency.TS_5)) ; 
frq_FS_5 =  collect(skipmissing(frequency.FS_5)) ; 

#  step 6 - *N-N* <=> N2*  +*       N-N rotation
frq_IS6 =  collect(skipmissing(frequency.FS5)) ; 
frq_TS6 =  collect(skipmissing(frequency.TS6)) ; 
frq_FS6 =  collect(skipmissing(frequency.FS6)) ; 

#  step 6_ - *NH-N* <=> N2*  +H*       N-N rotation
frq_IS_6 =  collect(skipmissing(frequency.FS_5)) ; 
frq_TS_6 =  collect(skipmissing(frequency.TS_6)) ; 
frq_FS_6 =  collect(skipmissing(frequency.FS_6)) ; 

#  step 7 - N2* <=> N2 (gas)  +*       N2 desorbtion
frq_IS7 =  collect(skipmissing(frequency.FS6)) ; 
frq_FS7 =  collect(skipmissing(frequency.FS7)) ; 

#  step 8 - H* + H* <=> H2 (gas) + 2*       N-N recombination
frq_IS8 =  collect(skipmissing(frequency.IS8)) ; 
frq_FS8 =  collect(skipmissing(frequency.FS8)) ; 

#------------------------------------------------------------------------------
Ef1  = Eng[1] # +(ZPE(frq_FS1)-ZPE(frq_IS1))

Ef2 = Eng[2]  # +(ZPE(frq_TS2)-ZPE(frq_IS2))
Eb2 = Eng[3]  # +(ZPE(frq_TS2)-ZPE(frq_FS2))

Ef3 = Eng[4]  # +(ZPE(frq_TS3)-ZPE(frq_IS3))
Eb3 = Eng[5]  # +(ZPE(frq_TS3)-ZPE(frq_FS3))

Ef4 = Eng[6]  # +(ZPE(frq_TS4)-ZPE(frq_IS4))
Eb4 = Eng[7]  # +(ZPE(frq_TS4)-ZPE(frq_FS4)) 
 
Ef5 = Eng[8]  # +(ZPE(frq_TS5)-ZPE(frq_IS5)) 
Eb5 = Eng[9]  # +(ZPE(frq_TS5)-ZPE(frq_FS5))  

Ef_5 = Eng[10]  # +(ZPE(frq_TS_5)-ZPE(frq_IS_5)) 
Eb_5 = Eng[11]  # +(ZPE(frq_TS_5)-ZPE(frq_FS_5)) 
 
Ef6 = Eng[12]  # +(ZPE(frq_TS6)-ZPE(frq_IS6))    
Eb6 = Eng[13]  # +(ZPE(frq_TS6)-ZPE(frq_FS6))

Ef_6 = Eng[14]  # +(ZPE(frq_TS_6)-ZPE(frq_IS_6))
Eb_6 = Eng[15]  # +(ZPE(frq_TS_6)-ZPE(frq_FS_6)) 
    
Ef7 = Eng[16]  # +(ZPE(frq_FS7)-ZPE(frq_IS7)) 

Ef8 = Eng[17]  # +(ZPE(frq_FS8)-ZPE(frq_IS8))


eV = 96.48530749926e3 #J/mol
arguments =(Ef1    ,   Ef2    ,    Eb2   ,     Ef3,    Eb3   ,    Ef4    ,   Eb4   ,    Ef5   ,    Eb5   ,   Ef_5,   Eb_5  ,   Ef6  ,   Eb6  , Ef_6  ,   Eb_6 ,   Ef7  ,  Ef8  )

    # mass of relevant reactants
m_NH3 = 17.03 * 1.66054e-27
m_N2 = 14 * 1.66054e-27
m_H2 = 2 * 1.66054e-27
m_H = 1.66054e-27
println("---Start ODE----\n")
# =============================================================================
#   FINISH OUTPUTFILE
# =============================================================================
th_NH3_0=0  # Coverages, Pre-adsorbed in ML
theta = 1-th_NH3_0
y0 = [th_NH3_0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0; 1-th_NH3_0 ; 0 ; 0; P]; #initial conditions

#path to file where coverage will be saved
batch_file_path =string(catalyst,"/",string("Batch_Cov_",catalyst,"_",string(P),"Pa_",t_endstr,"s_noN0.dat"))
rc_file_path=string(catalyst,"/",string("RC_",catalyst,"_",string(P),"Pa_",t_endstr,"s_noN0.dat"))
PE_file_path=string(catalyst,"/",string("pre-exp_",catalyst,"_",string(P),"Pa_",t_endstr,"s_noN0.dat"))
dG_file_path=string(catalyst,"/",string("dG0_",catalyst,"_",string(P),"Pa_",t_endstr,"s_noN0.dat"))
dG_RT_file_path=string(catalyst,"/",string("dG_RT_",catalyst,"_",string(P),"Pa_",t_endstr,"s_noN0.dat"))

#create file
header= "T Th_NH3 Th_NH2 Th_NH Th_N Th_NN Th_NNH Th_N2 Th_H Th P_N2 P_H2 P_NH3"

df=DataFrame(Temp=missing,Th_NH3= missing  ,Th_NH2= missing  ,Th_NH= missing  ,Th_N= missing  ,Th_NN= missing  ,Th_NNH= missing  ,Th_N2= missing  ,Th_H= missing  ,Th= missing  ,P_N2= missing  ,P_H2= missing  ,P_NH3= missing)

rate_const=DataFrame(Temp=missing,k_ads_1=missing , k_des_1=missing , k_f_2=missing , k_b_2=missing , k_f_3=missing , k_b_3=missing , k_f_4=missing , k_b_4=missing , k_f_5=missing , k_b_5=missing , k_f5=missing , k_b5=missing , k_f_6=missing , k_b_6=missing , k_f6=missing , k_b6=missing , k_f_7=missing , k_b_7=missing , k_f_8=missing , k_b_8=missing)

pre_exp= DataFrame(Temp=missing,preE_ads_1=missing , preE_des_1=missing , preE_f_2=missing , preE_b_2=missing , preE_f_3=missing , preE_b_3=missing , preE_f_4=missing , preE_b_4=missing , preE_f_5=missing , preE_b_5=missing , preE_f5=missing , preE_b5=missing , preE_f_6=missing , preE_b_6=missing , preE_f6=missing , preE_b6=missing , preE_f_7=missing , preE_b_7=missing , preE_f_8=missing , preE_b_8=missing)

dG= DataFrame(Temp=missing, dG_ads_1=missing , dG_des_1=missing , dG_f_2=missing , dG_b_2=missing , dG_f_3=missing , dG_b_3=missing , dG_f_4=missing , dG_b_4=missing , dG_f_5=missing , dG_b_5=missing , dG_f5=missing , dG_b5=missing , dG_f_6=missing , dG_b_6=missing , dG_f6=missing , dG_b6=missing , dG_f_7=missing , dG_b_7=missing , dG_f_8=missing , dG_b_8=missing)


dG_RT= DataFrame(Temp=missing, dG_RT_ads_1=missing , dG_RT_des_1=missing , dG_RT_f_2=missing , dG_RT_b_2=missing , dG_RT_f_3=missing , dG_RT_b_3=missing , dG_RT_f_4=missing , dG_RT_b_4=missing , dG_RT_f_5=missing , dG_RT_b_5=missing , dG_RT_f5=missing , dG_RT_b5=missing , dG_RT_f_6=missing , dG_RT_b_6=missing , dG_RT_f6=missing , dG_RT_b6=missing , dG_RT_f_7=missing , dG_RT_b_7=missing , dG_RT_f_8=missing , dG_RT_b_8=missing)

#rate=DataFrame(Temp=missing,rN2=missing , rH2=missing , rNH3_consumption=missing)
t0 = time()
@time for T in 240:20:1200
    # =============================================================================
    # 4/ rate constant calculation
    # =============================================================================   
    println("---Writing rate constanstant file----")
    #P0= 1e5
    
    # calculate all reaction rate constants
        #Dissociation of NH3
    k_ads_1 = k_ads(T, P0, A, m_NH3)/Beta
    k_des_1 = k_ads(T, P0, A, m_NH3)* EXP(Ef1,T)* (Qvib(T,frq_IS1))/Qvib(T,frq_FS1)* Qtran(T,m_NH3)*QrotNL(T,3,9.444^2*6.196)/Beta  
            
    k_f_2 = k_surf(T,  Ef2) * Qvib(T, frq_TS2) /Qvib(T, frq_IS2)/Beta
    k_b_2 = k_surf(T,  Eb2) * Qvib(T, frq_TS2) /Qvib(T, frq_FS2)/Beta  # dissociation of NH3
    
    k_f_3 = k_surf(T,  Ef3) * Qvib(T, frq_TS3) /Qvib(T, frq_IS3)/Beta
    k_b_3 = k_surf(T,  Eb3) * Qvib(T, frq_TS3) /Qvib(T, frq_FS3)/Beta  # dissociation of NH3
        
    k_f_4 = k_surf(T,  Ef4)  * Qvib(T, frq_TS4) /Qvib(T, frq_IS4)/Beta
    k_b_4 = k_surf(T,  Eb4)  * Qvib(T, frq_TS4) /Qvib(T, frq_FS4)/Beta # dissociation of NH
    
    #------------------------------------
    k_f_5 = k_surf(T, Ef5) * Qvib(T, frq_TS5) /Qvib(T, frq_IS5)/Beta #+Th_NH +Th_H
    k_b_5 = k_surf(T, Eb5) * Qvib(T, frq_TS5) /Qvib(T, frq_FS5)/Beta  # N-N coupling
    
    k_f5 = k_surf(T, Ef_5) * Qvib(T, frq_TS_5) /Qvib(T, frq_IS5)/Beta #!!!
    k_b5 = k_surf(T, Eb_5) * Qvib(T, frq_TS_5) /Qvib(T, frq_FS5)/Beta # NH-N coupling
    
    k_f_6 = k_surf(T, Ef6) * Qvib(T, frq_TS6) /Qvib(T, frq_IS6)/Beta
    k_b_6 = k_surf(T, Eb6) * Qvib(T, frq_TS6) /Qvib(T, frq_FS6)/Beta  # N-N coupling
    
    k_f6 = k_surf(T, Ef_6) * Qvib(T, frq_TS_6) /Qvib(T, frq_IS6)/Beta #!!!
    k_b6 = k_surf(T, Eb_6) * Qvib(T, frq_TS_6) /Qvib(T, frq_FS6)/Beta  #NH-H = NN + H
    
    k_f_7 = k_ads(T, P0, A, m_N2)* EXP(Ef7,T)*(Qvib(T,frq_FS7))/Qvib(T,frq_IS7)* Qtran(T,m_N2)*Qrot(T,2,1.9982)/Qrot(T,1,1.9982) 
    k_b_7 = k_ads(T, P0, A, m_N2)/Beta  #N2 adsorption
    
    #------------------------------------    
    k_f_8 = k_surf(T, Ef8+0.1*Joule ) *Qvib(T,frq_FS8)/Qvib(T,frq_IS8)*Qrot(T,2,60.853)*Qtran(T,m_H2)/Beta  #  H-H coupling and desorption
    k_b_8 = k_ads(T, P0, A, m_H2)/Beta
    # Rotational constant B source: http://www.lifesci.sussex.ac.uk/research/fluorine/p5qsp3l/sw_teaching/f1177_html/rotlab/node15.html

    println("---Writing rate constant file----")
    push!(rate_const,(T, k_ads_1, k_des_1, k_f_2, k_b_2, k_f_3, k_b_3, k_f_4, k_b_4, k_f_5, k_b_5, k_f5, k_b5, k_f_6, k_b_6, k_f6, k_b6, k_f_7, k_b_7, k_f_8, k_b_8),promote=true)

    # =============================================================================
    # 1/ functionine function
    # =============================================================================   
    println("---Run Batch at ",string(T),"K---")
    function dxdt!(dy, y,p, t)         # = @ode_def begin
        dxdt = zeros(12)
    
        #θ_NH3, θ_NH2, θ_NH, θ_N, θ_NN, θ_NNH, θ_N2, θ_H,  θ,  P_N2,  P_H2, P_NH3 = y
        Th_NH3, Th_NH2, Th_NH, Th_N, Th_NN, Th_NNH, Th_N2, Th_H,  Th,  P_N2,  P_H2, P_NH3 = y 
        #k_ads_1, k_des_1, k_f_2, k_b_2, k_f_3, k_b_3, k_f_4, k_b_4, k_f_5, k_b_5, k_f5, k_b5, k_f_6, k_b_6, k_f6, k_b6, k_f_7, k_b_7, k_f_8, k_b_8 = p
        
        rf1 = k_ads_1 * P_NH3/P0 * Th;  # adsorption and desorption of NH3
        rb1 = k_des_1 * Th_NH3 ;
        r1 =  3*(rf1-rb1) # #
        
        rf2 = k_f_2 * Th_NH3 * Th * N0; #good
        rb2 = k_b_2 * Th_NH2 * Th_H * N0; 
        r2 = 3*(rf2-rb2); # K2=rf2/rb2
        
        rf3 = k_f_3 * Th_NH2 * Th* N0; #good
        rb3 = k_b_3 * Th_NH * Th_H* N0 ; 
        r3 = 3*(rf3-rb3); #K3=rf3/rb3
        
        rf4 = k_f_4 * Th_NH * Th* N0;
        rb4 = k_b_4 * Th_N * Th_H* N0 ; 
        r4 = 2*(rf4-rb4); #K4=rf4/rb4
        
        rf5 = k_f_5 * Th_N^2 * N0 ; 
        rb5 = k_b_5 * Th_NN * N0 ; 
        r5 = 0.5*(rf5-rb5); 
        
        rf_5 = k_f5 * Th_NH *Th_N* N0;  # N* + NH* <=> *N-NH* 
        rb_5 = k_b5 * Th_NNH * N0 ; 
        r_5 = 1*(rf_5 - rb_5)
        
        rf6 = k_f_6 * Th_NN* N0;
        rb6 = k_b_6 * Th_N2 * Th * N0 ; 
        r6 = 0.5*(rf6-rb6); # K6=rf6/rb6
        
        rf_6 = k_f6 * Th_NNH^2 * N0 ; #  *N-NH* <=> N2* + H
        rb_6 = k_b6 * Th_N2 * Th_H * N0;  
        r_6 = 1*(rf_6 - rb_6)
        
        rf7 = k_f_7 * Th_N2;
        rb7 = k_b_7 * P_N2/P0 * Th  ; 
        r7 = 1.5*(rf7 - rb7)
        
        rf8 = k_f_8 * Th_H^2* N0 ;
        rb8 = k_b_8 * Th^2 * P_H2/P0 * N0; 
        r8 = 4.5*(rf8 - rb8)
        
        dy[1] =  r1 - r2               # dθ(NH3)/dt
        dy[2] =  r2 - r3               # dθ(NH2)/dt
        dy[3] =  r3 - r4 - r_5         # dθ(NH)/dt
        dy[4] =  r4 - 2*r5 - r_5       # dθ(N)/dt 2*r5
        
        dy[5] =  r5 - r6 # = 0         # dθ(N-N)/dt
        dy[6] =  r_5 - r_6             # dθ(N-NH)/dt
        dy[7] =  r6 + r_6 - r7         # dθ(N2)/dt
        dy[8] =  r2 + r3 + r4 + r_6 - 2*r8   # dθ(H)/dt
        dy[9] =  - (dy[1] + dy[2]+dy[3]+dy[4]+ dy[5]+ dy[6]+dy[7]+dy[8])
        
        dy[10] =  r7  # P_N2
        dy[11] =  r8  # P_H2
        dy[12] =  -r1  # P_NH3
        nothing
    end 
        
    # =============================================================================
    # 2/ Solve the differntial equation
    # =============================================================================
    println("---Start solving ODE----\n\n")
    # for k in range(len(th_NH3_0))
    tspan = (1e-5, t_end) #;(0.0,t_end)
    #p=[k_ads_1, k_des_1, k_f_2, k_b_2, k_f_3, k_b_3, k_f_4, k_b_4, k_f_5, k_b_5, k_f5, k_b5, k_f_6, k_b_6, k_f6, k_b6, k_f_7, k_b_7, k_f_8, k_b_8]
    prob = ODEProblem(dxdt!, (y0), (tspan))#, [k_ads_1, k_des_1, k_f_2, k_b_2, k_f_3, k_b_3, k_f_4, k_b_4, k_f_5, k_b_5, k_f5, k_b5, k_f_6, k_b_6, k_f6, k_b6, k_f_7, k_b_7, k_f_8, k_b_8]);
    #sol = solve(prob, ImplicitHairerWannerExtrapolation(),  reltol=(1e-6), abstol=(1e-16), dtmin = 1e-18)
    #end
    _prob = remake(prob)
    sol=solve(_prob, ImplicitHairerWannerExtrapolation(), reltol=1e-6, abstol=1e-18, 
            saveat=t_end,sensealg=ForwardDiffSensitivity(),  dtmin = 1e-20)
    #dtmax = 1.0, dtmin = 2e-6, dtmax = 1e-1 #ImplicitHairerWannerExtrapolation #QuadratureAdjoint
    #Th_NH3, Th_NH2, Th_NH, Th_N, Th_NN, Th_NNH, Th_N2, Th_H,  Th,  P_N2,  P_H2, P_NH3 = r.y 
    #   1       2       3      4      5      6     7     8    9     10     11
    t = sol.t
    #T = Beta * sol.t .+ T0
    
    # =============================================================================
    # 3/ converage and pressure as function of Temperature
    # =============================================================================   
    println("---Writing Batch file of coverage and pressure----")

    push!(df,(T, last(sol[1,:]), last(sol[2,:]), last(sol[3,:]), last(sol[4,:]), last(sol[5,:]), last(sol[6,:]), last(sol[7,:]), last(sol[8,:]), last(sol[9,:]), last(sol[10,:]), last(sol[11,:]), last(sol[12,:])), promote=true)
    
        
    # =============================================================================
    # 5/ calculate pre-exponential factors
    # =============================================================================   
    println("---Calculating Pre-exponential----")
    
    preE_ads_1= k_ads_1/EXP(0.1*Joule,T)
    preE_des_1= k_des_1/EXP(Ef1+0.1*Joule,T)
    
    preE_f_2= k_f_2/EXP(Ef2,T)
    preE_b_2= k_b_2/EXP(Eb2,T)
    
    preE_f_3= k_f_3/EXP(Ef3,T)
    preE_b_3= k_b_3/EXP(Eb3,T)
    
    preE_f_4= k_f_4/EXP(Ef4,T)
    preE_b_4= k_b_4/EXP(Eb4,T)
    
    preE_f_5= k_f_5/EXP(Ef5,T)
    preE_b_5= k_b_5/EXP(Eb5,T)
    
    preE_f5= k_f5 /EXP(Ef_5,T)
    preE_b5= k_b5 /EXP(Eb_5,T)
    
    preE_f_6= k_f_6/EXP(Ef6,T)
    preE_b_6= k_b_6/EXP(Eb6,T)
    
    preE_f6= k_f6 /EXP(Ef_6,T)
    preE_b6= k_b6 /EXP(Eb_6,T)
    
    preE_f_7= k_f_7/EXP(Ef7+0.1*Joule,T)
    preE_b_7= k_b_7/EXP(0.1*Joule,T)
    
    preE_f_8= k_f_8/EXP(Ef8+0.1*Joule ,T) #+ alpha_H*Th_H
    preE_b_8= k_b_8/EXP(0.1*Joule,T)
        
    
    println("---Append Pre-exponential to data frame----")
    push!(pre_exp,(T, preE_ads_1, preE_des_1, preE_f_2, preE_b_2, preE_f_3, preE_b_3, preE_f_4, preE_b_4, preE_f_5, preE_b_5, preE_f5, preE_b5, preE_f_6, preE_b_6, preE_f6, preE_b6, preE_f_7, preE_b_7, preE_f_8, preE_b_8),promote=true)
    
    
    # =============================================================================
    # 6/ calculate Free Gibbs energy
    # =============================================================================   
    println("---Calculating Gibbs free energy/RT AND Gibbs free energy in J/mol----")
    cst = (kb/h*T)
    
    dG_RT_ads_1 = -log(EXP(0.1*Joule,T)); dG_ads_1 = R*T* dG_RT_ads_1
    dG_RT_des_1 = -log(k_des_1/k_ads(T, P0, A, m_NH3)*EXP(0.1*Joule,T)); dG_des_1 = R*T* dG_RT_des_1
    
    dG_RT_f_2 = -log(k_f_2/cst); dG_f_2 = R*T* dG_RT_f_2
    dG_RT_b_2 = -log(k_b_2/cst); dG_b_2 = R*T* dG_RT_b_2
    
    dG_RT_f_3 = -log(k_f_3/cst); dG_f_3 = R*T* dG_RT_f_3
    dG_RT_b_3 = -log(k_b_3/cst); dG_b_3 = R*T* dG_RT_b_3
    
    dG_RT_f_4 = -log(k_f_4/cst); dG_f_4 = R*T* dG_RT_f_4
    dG_RT_b_4 = -log(k_b_4/cst); dG_b_4 = R*T* dG_RT_b_4
    
    dG_RT_f_5 = -log(k_f_5/cst); dG_f_5 = R*T* dG_RT_f_5
    dG_RT_b_5 = -log(k_b_5/cst); dG_b_5 = R*T* dG_RT_b_5
    
    dG_RT_f5 = -log(k_f5/cst); dG_f5 = R*T* dG_RT_f5
    dG_RT_b5 = -log(k_b5/cst); dG_b5 = R*T* dG_RT_b5
    
    dG_RT_f_6 = -log(k_f_6/cst); dG_f_6 = R*T* dG_RT_f_6
    dG_RT_b_6 = -log(k_b_6/cst); dG_b_6 = R*T* dG_RT_b_6
    
    dG_RT_f6 = -log(k_f6/cst); dG_f6 = R*T* dG_RT_f6
    dG_RT_b6 = -log(k_b6/cst); dG_b6 = R*T* dG_RT_b6
    
    dG_RT_f_7 = -log(k_f_7 / k_ads(T, P0, A, m_N2)*EXP(0.1*Joule,T)); dG_f_7 = R*T* dG_RT_f_7
    dG_RT_b_7 = -log(EXP(0.1*Joule,T)); dG_b_7 = R*T* dG_RT_b_7
    
    dG_RT_f_8 = -log(k_f_8/cst); dG_f_8 = R*T* dG_RT_f_8
    dG_RT_b_8 = -log(EXP(0.1*Joule,T)); dG_b_8 = R*T* dG_RT_b_8

    
    println("---Append DG/RT  to data_frame----")
    push!(dG_RT,(T,dG_RT_ads_1 , dG_RT_des_1 , dG_RT_f_2 , dG_RT_b_2 , dG_RT_f_3 , dG_RT_b_3 , dG_RT_f_4 , dG_RT_b_4 , dG_RT_f_5 , dG_RT_b_5 , dG_RT_f5 , dG_RT_b5 , dG_RT_f_6 , dG_RT_b_6 , dG_RT_f6 , dG_RT_b6 , dG_RT_f_7 , dG_RT_b_7 , dG_RT_f_8 , dG_RT_b_8 ), promote=true)
    
    
    println("---Append Gibbs free enrgy to data_frame----")
    push!(dG,(T, dG_ads_1 , dG_des_1 , dG_f_2 , dG_b_2 , dG_f_3 , dG_b_3 , dG_f_4 , dG_b_4 , dG_f_5 , dG_b_5 , dG_f5 , dG_b5 , dG_f_6 , dG_b_6 , dG_f6 , dG_b6 , dG_f_7 , dG_b_7 , dG_f_8 , dG_b_8), promote=true)
    
    
    
    println("\n~~~~~~~~~_Calculation at ",string(T),"K IS DONE_~~~~~~~~~~\n")
end

print("\n duration time = ", string(time() - t0))
 

#Write batch file
writedlm(batch_file_path, Iterators.flatten(([names(identity.(dropmissing(df)) )], eachrow(identity.(dropmissing(df)) ))), ' ')

#Write rate constant file
writedlm(rc_file_path, Iterators.flatten(([names(identity.(dropmissing(rate_const)) )], eachrow(identity.(dropmissing(rate_const)) ))), ' ')

#Write Pre-exponential file
writedlm(PE_file_path, Iterators.flatten(([names(identity.(dropmissing(pre_exp)) )], eachrow(identity.(dropmissing(pre_exp)) ))), ' ')

#Write DG/RT file
writedlm(dG_RT_file_path, Iterators.flatten(([names(identity.(dropmissing(dG_RT)) )], eachrow(identity.(dropmissing(dG_RT)) ))), ' ')

#Write Gibbs free energy file
writedlm(dG_file_path, Iterators.flatten(([names(identity.(dropmissing(dG)) )], eachrow(identity.(dropmissing(dG)) ))), ' ')

#----------------------------------------
println("\n\nWhatever you do, work at it with all your heart,\nas working for the Lord,\nnot for human masters. Col 2:23\n\n\n")
