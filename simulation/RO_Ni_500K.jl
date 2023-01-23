#!/usr/bin/env julia

using GRUtils
using DifferentialEquations
using DiffEqSensitivity, ForwardDiff
using LaTeXStrings

#using Gtk
using CSV
using DataFrames
using DelimitedFiles

# What is your Temperature?

print("What is your Catalyst?\n");
catalyst = "Ni" #readline()
print("Your catalyst is :", catalyst)
T = 500 #readline();

#k_rate = open_dialog("Pick a file", GtkNullContainer(), String[])
#k =  last(DataFrame(CSV.File(k_rate; header=1, delim=" ")),1);

kf1=float(1.5632E+08) ; kb1= float(1.7238E+13) ; K1 = kf1/kb1 ;
kf2=float(2.4902E+02) ; kb2= float(2.4638E+03) ; K2 = kf2/kb2 ;
kf3=float(4.2750E+06) ; kb3= float(7.5569E+04) ; K3 = kf3/kb3 ;
kf4=float(1.1662E+02) ; kb4= float(1.0602E+05) ; K4 = kf4/kb4 ;
kf5=float(6.6260E-04) ; kb5= float(6.3836E-03) ; K5 = kf5/kb5 ;
kf6=float(3.2208E-06) ; kb6= float(2.9432E+06) ; K6 = kf6/kb6 ;
kf7=float(2.8969E+14) ; kb7= float(2.3095E+10) ; K7 = kf7/kb7 ;
kf8=float(1.4556E+08) ; kb8= float(5.2296E+02) ; K8 = kf8/kb8 ;
kf9=float(1.7299E+14) ; kb9= float(1.7241E+08) ; K9 = kf9/kb9 ;
kf10=float(6.7091E+10); kb10= float(4.5615E+08)  ; K10 = kf10/kb10 ; 

#K = Array([K1,K2, K3, K4, K5, K6, K7, K8, K9, K10])
kf = Array([kf1,kf2, kf3, kf4, kf5, kf6, kf7, kf8, kf9, kf10])
kr = Array([ kb1, kb2,  kb3,  kb4,  kb5,  kb6,  kb7,  kb8,  kb9,  kb10])
#press0 = Array([0.000, 0.0001, 1])
#press0 = Array([0.000, 0.00004, 1])
press0 = Array([0.000, 0.00, 1])

function calc_rate2(y, press)
    
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
        
    coef=Array([3, 3, 3, 2, 0.5, 1, 0.5, 1, 1.5, 4.5])    
    fwd_factor = Array([P_NH3* Th,
                Th_NH3 * Th, 
                Th_NH2 * Th, 
                Th_NH * Th, 
                Th_N^2, 
                Th_NH *Th_N, 
                Th_NN,
                Th_NNH,
                Th_N2,
                Th_H^2])
    rev_factor = Array([Th_NH3, 
                Th_NH2 * Th_H, 
                Th_NH * Th_H, 
                Th_N * Th_H, 
                Th_NN,
                Th_NNH, 
                Th_N2 * Th,
                Th_N2 * Th_H,
                P_N2* Th,
                Th^2 * P_H2])
    r = (coef .* kf .* fwd_factor .- coef .*kr .* rev_factor)
    return r
end

function calc_overall_rate2(y, press)
    P_N2 = press[1]
    P_H2 = press[2]
    P_NH3= press[3]
    r1 = 3*(kf[1] * y[9] * press[3] - kr[1] * y[1])
    return -r1
end

function state_eqn2(y, press, t)
    r = calc_rate2(y, press)
    eq1 = (r[1]- r[2])              # dθ(NH3])/dt
    eq2 = (r[2]- r[3])              # dθ(NH2])/dt
    eq3 = (r[3]- r[4]- r[6])         # dθ(NH])/dt
    eq4 = (r[4]- 2*r[5]- r[6])       # dθ(N])/dt 2*r[5 
    eq5 = (r[5]- r[7])              # dθ(N]-N])/dt
    eq6 = (r[6]- r[8])              # dθ(N]-NH])/dt
    eq7 = (r[7]+ r[8]- r[9])         # dθ(N2])/dt
    eq8 = (r[2]+ r[3]+ r[4]+ r[8]- 2*r[10])   # dθ(H])/dt
    eq9 = - (eq1+ eq2+ eq3+ eq4+ eq5+ eq6+ eq7+ eq8) # dθ/dt 
    Array([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9])
end

time_end = 100000
tspan = (0., time_end)
th_NH3_0 =  0.0
theta = 1-th_NH3_0  # Initial free sites 
y0 = [th_NH3_0; 0.0; 0.00; 0.00; 0.0; 0.0; 0.0; 0.0; theta]
prob2 = ODEProblem(state_eqn2, y0, tspan, press0);

function rate_wrapper2(p)
    p = exp.(p)
    nt = 10000
    _prob = remake(prob2, p=p)
    theta_sol = Array(solve(_prob, Rodas4P(), saveat=LinRange(10, time_end, nt),
            sensealg=ForwardDiffSensitivity(), reltol=1e-13,abstol=1e-13,dtmin = 1e-25))' # [nt, nu]
    p_matrix = reshape(p, 1, size(p, 1)) # [1, np]
    p_repeat = repeat(p_matrix, nt, 1) # [nt, np]
    r = Array{Real, 1}(undef, nt)
    for i in 1:nt
        r[i] = calc_overall_rate2(theta_sol[i, :], p_repeat[i, :])
    end
    return log.(-r)
end
drdp_N2 = ForwardDiff.jacobian(rate_wrapper2, log.(press0))
DrC = round.(drdp_N2[end, :]; digits=3)

print(string(T, DrC))
