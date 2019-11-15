
#Pkg.add("DiffEqSensitivity ")
include("smc_model.jl")
include("smc_model_parameters.jl")
using DifferentialEquations
using DiffEqSensitivity 
using Plots
gr()
using DataFrames

q0 = Control_params();
p0 = Fixed_params();
s0 = State();
q0_dict = type2dict(q0);
s0_dict = type2dict(s0);
(q0_elevated,s0_elevated) = elevated(q0_dict);


N=1;
t = collect(range(0, stop=10, length=200))
q_range=[[30,40],[1,200],[0,1],[0,0.1],[0,0.01],[0,0.1],[0.01,1],[200,1000],[1000,7000],[1,100],[0,1],[1000,5000],[1,100],[0,10],[0,10],[0,1],[0,0.0001],[0,1],[0,0.1]]
sobol_sensitivity_0(s0,q0,t,q_range,N,1)

prob.p



(output, problem, solution) = solveODEs(
    s0_elevated,q0_elevated,
    timespan=(0.0,60*30),
    fix = ["Ca_mit","Ca_ecs"]
    ,dtmax=100);

prob1 = remake(problem;u0=s0_elevated,p=q0_elevated);
prob1.u0

p1 = plot(output[:t]/60,output[:Ca_cyt],label="Ca_cyt")
p3 = plot(output[:t]/60,output[:other_Ca_cyt_tot],label="Ca_cyt_tot")

p2 = plot(output[:t]/60,output[:Ca_er],label="Ca_er")
p4 = plot(output[:t]/60,output[:other_Ca_er_tot],label="Ca_er_tot")

plot(p1,p2,p3,p4)



p1 = plot(output[:t]/60,output[:J_serca],label="J_serca")
p2 = plot(output[:t]/60,output[:J_ryr],label="J_ryr")
p3 = plot(output[:t]/60,output[:J_ipr],label="J_ipr")
p4 = plot(output[:t]/60,output[:J_leak_er],label="J_er_leak")
plot(p1,p2,p3,p4,layout=(2,2))

q_range=[[30,40],[1,200],[0,1],[0,0.1],[0,0.01],[0,0.1],[0.01,1],[200,1000],[1000,7000],[1,100],[0,1],[1000,5000],[1,100],[0,10],[0,10],[0,1],[0,0.0001],[0,1],[0,0.1]]

q = [(q_range[j][2] -q_range[j][1])*rand() + q_range[j][1] for j in 1:length(q_range)]

 q0_dict[Symbol("SERCA")] = q[1]; ## ER 
        q0_dict[Symbol("SERCA_s")] = q[2];   ## ER 
        q0_dict[Symbol("Qvocc")] = q[3]; #
    # Leaks
       q0_dict[Symbol("g_leak_mit")] = q[4];
       q0_dict[Symbol("g_leak_ecs")] = q[5];
       q0_dict[Symbol("k_leak_er")] = q[6];
    # steady state values
       q0_dict[Symbol("Ca_cyt_infty")] = q[7];# muM
       q0_dict[Symbol("Ca_er_infty")] = q[8];# muM
       q0_dict[Symbol("ATP_infty")] = q[9];
       q0_dict[Symbol("ADP_infty")] = q[10];
       q0_dict[Symbol("Qryr")] = q[11];  # Junming
       q0_dict[Symbol("Qip3r")] = q[12]; # /s
       q0_dict[Symbol("MyoTot")] = q[13]; # muM check this! @TODO
       q0_dict[Symbol("Qncx")] = q[14]; # computed using Johny
       q0_dict[Symbol("Qpmca")] = q[15]; # Johny
       q0_dict[Symbol("Vnclx")] = q[16];#0.026673080107946  # Wacquier 2016
       q0_dict[Symbol("Vmcu")] = q[17]; # Wacquier 2016
       q0_dict[Symbol("kncx2")] = q[18];# Demir muA/muM
       q0_dict[Symbol("L")] = 0.0;
        (q0_elevated,s0_elevated) = elevated(q0_dict);

 (q0_elevated,s0_elevated) = elevated(q0_dict);

typeof(s0_elevated)

    test = ODEs(s0,0,q0,fix=["Ca_mit","Ca_ecs"])
    sdict = type2dict(s0)
    ode_vars = keys(test);
    u0 = broadcast(v->sdict[Symbol(v)], ode_vars);

    function ode_problem_fun(du,u,p,t)
        ydict = Dict(Symbol(k)=>v for (k,v) in zip(ode_vars,u))
        ydot = ODEs(reconstruct(s,ydict),t,q,fix=["Ca_mit","Ca_ecs"])
        for (j,v) in enumerate(ode_vars)
            du[j] = ydot[v]
        end
        return du
    end
    timespan = (t[1],t[end]);
    prob = ODEProblem(ode_problem_fun,u0,timespan,q0)

b=type2dict(prob.p)
b.count

prob1 = remake(prob;u0=type2dict(s0_elevated),p=type2dict(q0_elevated);

prob1.p

@assert length(prob.p) == length(p_range)

length(prob.p)

prob.p

prob.p

length.(collect(keys(b)))

map(length, keys(b))

type2dict(prob.p).count

Array(solve(prob1,alg_hints=[:stiff];saveat=t))

typeof(q_evaluted)

typeof(s0_evalued)


