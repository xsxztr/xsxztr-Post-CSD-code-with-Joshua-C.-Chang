using DifferentialEquations
using Parameters
using DataStructures
using NLsolve
using DataFrames
using WAV
using Dierckx 
using FFTW
using DSP
using Statistics

include("smc_model_parameters.jl")

s_default = State();
p_default  = Fixed_params();
q_default = Control_params();

"""
    Compute RHS of the ODEs. Returns Dicts as in the python code

"""
function ODEs(state::State = s_default, t = 0.0, q::Control_params = q_default; onlyODE = true, fix = String[], drop = String[],rd_cer=1.0,rd_cmi=1.0,rd_ces=1.0)
    output = OrderedDict{String,DefaultOrderedDict}();
    dot = DefaultOrderedDict{String,Float64}(0.0);
    J = DefaultOrderedDict{String,Float64}(0.0);
    JK = DefaultOrderedDict{String,Float64}(0.0);
    JNa = DefaultOrderedDict{String,Float64}(0.0);
    other = DefaultOrderedDict{String,Float64}(0.0);
    I = DefaultOrderedDict{String,Float64}(0.0);

    p = p_default
    s = state

    rho_c = (1+p.S_CM*p.K_d/(p.K_d+s.Ca_cyt)^2 + 
        p.B_F*p.K_dB/(p.K_dB+s.Ca_cyt)^2)

    rho_er = (1+p.K_CalrC*p.CalrC / 
        (p.K_CalrC+s.Ca_er)^2+p.K_CalrP * 
        p.CalrP/(p.K_CalrP+s.Ca_er)^2)

    other["rho_c"] = rho_c
    other["rho_er"] = rho_er

    other["Ca_cyt_tot"] = (s.Ca_cyt * 
        (1+p.S_CM/(p.K_d + s.Ca_cyt) +
            p.B_F/(p.K_dB+s.Ca_cyt)))
    other["Ca_er_tot"] = (s.Ca_er * 
        (1+p.CalrC/(p.K_CalrC+s.Ca_er) +
            p.K_CalrP/(p.K_CalrP+s.Ca_er)))

    alpha = 2*p.F*p.Vol_cyt*1e6 # 

    """
    ER Higgins2006
    """

    Pt = q.SERCA
    K1 = sqrt(0.7)
    K3 = sqrt(1.111111e-5)
    k2 = 0.6 * q.SERCA_s
    km2 = 0.97 * q.SERCA_s
    k4 = 0.4*q.SERCA_s
    km4 = 1.2e-3*q.SERCA_s

    J["serca"] = (-Pt*2*
        (-K1^2*K3^2*km2*km4*s.Ca_er^2 + k2*k4*s.Ca_cyt^2)/(
            s.Ca_er^2*s.Ca_cyt^2*K3^2*(k2+km2)+ 
            s.Ca_cyt^2*(k4+k2)+ 
                s.Ca_er^2*K1^2*K3^2*(km2+km4)+K1^2*(k4+km4)
        ))
    # IP3 dynamics
    L = q.L
    #rhor = L*s.Rs/(L+p.Rs_k1)/p.G_zeta/p.G_RT
    rhor = L/(L+p.Rs_k1)
    rh = p.IP3_eta*1/(1+p.IP3_kc/s.Ca_cyt)*s.G

    nu_PLC_delta = (p.V_delta/(1+s.IP3/p.k_delta) * 
        s.Ca_cyt^2/(s.Ca_cyt^2+p.K_PLC^2))
    nu_5P = p.r_5p*s.IP3
    nu_3K = (p.v_3k*(s.Ca_cyt^4/(s.Ca_cyt^4 +
                p.IP3_3K_K_D^4))*s.IP3/(s.IP3+p.IP3_3K_K_3))
    #dot["Rs"] = p.Rs_kr*p.G_RT-(p.Rs_kr+p.Rs_kp*L/(p.Rs_k1+L))*s.Rs-p.Rs_kr*s.Rsp
    #dot["Rsp"] = L*(p.Rs_kp*s.Rs/(p.Rs_k1+L)-p.Rs_ke*s.Rsp/(p.Rs_k2+L))
    dot["G"] = p.G_ka*(p.G_delta+rhor)*(p.G_T-s.G)-p.G_kd*s.G
    dot["IP3"] = rh*s.PIP2 - p.IP3_kdeg*s.IP3
    #dot["IP3"] = rh*s.PIP2 - nu_5P - nu_3K
    dot["PIP2"] = -(rh+p.PIP2_rr )*s.PIP2-p.PIP2_rr*s.IP3+p.PIP2_rr*p.PIP2T

    # IPR
    
    k1 = p.b1/p.a1/s.IP3
    k3 = p.b3/p.a3/s.IP3
    X11 = 1-s.IP3RX10-s.IP3RX00-s.IP3RX01
    dot["IP3RX00"] = ((p.b4*k3+p.b2)*s.IP3RX01/(1+k3)+p.b5*s.IP3RX10 - 
        (p.a4*k1+p.a5+p.a2)*s.Ca_cyt*s.IP3RX00/(1+k1))
    dot["IP3RX10"] = (p.b2*X11+p.a5 * 
        s.Ca_cyt*s.IP3RX00 / (1+k1) 
        -(p.a2*s.Ca_cyt+p.b5)*s.IP3RX10)
    dot["IP3RX01"] = ((p.a4*k1+p.a2)*s.Ca_cyt*s.IP3RX00 /(1+k1) 
        +p.b5*X11 -(p.b4*k3+p.b2+p.a5*s.Ca_cyt)*s.IP3RX01/(1+k3))
    Pip3r = (s.IP3RX10^4+4*s.IP3RX10^3*(1-s.IP3RX10))
    J["ipr"] = q.Qip3r*Pip3r*(s.Ca_er-s.Ca_cyt)
    I["ipr"] = J["ipr"] * alpha

    # RyR
    R00 = 1-(s.RyRR10+s.RyRR11+s.RyRR01)
    PRyR = s.RyRR10^2
    dot["RyRR10"] = (p.Kr1*s.Ca_cyt^2*R00 + p.Kmr2 * 
        s.RyRR11 - s.RyRR10*(p.Kmr1 + p.Kr2*s.Ca_cyt))
    dot["RyRR11"] = (p.Kr2*s.Ca_cyt*s.RyRR10 + p.Kr1 * 
        s.Ca_cyt^2*s.RyRR01 - (p.Kmr2+p.Kmr1)*s.RyRR11)
    dot["RyRR01"] = (p.Kmr1*s.RyRR11 + p.Kr2*s.Ca_cyt * 
        R00 - (p.Kmr2 + p.Kmr1*s.Ca_cyt^2)*s.RyRR01)

    J["ryr"] = q.Qryr*PRyR*(s.Ca_er-s.Ca_cyt)
    I["ryr"] = alpha*J["ryr"]
    #leak
    J["leak_er"] = q.k_leak_er*(s.Ca_er-s.Ca_cyt)
    """
    Mito: MCU/NCLX from Wacquier (Keep these values)
    """
    ECa_mito = (log(s.Ca_mit) - log(max(s.Ca_cyt, 1e-15)))/p.FRT/2

    J["leak_mit"] = -q.g_leak_mit * (s.Phi_mit-ECa_mito)
    J["mcu"] = -(q.Vmcu*(s.Ca_cyt/p.K1)*(1+s.Ca_cyt/p.K1)^3*exp(
        p.p1*s.Phi_mit)/((1+s.Ca_cyt/p.K1)^4+(p.L/(1+abs(s.Ca_cyt)/p.K2)^2.3)))
    I["mcu"] = J["mcu"]*alpha

    J["nclx_mit"] = (q.Vnclx * 
        (s.Ca_mit/s.Ca_cyt)*exp(p.p2*s.Phi_mit))

    """
    ECS
    """
    
    phiF = exp(p.etancx * s.Phi_ecs*p.FRT)
    phiR = exp((p.etancx - 1)*s.Phi_ecs*p.FRT)

    #= Williams
    Incx = (q.Incxbar*((s.Na_cyt/s.Na_ecs)^3*phiF - s.Ca_cyt*phiR/s.Ca_ecs) 
        / (1+(p.Kncxna/s.Na_ecs)^3) 
        / (1+p.Kncxca/s.Ca_ecs) 
        / (1+p.kncxsat*phiR))
    J["ncx_ecs"] = -Incx*p.A_m/(p.F*p.Vol_cyt)
    I["ncx_ecs"] = J["ncx_ecs"] * alpha
    JNa["ncx_ecs"] = -J["ncx_ecs"]

    =#

    # Demir

    Incx = (q.kncx2*(s.Na_cyt^3*s.Ca_ecs*phiF - s.Na_ecs^3*s.Ca_cyt*phiR)
        /(1+p.dncx*(s.Ca_cyt*s.Na_ecs^3 + s.Ca_ecs*s.Na_cyt^3))
        *(1.0/(1+(p.kncx1/s.Ca_cyt)^2))
        )
    I["ncx_ecs"] = Incx
    J["ncx_ecs"] = Incx/(2*p.F*p.Vol_cyt)
    JNa["ncx_ecs"] = -J["ncx_ecs"]

    # Weber NCX

    ECa_ecs = (log(s.Ca_ecs) - log(max(s.Ca_cyt, 1e-12)))/p.FRT/2
    J["leak_ecs"] = -q.g_leak_ecs * (s.Phi_ecs-ECa_ecs)
    J["pmca"] = -q.Qpmca*s.Ca_cyt/(s.Ca_cyt+p.kpmca)
    

    dbl = 1/(1+exp(-s.Phi_ecs/8.3))
    fbl = 1/(1+exp((s.Phi_ecs+42)/9.1))
    J["vocc"] = (-q.Qvocc*dbl*fbl * 
        (s.Phi_ecs - ECa_ecs))
    I["vocc"] = J["vocc"]*alpha

    """
    Myosin-actin
    """
    M = 1 - s.MyoAMp - s.MyoMp - s.MyoAM
    dot["MyoMp"] = (p.MyoK4*s.MyoAMp + p.MyoK1 * 
        (s.Ca_cyt)^3 * 
        M - (p.MyoK2+p.MyoK3)*s.MyoMp)

    dot["MyoAMp"] = (p.MyoK3*s.MyoMp + p.MyoK6 * 
        (s.Ca_cyt)^3 * 
        s.MyoAM - (p.MyoK4+p.MyoK5)*s.MyoAMp)

    dot["MyoAM"] = (p.MyoK5*s.MyoAMp - 
        (p.MyoK7+p.MyoK6*(s.Ca_cyt)^3)*s.MyoAM)
    
    other["Fr"] = s.MyoAMp + s.MyoAM

    """
    Assemble the total fluxes and ODEs
    """
    J["er_cyt"] = rd_cer*(J["serca"] + J["ipr"] + J["ryr"] + J["leak_er"])
    J["mit_cyt"] = rd_cmi*(J["leak_mit"] + J["nclx_mit"] + J["mcu"])
    J["ecs_cyt"] = rd_ces*(J["leak_ecs"] + J["pmca"] + J["ncx_ecs"] + J["vocc"])

    # myosin reaction
    J["cyt_cyt"] = -3*(dot["MyoMp"] + dot["MyoAMp"])*q.MyoTot

    dot["Ca_cyt"] = (J["er_cyt"] + J["mit_cyt"]+
                     J["ecs_cyt"] + J["cyt_cyt"])/rho_c
    dot["Ca_er"] = -(J["er_cyt"])*p.RVer/rho_er
    dot["Ca_ecs"] = -(J["ecs_cyt"])
    dot["Ca_mit"] = -(J["mit_cyt"])*p.RVmit

    dot["Ca_ecs_source"] = dot["Ca_ecs"]
    dot["Ca_mit_source"] = dot["Ca_mit"]

    output["J"] = J
    output["dot"] = dot
    output["JK"] = JK
    output["JNa"] = JNa
    output["other"] = other

    if length(drop)>0
        for val in drop
            delete!(dot,String(val))
        end
    end

    if length(fix)>0
        for val in fix
            dot[val] = 0.0
        end
    end

    if !onlyODE
        return output
    end
    return dot
end

function create_ode_problem_fun(s,q)

end

function solveODEs(s,q; timespan = (0.0,3600),saveat=10,dtmax=1,fix=["Ca_mit","Ca_ecs"],rd_cer=1.0,rd_cmi=1.0,rd_ces=1.0)
    test = ODEs(s,0,q,fix=fix)
    sdict = type2dict(s)
    ode_vars = keys(test);
    u0 = broadcast(v->sdict[Symbol(v)], ode_vars)

    function ode_problem_fun(du,u,p,t)
        ydict = Dict(Symbol(k)=>v for (k,v) in zip(ode_vars,u))
        ydot = ODEs(reconstruct(s,ydict),t,q,fix=fix,rd_cer=rd_cer,rd_cmi=rd_cmi,rd_ces=rd_ces)
        for (j,v) in enumerate(ode_vars)
            du[j] = ydot[v]
        end
        return du
    end

    prob = ODEProblem(ode_problem_fun,u0,timespan)
    solution = solve(prob,alg_hints=[:stiff],dtmax=dtmax,saveat=saveat);

    u = solution.u[1]
    t = solution.t[1]
    ydict = OrderedDict(Symbol(k)=>v for (k,v) in zip(ode_vars,u))
    ydot = ODEs(reconstruct(s,ydict),
                    t,q,fix=["Ca_mit","Ca_ecs"],
                    onlyODE=false,rd_cer=rd_cer,rd_cmi=rd_cmi,rd_ces=rd_ces)
    vals = merge(ydict,Dict( Symbol("J_"*String(k)) => v for (k,v) = ydot["J"] ),
        Dict( Symbol("J_"*String(k)) => v for (k,v) = ydot["J"]),
        Dict( Symbol("other_"*String(k)) => v for (k,v) = ydot["other"]),
    )
    df = DataFrame(vals)
    
    for (t,u) in zip(solution.t[2:end],solution.u[2:end])
        ydict = Dict(Symbol(k)=>v for (k,v) in zip(ode_vars,u))
        ydot = ODEs(reconstruct(s,ydict),
                t,q,fix=["Ca_mit","Ca_ecs"],
                onlyODE=false,rd_cer=rd_cer,rd_cmi=rd_cmi,rd_ces=rd_ces)
        vals = merge(ydict,Dict( Symbol("J_"*String(k)) => v for (k,v) = ydot["J"] ),
            Dict( Symbol("J_"*String(k)) => v for (k,v) = ydot["J"]),
            Dict( Symbol("other_"*String(k)) => v for (k,v) = ydot["other"]),
        )
        push!(df,vals)
    end
    df.t = solution.t
    
    return (df,prob,solution)

end

function balance(sdict, qdict; guess = DefaultOrderedDict{String,Float64}(0.0), excludedODEs = String[], includedJs = String[])
    # 
    test = ODEs(onlyODE=true, drop=excludedODEs)
    s0 = type2dict(s_default)
    q0 = type2dict(q_default)
    dot_vars = keys(test)

    #state_vars = broadcast(v->string(v),keys(s0))
    #param_vars = broadcast(v->string(v),keys(q0))
    state_vars = keys(s0)
    param_vars = keys(q0)
    free_state_vars = setdiff(state_vars,keys(sdict))
    free_params = setdiff(param_vars,keys(qdict))

    dim = length(free_params) + length(free_state_vars)

    rank = length(dot_vars) + length(includedJs)

    if rank != dim
        print(string(rank) * " equations doesn't match " *
        string(dim) * " degrees of freedom")
        return nothing, nothing
    end;

    base_state = copy(s0)  # copy the inputted information
    base_parms = copy(q0)

    for v in keys(sdict)
        base_state[v] = sdict[v]
    end
    for v in keys(qdict)
        base_parms[v] = qdict[v]
    end
    
    function equations!(F,p)
        state = copy(base_state)
        parms = copy(base_parms)

        for (k, v) in zip(free_state_vars, p[1:length(free_state_vars)])
            state[k] = v
        end

        for (k, v) in zip(free_params, p[(length(free_state_vars)+1):end])
            parms[k] = v
        end

        out = ODEs(reconstruct(State(),state), 0, reconstruct(Control_params(),parms), onlyODE = false , drop=excludedODEs)
        for (j,v) in enumerate(dot_vars)
            F[j] = out["dot"][string(v)]
        end
        for (j,v) in enumerate(includedJs)
            F[j + length(dot_vars)] = out["J"][string(v)]
        end
    end

    if length(guess)>0
        initial_vals = broadcast(v->  v in keys(guess) ? guess[v] : s0[v], free_state_vars)
        initial_vals = vcat(initial_vals, broadcast( v-> v in keys(guess) ? guess[v] : 0[v],
                        free_params))
    else
        initial_vals = broadcast(v -> s0[v], free_state_vars)
        initial_vals = vcat(initial_vals, broadcast(v -> Float64(q0[v]), free_params))
    end

    out = nlsolve(equations!, initial_vals)

    return Dict(k => v for (k, v) = zip(free_state_vars, out.zero[1:length(free_state_vars)])), Dict(k => v for (k, v) = zip(free_params, out.zero[(length(free_state_vars)+1):end]))
    
    
end

function initial_condition(y, t, integrator)
    t
end



module balancetest
    function run()
        q0 = Control_params();
        p0 = Fixed_params();
        s0 = State();
        q0_dict = type2dict(q0);
        s0_dict = type2dict(s0);
        delete!(q0_dict,Symbol("g_leak_mit"))
        q0_dict[Symbol("k_leak_er")] = 0
        q0_dict[Symbol("g_leak_ecs")] = 0
        delete!(q0_dict,Symbol("SERCA"))
        delete!(q0_dict,Symbol("Qpmca"))
        q0_dict[Symbol("Qryr")]  = 2
        q0_dict[Symbol("Qip3r")]  = 4
        q0_dict[Symbol("Vnclx")]*=4
        q0_dict[Symbol("MyoTot")] = 10

        delete!(s0_dict,Symbol("IP3"))
        delete!(s0_dict,Symbol("IP3RX00"))
        delete!(s0_dict,Symbol("IP3RX10"))
        delete!(s0_dict,Symbol("IP3RX01"))
        delete!(s0_dict,Symbol("RyRR10"))
        delete!(s0_dict,Symbol("RyRR11"))
        delete!(s0_dict,Symbol("RyRR01"))
        delete!(s0_dict,Symbol("MyoMp"))
        delete!(s0_dict,Symbol("MyoAM"))
        delete!(s0_dict,Symbol("MyoAMp"))
        s0_dict[Symbol("Ca_er")] = 500
        excludedODEs = map(v -> Symbol(v), ["Ca_ecs","Ca_er","Ca_cyt","Phi_ecs","Ca_mit_source","Ca_ecs_source"])
        includedJs = map(v -> Symbol(v), ["er_cyt","ecs_cyt"])
        balance(s0_dict,q0_dict,excludedODEs=excludedODEs,includedJs = includedJs)
    end

end



function elevated(q0_dict; Ca_mit = 0.25)
   # q0_dict = type2dict(q);
    s_default = State();
    s0_dict = type2dict(s_default);
        
    delete!(q0_dict,Symbol("g_leak_mit"));
    delete!(q0_dict,Symbol("g_leak_ecs"));
    delete!(q0_dict,Symbol("SERCA"));
   # q0_dict[Symbol("Qryr")]  = .2;
   # q0_dict[Symbol("SERCA_s")]  = 100;
   # q0_dict[Symbol("Qip3r")]  = 2000;
   # q0_dict[Symbol("MyoTot")] = 10.0;
   # q0_dict[Symbol("Vnclx")] *= 5;
   # q0_dict[Symbol("Qpmca")] *= 100;
   # q0_dict[Symbol("kncx2")] *= .10
   # q0_dict[Symbol("Qvocc")] *= 10;
    #q0_dict[Symbol("L")] = 0.0;
    delete!(s0_dict,Symbol("IP3"));
    delete!(s0_dict,Symbol("IP3RX00"));
    delete!(s0_dict,Symbol("IP3RX10"));
    delete!(s0_dict,Symbol("IP3RX01"));
    delete!(s0_dict,Symbol("RyRR10"));
    delete!(s0_dict,Symbol("RyRR11"));
    delete!(s0_dict,Symbol("RyRR01"));
    delete!(s0_dict,Symbol("MyoMp"));
    delete!(s0_dict,Symbol("MyoAM"));
    delete!(s0_dict,Symbol("MyoAMp"));
    delete!(s0_dict,Symbol("G"));
    delete!(s0_dict,Symbol("PIP2"));
    delete!(s0_dict,Symbol("Rs"));
    delete!(s0_dict,Symbol("Rsp"));
        
    excludedODEs = map(v -> Symbol(v), ["Ca_ecs","Ca_er","Ca_cyt","Phi_ecs","Ca_mit_source","Ca_ecs_source"]);
    includedJs = map(v -> Symbol(v), ["er_cyt","ecs_cyt"]);
    s0_balanced_dict, q0_balanced_dict = balance(s0_dict,q0_dict,excludedODEs=excludedODEs,includedJs = includedJs);
    s0_balanced = reconstruct(s0,merge(s0_dict,s0_balanced_dict));
    q0_balanced = reconstruct(q0,merge(q0_dict,q0_balanced_dict));
       # ODEs(s0_balanced,0,q0_balanced,fix=["Ca_mit","Ca_ecs"],onlyODE=false);
    s0_elevated = reconstruct(s0_balanced, Ca_mit = Ca_mit);
    q0_elevated = reconstruct(q0_balanced, L=0.0);
    s0_elevated_dict = type2dict(s0_elevated);
    return (q0_elevated,s0_elevated);
end

function give_rand_qs(q_range,q0_dict,s0_dict,q_fixed=nothing,indices=nothing)
    
   q = [(q_range[j][2] -q_range[j][1])*rand() + q_range[j][1] for j in 1:length(q_range)];
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
   
   (q0_elevated,s0_elevated) = elevated(q0_dict,s0_dict);
    return (q0_elevated,s0_elevated)
end


function sobol_sensitivity_0(s,q,t,p_range,N,ord=2)
    test = ODEs(s,0,q,fix=["Ca_mit","Ca_ecs"])
    sdict = type2dict(s)
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
    prob = ODEProblem(ode_problem_fun,u0,timespan,q)

    f = function (q)
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
        
       
        test1 = ODEs(s0_elevated,0,q0_elevated,fix=["Ca_mit","Ca_ecs"])
        sdict1 = type2dict(s0_elevated)
        ode_vars1 = keys(test1);
        u1 = broadcast(v->sdict[Symbol(v)], ode_vars1)
        prob1 = remake(prob;u0=u1,p=q0_elevated);
        Array(solve(prob1,alg_hints=[:stiff];saveat=t))
    end
    @assert type2dict(prob.p).count == length(p_range)
    return sobol_sensitivity(f,p_range,N,ord)
end
 


function fluxratesample(N,s0_elevated,q0_elevated,rd_range; timespan = (0.0,60*20),dtmax=1,fix=["Ca_mit","Ca_ecs"])
   # A = Vector{Float64}(undef, N);
    rd_val = Matrix{Float64}(undef,N,3);
    p_cyt_data = Matrix{Float64}(undef,N,4);
    p_mitcyt_data = Matrix{Float64}(undef,N,4);
    p_esccyt_data = Matrix{Float64}(undef,N,4);
    
    for i in 1:N
        #rd = [(rd_range[j][2] -rd_range[j][1])*rand() + rd_range[j][1] for j in 1:length(rd_range)];
        rd = [10^((log10(rd_range[j][2]) -log10(rd_range[j][1]))*rand() + log10(rd_range[j][1])) for j in 1:length(rd_range)]
        (output, problem, solution) = solveODEs(s0_elevated,q0_elevated,timespan=timespan,fix = fix
         ,dtmax=dtmax,rd_cer=rd[1],rd_cmi=rd[2],rd_ces=rd[3]); 
    #    A[i]=maximum(output[:Ca_cyt]);
        rd_val[i,1]=rd[1];
        rd_val[i,2]=rd[2];
        rd_val[i,3]=rd[3];
    
        fs = 1;
        tti=0:1/fs: output[:t][end];
        Ni =length(tti);
     # cyc 
        spl = Spline1D(output[:t], output[:Ca_cyt]; k=2); 
        signal_o= spl(tti);
        
        (freq_dom,p_max,p_min,p_mean)=estimatefr(signal_o,fs);  
        
        p_cyt_data[i,1] = freq_dom;
        p_cyt_data[i,2] = p_max;
        p_cyt_data[i,3] = p_min;
        p_cyt_data[i,4] = p_mean;
        
        
      # j_mitcyt 
        spl_mitcyt = Spline1D(output[:t], output[:J_mit_cyt]; k=2);
        signal_mitcyt_o= spl_mitcyt(tti);
        (freq_mitcyt_dom,p_mitcyt_max,p_mitcyt_min,p_mitcyt_mean)=estimatefr(signal_mitcyt_o,fs); 
        p_mitcyt_data[i,1] = freq_mitcyt_dom;
        p_mitcyt_data[i,2] = p_mitcyt_max;
        p_mitcyt_data[i,3] = p_mitcyt_min;
        p_mitcyt_data[i,4] = p_mitcyt_mean;
        
       #J_esccyt 
        spl_esccyt = Spline1D(output[:t], output[:J_ecs_cyt]; k=2);
        signal_esccyt_o= spl_esccyt(tti);
        (freq_esccyt_dom,p_esccyt_max,p_esccyt_min,p_esccyt_mean)=estimatefr(signal_esccyt_o,fs); 
        p_esccyt_data[i,1] = freq_esccyt_dom;
        p_esccyt_data[i,2] = p_esccyt_max;
        p_esccyt_data[i,3] = p_esccyt_min;
        p_esccyt_data[i,4] = p_esccyt_mean;
        
    end
return  (p_cyt_data,p_mitcyt_data,p_esccyt_data,rd_val)
    
end

function estimatefr(signal_o,fs;start=500)
    Ns = length(signal_o);
    if std(signal_o)!= 0.
        signal_n = (signal_o -mean(signal_o)*ones(Float16, (Ns,1)))/std(signal_o);
    else
        signal_n = (signal_o -mean(signal_o)*ones(Float16, (Ns,1)));
    end
    signal = signal_n[start:end];
    N = length(signal);
    r = periodogram(signal;nfft=N,fs=fs);
    p = r.power;
    p_sort =sortperm(abs.(p));
    freqArray = r.freq;
    freq_dom = freqArray[p_sort[end]];
    
    if freq_dom != 0
        period = 1/freq_dom;
        per_l = ceil(Int,period/fs);
    else
        per_l =  ceil(Int,N/4);
    end
    sig_p = signal_o[max(500,(end-2*per_l)):end];
 #   t_p = tti[end-2*per_l:end];
    p_max = maximum(sig_p);
    p_min = minimum(sig_p);
    #p_A = p_max-p_min;
    p_mean = mean(sig_p);

 return (freq_dom,p_max,p_min,p_mean)
    
end




####################### sensitivity part

function getf(q,prob,t)
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
        
       
        test1 = ODEs(s0_elevated,0,q0_elevated,fix=["Ca_mit","Ca_ecs"])
        sdict1 = type2dict(s0_elevated)
        ode_vars1 = keys(test1);
        u1 = broadcast(v->sdict1[Symbol(v)], ode_vars1)
        prob1 = remake(prob;u0=u1,p=q0_elevated);
        Array(solve(prob1,alg_hints=[:stiff];saveat=t))
end



function getf_new(rd,q,s,t)       
        s_r = reconstruct(s, Ca_mit = rd[4]*0.25)
        if length(rd)>4
            q_r= reconstruct(q, L=rd[5]*0.4)
        else
            q_r =q
        end
        test1 = ODEs(s_r,0,q_r,fix=["Ca_mit","Ca_ecs"])
        sdict1 = type2dict(s_r)
        ode_vars1 = keys(test1);
        u1 = broadcast(v->sdict1[Symbol(v)], ode_vars1)
        function ode_problem_fun_r(du,u,p,t)
            ydict = Dict(Symbol(k)=>v for (k,v) in zip(ode_vars1,u))
            ydot = ODEs(reconstruct(s_r,ydict),t,q_r,fix=fix,rd_cer=rd[1],rd_cmi=rd[2],rd_ces=rd[3])
            for (j,v) in enumerate(ode_vars1)
                du[j] = ydot[v]
            end
            return du
        end
    timespan = (t[1],t[end]);
    fix=["Ca_mit","Ca_ecs"]
    prob1 = ODEProblem(ode_problem_fun_r,u1,timespan)
    Array(solve(prob1,alg_hints=[:stiff],dtmax=10;saveat=t))
end



function give_rand_p(p_range,p_fixed=nothing,indices=nothing)

    if p_fixed == nothing

        p = [(p_range[j][2] -p_range[j][1])*rand() + p_range[j][1] for j in 1:length(p_range)]

    else

        p =  zeros(length(p_range))

        j = 1

        for i in 1:length(p_range)

            if i in indices

                p[i] = p_fixed[j]

                j += 1

            else

                p[i] = (p_range[i][2] -p_range[i][1])*rand() + p_range[i][1]

            end

        end

    end

    p

end

function give_rand_p_new(p_range,p_fixed=nothing,indices=nothing)

    if p_fixed == nothing

        p = [10^((log10(p_range[j][2]) -log10(p_range[j][1]))*rand() + log10(p_range[j][1])) for j in 1:length(p_range)]
    else

        p =  zeros(length(p_range))

        j = 1

        for i in 1:length(p_range)

            if i in indices

                p[i] = p_fixed[j]

                j += 1

            else

                p[i] = 10^((log10(p_range[i][2]) -log10(p_range[i][1]))*rand() + log10(p_range[i][1]))
 
            end

        end

    end

    p

end




function calc_mean_var(prob,t,p_range,N)
    q =give_rand_p(p_range);
    f = getf(q,prob,t);
    y1 = Array(f)

    y0 = zero(y1)

    v = zero(y1)

    for i in 1:N
        q =give_rand_p(p_range);
        f = getf(q,prob,t);
        y1 = Array(f)

        @. y0 += y1

        @. v += y1^2

    end

    y0 = @. y0/N

    y0_sq = [i.^2 for i in y0]

    v = @. v/N - y0_sq

    y0,v

end

function calc_mean_var_new(q,s,t,p_range,N)
    rd = give_rand_p_new(p_range);
    f = getf_new(rd,q,s,t);
    y1 = Array(f)

    y0 = zero(y1)

    v = zero(y1)

    for i in 1:N
        rd = give_rand_p_new(p_range);
        f = getf_new(rd,q,s,t);
        y1 = Array(f)

        @. y0 += y1

        @. v += y1^2

    end

    y0 = @. y0/N

    y0_sq = [i.^2 for i in y0]

    v = @. v/N - y0_sq

    y0,v

end



function first_order_var(prob,t,p_range,N,y0)

    ys = Array{typeof(y0)}(undef,length(p_range))

    for i in 1:length(p_range)

        y = zero(y0)

        for j in 1:N

            p2 = give_rand_p(p_range)

            p1 = give_rand_p(p_range,[p2[i]],[i])
            
            f1 = getf(p1,prob,t)
            f2 = getf(p2,prob,t)

            yer =  Array(f1) .* Array(f2)

            @. y += yer

        end

        y = @. y/N - y0^2

        ys[i] = copy(y)

    end

    ys

end


function first_order_var_new(q,s,t,p_range,N,y0)

    ys = Array{typeof(y0)}(undef,length(p_range))

    for i in 1:length(p_range)

        y = zero(y0)

        for j in 1:N

            p2 = give_rand_p_new(p_range)

            p1 = give_rand_p_new(p_range,[p2[i]],[i])
            
            f1 = getf_new(p1,q,s,t)
            f2 = getf_new(p2,q,s,t)

            yer =  Array(f1) .* Array(f2)

            @. y += yer

        end

        y = @. y/N - y0^2

        ys[i] = copy(y)

    end

    ys

end

function second_order_var(prob,t,p_range,N,y0)

    ys = Array{typeof(y0)}(undef,Int((length(p_range)*(length(p_range)-1))/2))

    curr = 1

    for i in 1:length(p_range)

        for j in i+1:length(p_range)

            y = zero(y0)

            for k in 1:N

                p2 = give_rand_p(p_range)

                p1 = give_rand_p(p_range,[p2[i],p2[j]],[i,j])
                f1 = getf(p1,prob,t)
                f2 = getf(p2,prob,t)

                y .+=  Array(f1) .* Array(f2)

            end

            y = @. y/N - y0^2

            ys[curr] = copy(y)

            curr += 1

        end

    end

    ys_frst_order = first_order_var(prob,t,p_range,N,y0)

    j = 1

    for i in 1:length(p_range)

        for k in i+1:length(p_range)

            ys[j] = @. ys[j] - ( ys_frst_order[i] + ys_frst_order[k] )

            j += 1

        end

    end

    ys

end


function second_order_var_new(q,s,t,p_range,N,y0)

    ys = Array{typeof(y0)}(undef,Int((length(p_range)*(length(p_range)-1))/2))

    curr = 1

    for i in 1:length(p_range)

        for j in i+1:length(p_range)

            y = zero(y0)

            for k in 1:N

                p2 = give_rand_p_new(p_range)

                p1 = give_rand_p_new(p_range,[p2[i],p2[j]],[i,j])
                f1 = getf_new(p1,q,s,t)
                f2 = getf_new(p2,q,s,t)

                y .+=  Array(f1) .* Array(f2)

            end

            y = @. y/N - y0^2

            ys[curr] = copy(y)

            curr += 1

        end

    end

    ys_frst_order = first_order_var_new(q,s,t,p_range,N,y0)

    j = 1

    for i in 1:length(p_range)

        for k in i+1:length(p_range)

            ys[j] = @. ys[j] - ( ys_frst_order[i] + ys_frst_order[k] )

            j += 1

        end

    end

    ys

end

if length(ARGS)>0 && ARGS[1] == "--run"
    using balancetest
end