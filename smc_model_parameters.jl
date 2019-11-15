using Parameters

@with_kw struct Control_params @deftype Float64
    SERCA = 17.548
    SERCA_s = 5
    Qvocc= 1.947023958788220*0.01
    # Leaks
    g_leak_mit= 0.0020741774577606935
    g_leak_ecs= 0.002234985778040491
    k_leak_er= 0
    # steady state values
    Ca_cyt_infty= 0.1 # muM
    Ca_er_infty= 500 # muM
    ATP_infty= 3000
    ADP_infty= 10
    Qryr= 2  # Junming
    Qip3r = 4 # /s
    MyoTot= 10 # muM check this! @TODO
    Qncx= 1.85 # computed using Johny
    Qpmca= 0.103211*0.1  # Johny
    Vnclx= 0.35*0.0733 #0.026673080107946  # Wacquier 2016
    Vmcu= 0.00006*0.0733  # Wacquier 2016
    kncx2 = 1.248e1 *0.1# Demir muA/muM
    L = 0.1
end;

@with_kw struct State @deftype Float64
    Ca_cyt=0.1  # muM
    Ca_er=500  # muM
    Ca_ecs=1300  # muM
    Ca_mit=.1  # muM
    Ca_mit_source=0
    Ca_ecs_source=0
    Na_cyt= 1000
    Na_ecs= 138e3
    K_cyt= 140e3
    K_ecs= 3.5e3
    # mitochondrial membrane potential
    Phi_mit= 160  # mV
    Phi_ecs= -54.03 # mV # plasma membrane potential
    RyRR10= 0.003278103178479715
    IP3RX00= 0.4843494729034235
    IP3RX01= 0.22898031672932792
    IP3RX10= 0.1482017534536128
    RyRR11= 5.853461104503355e-07
    MyoMp= 0.018698802912942764
    MyoAMp= 0.014183208692463749
    IP3= 0.46938692642554697
    RyRR01= 0.0001778462416021997
    MyoAM= 0.06061200295918522
    # Energetics
    ATP= 3000
    ADP= 300
    Glu= 1
    pH_cyt= 7.1
    pH_mit= 7.8
    G = .00016682728821818067e-1
    PIP2 = 4e7/(6.022e23*0.7e-12)*1e6 # muM 
  #  Rs = 2e4/(6.022e23*0.7e-12)*1e6 *0.85 #muM
   # Rsp = 0.0 #muM  
end;

@with_kw struct Fixed_params @deftype Float64
    # RyR
    Kr1 = 2.5   # 2.5 (muM)^-2 s^-1
    Kr2= 1.5  # /muM/s
    Kmr1= 7.6  # /s
    Kmr2= 84  # s
    Vol_cyt = 0.7e-12 # L
    V_delta= 0.05
    K_PLC= 0.1
    k_delta= 1.5
    r_5p= 0.04
    v_3k=  2
    IP3_3K_K_D= 0.7
    IP3_3K_K_3= 1
    Vs= 30
    nVs= 20
    p2= 0.016  # Wacquier 2016
    kmx= 0.008  # 1/S Wacquier 2016
    p3= 0.05  # mv^-1  Wacquier 2016
    RV1= 2000
    RVer= 10  # Vcyt/Ver
    RVmit= 1/0.0733  # Vcyt/V_mit wacqiuer
    F= 96480 #C/ mol
    R= 8315
    Temp= 310.16
    FRT= F/R/Temp
    K1= 6  # Wacquier 2016
    K2= 0.38  # Wacquier 2016
    L= 50  # Wacquier 2016
    p1= 0.1  # Wacquier 2016
    # parameters for buffer
    K_d= 260E-3  # Johny et al. 2015
    K_dB= 530E-3  # Johny et al. 2015
    S_CM= 30  # Johny et al. 2015
    B_F= 30  # Johny et al. 2015
    K_CalrC= 2E3
    CalrC= 36E3/5
    K_CalrP= 10
    CalrP= 3.6E3/5
    # PMCA channel
    kpmca= 0.15  # muM Johny et al. 2017
    kpmca_ATP= 1.05e3 # http=//jgp.rupress.org/content/138/4/381
    knaca= 5.82E5  # 2.91E5  # \muM/s Lyer.  2005 Biophysical 87
    Nao= 138E3  # muM  Lyer.  2005 Biophysical 87
    Nai= 10E3  # mu M  Lyer.  2005 Biophysical 87
    K_cyt_infty= 140E3 # mu M
    K_ecs_infty= 5000
    Cao= 2E3  # muM  Lyer.  2005 Biophysical 87
    kmna= 87.5E3  # muM  Lyer.  2005 Biophysical 87
    kmca= 1.38E3  # muM  Lyer.  2005 Biophysical 87
    ksat= 0.2
    eta= 0.35
    Vcm= -54.03
    # IP3R 4 STATE DYNAMICS PARAMETERS
    a1 = 167.6 # 1/muM/s
    a2 = 3.81#1/muM/s
    a3 = 413.4 #1/muM/s
    a4 = 0.3101 #1/muM/s
    a5 = 53.9 #1/muM/s
    b1 = 228 #1/s
    b2 = 0.409#1/s
    b3 = 188.5#1/s
    b4 = 0.096#1/s
    b5 = 4.52#1/s 
    #################################
    ######### Myosin
    MyoK1= 17 # (muM)^-3 s^-1
    MyoK2= 0.5 # s^-1
    MyoK3= 0.4 # s^-1
    MyoK4= 0.1 # s^-1
    MyoK5= 0.5 # s^-1
    MyoK6= 17 # (muM)^-3 s^-1
    MyoK7= 0.1 # s^-1
    C= 20e-12 # F from Johny
    # Johny BK channel
    k_bk= 1.85*11000
    # NCX from Williams
    Kncxca= 1380 # muM
    Kncxna= 8.75e4 #muM
    kncxsat= .1 # dimensionless
    etancx= 0.35 # dimensionless
    A_m = 1.5340e-4
    # NCX from Demir
    kncx1 = 0.125 # muM
    dncx = 1e8
    # IP3 model
    kdegDAG= 1 # 1/s
   
    IP3_eta= 2.781e-5*(6.023e23*.7e-12)*1e-6 # 1/(muM s)# Lemon
    IP3_kc= .400 # muM# Lemon
    IP3_kdeg= 1.25 #1/s  # Lemon
    G_zeta = 0.85  # Lemon
    G_RT = 2e4/(6.022e23*0.7e-12)*1e6 #muM # Lemon
    G_kd= 0.15 # 1/s# Lemon
    G_ka= 0.017 # 1/s# Lemon
    G_delta= 1.235e-3# Lemon
    
    G_T= 0.23722541158608913 # muM# Lemon
    PIP2T = 5e7/(6.022e23*0.7e-12)*1e6 # muM # Johny
   
    Rs_kr = 1.75e-4 # 1/s # Lemon
    Rs_kp = 0.03 #1/s# Lemon
    Rs_ke = 6.0e-3 #1/s# Lemon
    Rs_k1 = 60.0 #mu M # Johny
    Rs_k2 = 100.0 # muM # Lemon
    PIP2_rr = 0.1 # 1/s # Lemon

end;
