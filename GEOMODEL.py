def adjust_model(p='', prm_phi='', srm_phi='', srm_T='', sq_db='', sq_angle='', sq_freq=''):
    
    # Function that creates a base model of GEO600 with set values and updates that model for all inputed parameter values
    # Returns base model if no parameter values are inputed
    
    # Import packages
    import finesse
    finesse.configure(plotting=True)
    import numpy as np
    
    # Constant Variables in GEO600 Model 
    laserWav = 1064e-9
    mRes = 5.6
    mirrorLoss = 130e-6
    fRes = 0.5
    QRes = 10**8
    
    # Implement and Create GEO600 Model
    model = finesse.Model()
    model.parse(
        f"""
        
        # Main Laser:
        
        l laser P=50 f=0 phase=0
        s s1 laser.p1 prm.p1 L=0
        
        ##########################################
        
        # Power-Recycling Mirror:
        
        m prm T=900e-6 L={mirrorLoss} phi=90         # Power-recycling mirror
        s prcav prm.p2 BS.p1 L=1.1463                # Power-recycling cavity
        
        ##########################################
        
        # Beam splitter:
        
        bs BS T=0.513872 L={mirrorLoss} alpha=42.834 
        
        ##########################################
        
        # North Arm:
        
        s sy BS.p2 bsy.p1 L=598.5682
        bs bsy T=8.3e-6 L={mirrorLoss}               # Far North mirror
        s sy2 bsy.p2 my.p1 L=597.0241
        m my T=13e-6 L={mirrorLoss} phi=90.017       # Central North mirror
        
        # Central North mirror set 50 picometers off dark fringe for DC readout
        
        ##########################################
        
        # East Arm:
        
        s sx BS.p3 bsx.p1 L=598.4497
        bs bsx T=8.3e-6 L={mirrorLoss}               # Far East mirror
        s sx2 bsx.p2 mx.p1 L=597.0630
        m mx T=13e-6 L={mirrorLoss}                  # Central East mirror
        
        ##########################################
        
        # Signal-Recycling Mirror:
        
        s srcav BS.p4 srm.p1 L=1.109                 # Signal-recycling cavity 
        m srm T=0.09995 L=50e-6 phi=90                # Signal-recycling mirror
        
        # Adjust the phi parameter (between 0 and 90) of the srm mirror to detune GEO600
        
        ##########################################
        
        # Squeezer:
        
        s darkP srm.p2 FI.p1
        dbs FI #Faraday isolator
        s lsqz sq1.p1 FI.p2
        sq sq1 db=6 angle=103.6                      # Squeezer
        
        # A squeezed source is injected into a Faraday Isolater (directional beam splitter in Finesse)
        
        ##########################################
        
        # Mirror Suspension:
        
        pendulum py my.mech mass={mRes} fz={fRes} Qz={QRes} fyaw=0.4 fpitch=2.1
        pendulum px mx.mech mass={mRes} fz={fRes} Qz={QRes} fyaw=0.4 fpitch=2.1
        pendulum pbsy bsy.mech mass={mRes} fz={fRes} Qz={QRes} fyaw=0.4 fpitch=2.1
        pendulum pbsx bsx.mech mass={mRes} fz={fRes} Qz={QRes} fyaw=0.4 fpitch=2.1
        pendulum ppr prm.mech mass=2.92 fz=0.57 Qz={QRes} fyaw=0.76 fpitch=1.5
        pendulum psr srm.mech mass=2.92 fz=0.57 Qz={QRes} fyaw=0.76 fpitch=1.5
        pendulum pbs BS.mech mass=9.3 fz=0.4 Qz={QRes} fyaw=0.7 fpitch=1.35
        
        ##########################################
        
        # Set starting signal frequency:
        
        fsig(1)
        
        """)
    
    # Update the model for every inputed parameter value
    if p!='':
        model.laser.P=p
    if prm_phi!='':
        model.prm.phi=prm_phi
    if srm_phi!='':
        model.srm.phi=srm_phi
    if srm_T!='':
        model.srm.T=srm_T
    if sq_db!='':
        model.sq1.db=sq_db
    if sq_angle!='':
        model.sq1.angle=sq_angle
    if sq_freq!='':
        model.sq1.f=sq_freq
    
    # Returns the Finesse Model
    return model