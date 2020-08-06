## This code is used to produce predicted trait values in Xu et al. Predictability of leaf traits with climate and elevation: a case study in Gongga Mountain, China (under review, Tree Physiology)

## import site, climate and trait input data
trait <- read.csv('/data/trait.csv', header = T)
site <- read.csv('/data/site.csv', header = T)
climate <- read.csv('/data/climate.csv', header = T)

trait <- merge(trait, site, by='Site.ID', all=T)
all <- merge(trait, climate, by='Site.ID', all=T)
all.dec <- subset(all, Leaf.phenology == 'deciduous')
climate <- merge(climate, site, by='Site.ID', all=T)

##-------- predict X ---------
## 1. parameter
beta <- 146 # Wang et al. (2017)

## 2. define functions
{
  # calculate air pressure in Pa
  calc_patm <- function( elv ){
    #-----------------------------------------------------------------------
    # Input:    - elevation, m (elv)
    # Output:   - float, atmospheric pressure at elevation 'elv', Pa (patm)
    # Features: Returns the atmospheric pressure as a function of elevation
    #           and standard atmosphere (1013.25 hPa)
    # Depends:  - connect_sql
    #           - flux_to_grid
    #           - get_data_point
    #           - get_msvidx
    # Ref:      Allen et al. (1998)
    #-----------------------------------------------------------------------
    
    # Define constants:
    kPo <- 101325   # standard atmosphere, Pa (Allen, 1973)
    kTo <- 298.15   # base temperature, K (Prentice, unpublished)
    kL <- 0.0065    # temperature lapse rate, K/m (Allen, 1973)
    kG <- 9.80665   # gravitational acceleration, m/s^2 (Allen, 1973)
    kR <- 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
    kMa <- 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
    
    # Convert elevation to pressure, Pa:
    patm <- kPo*(1.0 - kL*elv/kTo)**(kG*kMa/(kR*kL))
    
    return (patm)
  }
  
  # calculate K (MM coefficient of Rubisco) in Pa
  calc_k <- function(temp, patm) {
    #-----------------------------------------------------------------------
    # Input:    - float, air temperature, deg C (temp)
    #           - float, atmospheric pressure, Pa (patm)
    # Output:   float, Pa (mmk)
    # Features: Returns the temperature & pressure dependent Michaelis-Menten
    #           coefficient, K (Pa).
    # Ref:      Bernacchi et al. (2001), Improved temperature response 
    #           functions for models of Rubisco-limited photosynthesis, 
    #           Plant, Cell and Environment, 24, 253--259.
    #-----------------------------------------------------------------------
    
    # Define constants
    kc25 <- 39.97      # Pa, assuming 25 deg C & 98.716 kPa
    ko25 <- 2.748e4    # Pa, assuming 25 deg C & 98.716 kPa
    dhac <- 79430      # J/mol
    dhao <- 36380      # J/mol
    kR   <- 8.3145     # J/mol/K
    kco  <- 2.09476e5  # ppm, US Standard Atmosphere
    
    vc <- kc25*exp(dhac*(temp - 25.0)/(298.15*kR*(temp + 273.15)))
    vo <- ko25*exp(dhao*(temp - 25.0)/(298.15*kR*(temp + 273.15)))
    k  <- vc*(1 + kco*(1e-6)*patm/vo)
    
    return(k)
    
  }
  
  # calculate Gstar (CO2 compensation point) in Pa
  calc_gstar_gepisat <- function( temp ) {
    #-----------------------------------------------------------------------
    # Input:    float, air temperature, degrees C (tc)
    # Output:   float, gamma-star, Pa (gs)
    # Features: Returns the temperature-dependent photorespiratory 
    #           compensation point, Gamma star (Pascals), based on constants 
    #           derived from Bernacchi et al. (2001) study.
    # Ref:      Bernacchi et al. (2001), Improved temperature response 
    #           functions for models of Rubisco-limited photosynthesis, 
    #           Plant, Cell and Environment, 24, 253--259.
    #-----------------------------------------------------------------------
    
    # Define constants
    gs25 <- 4.220    # Pa, assuming 25 deg C & 98.716 kPa)
    dha  <- 37830    # J/mol
    kR   <- 8.3145   # J/mol/K
    
    gs <- gs25 * exp( dha * ( temp - 25.0 ) / ( 298.15 * kR * ( temp + 273.15 ) ) )
    
    return( gs )
    
  }
  # conver CO2 from ppm to Pa
  co2_to_ca <- function( co2, patm ){
    #-----------------------------------------------------------------------
    # Input:    - float, annual atm. CO2, ppm (co2), use 400(ppm)
    #           - float, monthly atm. pressure, Pa (patm)
    # Output:   - ca in units of Pa
    # Features: Converts ca (ambient CO2) from ppm to Pa.
    #-----------------------------------------------------------------------
    ca   <- ( 1.e-6 ) * co2 * patm         # Pa, atms. CO2
    return( ca )
  }
}

## 3. read monthly input data 
Tg <- climate$Tg # mean temperature during growing season in degreeC
elv.m <- climate$Elevation # elevation in meter
D0 <- climate$D0 # vapour pressure deficit in kPa
CO2.ppm <- 400 # ppm
TdJ <- climate$TdJ # daytime temperature in July
D.TdJ <- climate$D.TdJ

## 4. predict X
K1 <- calc_k(Tg,calc_patm(elv.m))
Gstar1 <- calc_gstar_gepisat(Tg)
f1 <- exp(-0.0227*(Tg-25)) # the viscosity of water relative to its value at 25˚C
K2 <- calc_k(TdJ,calc_patm(elv.m))
Gstar2 <- calc_gstar_gepisat(TdJ)
f2 <- exp(-0.0227*(TdJ-25)) # the viscosity of water relative to its value at 25˚C
pa <- co2_to_ca(CO2.ppm,calc_patm(elv.m))

E1 <- sqrt(beta*(K1 + Gstar1)/(1.6*f1))
E2 <- sqrt(beta*(K2 + Gstar2)/(1.6*f2)) 
D.sqrt1 <- sqrt(D0*1000)
D.sqrt2 <- sqrt(D.TdJ*1000)

pre.X.Tg <- ((E1+(Gstar1*D.sqrt1/pa))/(E1+D.sqrt1)) # under Tg
pre.X.TdJ <- ((E2+(Gstar2*D.sqrt2/pa))/(E2+D.sqrt2)) # under TdJ


#
##-------- predict Vcmax25 ---------
## 1. parameter
phi0 <- 0.085 # mol C/ mol photon
c <- 0.41
R <- 8.3145 # J mol-1 K-1
Ha.v <- 72000

## 2. define functions
# Arrhenius function modified terms
Arrhenius.t <- function(dS,Tl,Tref){
  Hd <- 200000 # J mol-1  
  R <- 8.3145 # J mol-1 K-1
  #        x<-A/C
  x <-(1+exp((Tl*dS-Hd)/(Tl*R)))/(1+exp((Tref*dS-Hd)/(Tref*R)))
  return(x)
}

## 3. read input data 
Iabs <- all.dec$R0
X <- all.dec$X
elv <- all.dec$Elevation
Tg2 <- all.dec$Tg
TdJ2 <- all.dec$TdJ

Tref.Vcmax<- 25 + 273.15
Tl.Tg <- Tg2 + 273.15 
Tl.TdJ <- TdJ2 + 273.15 
dSv.Tg <- 668.39 - 1.07*Tl.Tg
dSv.TdJ <- 668.39 - 1.07*Tl.TdJ

## 4. estimate Vcmax25
ca <- co2_to_ca(CO2.ppm,calc_patm(elv))
ci <- X*ca
K3 <- calc_k(Tg2,calc_patm(elv))
Gstar3 <- calc_gstar_gepisat(Tg2)
K4 <- calc_k(TdJ2,calc_patm(elv))
Gstar4 <- calc_gstar_gepisat(TdJ2)
m1 <- (ci - Gstar3)/(ci + 2*Gstar3)
m2 <- (ci - Gstar4)/(ci + 2*Gstar4)
pre.Vc.Tg <- phi0*Iabs*(ci + K3)/(ci + 2*Gstar3)*sqrt(1-(c/m1)^(2/3))
pre.Vc.TdJ <- phi0*Iabs*(ci + K4)/(ci + 2*Gstar4)*sqrt(1-(c/m2)^(2/3))

# convert to standard temperature
pre.Vc25.Tg <- pre.Vc.Tg/exp(Ha.v*(Tl.Tg - Tref.Vcmax)/(Tref.Vcmax*R*Tl.Tg))*Arrhenius.t(dSv.Tg,Tl.Tg,Tref.Vcmax)
pre.Vc25.Tg7 <- pre.Vc.TdJ/exp(Ha.v*(Tl.TdJ - Tref.Vcmax)/(Tref.Vcmax*R*Tl.TdJ))*Arrhenius.t(dSv.TdJ,Tl.TdJ,Tref.Vcmax)


#
##-------- predict Ma ---------
## 1. read input data
RLAI <- climate$RLAI # leaf-area-index weighted par
f <- climate$f # the ratio of growing season length to the number of days in the year
ap <- climate$alpha_p # AET/PET

## 2. predict Ma
pre.Ma.Tg <- exp(1.22*log(RLAI) + 0.78*log(f) - 0.06*Tg - 0.6*log(ap) + 1.7)
pre.Ma.TdJ <- exp(1.22*log(RLAI) + 0.78*log(f) - 0.06*TdJ - 0.6*log(ap) + 1.7)


#
##-------- predict Narea ---------
## 1. read input data
Vc25 <- all.dec$Vcmax25
Ma <- all.dec$Ma

## 2. predict Narea
# using Ma and Vcmax25 directly
pre.Narea1 <- 0.0215*Ma + 0.0026*Vc25

# using two-step in Dong et al. (2017)
Ns <- 10^(-2.67)*Ma^0.99
Nr <- 0.003135*Vc25
pre.Narea2 <- Ns + 7.23*Nr