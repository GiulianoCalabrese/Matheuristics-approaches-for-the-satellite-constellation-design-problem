"""
Structures de donnée modélisant les objects du problème. 
Contient aussi les fonctions permettant de construire ces structures. 
"""

"""
Calculer le rayon d'orbite d'une trajectoire à partir de l'indice du rayon d'orbite `n` et 
du nombre de période terreste de simulation `m_simu`.
"""
calculRayonOrbite(n::Real, m_simu::Real) = (μ*(m_simu*T_terre/(2π*n))^2)^(1/3)

"""
Calculer le nombre minimal d'instants à considérer par intervalles de surveillance telque 
l'écart angulaire correspondant à l'instersection des deux disques de vision d'un satellites 
pris à deux instants successifs soit au moins égale à `p*θ`, où θ est le rayon angulaire du 
du disque terrestre vu par un satellite de demi-angle d'ouverture `angleOuverture` sur 
l'orbite la plus basse. 
"""
function calculNbIntants(m_simu::Int, n_0::Int, angleOuverture::Real, n_max::Real, p::Real)
   a = calculRayonOrbite(n_max, m_simu)
   θ = calcul_θ(a, angleOuverture)
   return floor(Int, 2π*n_max/((2-p)*θ*n_0))
end

"""
Une instance du type composé `DonnesSCDP` représente une instance du SCDP, et contient 
donc les données de cette instance. Une instance contient aussi les paramètres de simulation
de la dynamique d'une constellation quelconque solution de ce problème.

# Champs 
- `coordsCibles`: Coordonnées des cibles à surveiller. De la forme [[latidude_1,
   longitude_1], [latidude_2, longitude_2], ...] où les grandeurs sont exprimées en radians
   et latitude est la latidude géodésique.  
- `angleOuverture` : angle d'ouverture des satellites de la constellation en radians. 
- `instantInitial` : date à laquelle les satellites de la constellation sont mis en orbite. 
   C'est donc la date à laquelle commence la simulation.
- `m_simu` : La simulation est calculée sur une durée de `m_simu*T_terre`. 
- `n_0` : La durée de simulation est découpée en intervalles de temps de durée égale. Chaque
   cible doit être survolée au moins une fois sur chacun de ces intervalles. n_0 est le 
   nombre de ces intervalles (sur la durée de simulation).  

- `nbCibles` : Nombre de cibles à surveiller. 
- `nbMaxSat` : Il existe une solution cette instance avec `nbMaxSat`. 
- `n_min`, `n_max` : Le rayon d'orbite d'un satellite est définit par un entier n telque 
   n*T = m_simu*T_terre où T est la période de rotation du satellite. n est le nombre de
   période du satellite sur la durée de la simulation. La troisième loi de Kepler permet de
   déduire le rayon d'orbite. Enfin, la constrainte LEO donne des bornes pour le rayon
   d'orbite, dont on déduit des bornes pour n qui sont `n_min` et `n_max`.
- `nbInstants` : Nombres d'instants par intervalles de surveillance qui seront considérés
   lors de la simulation. 
- `dt` : Pas de temps en secondes.
- `ts` : vecteur de temps sur lequel l'analyse est conduite, en seconde écoulées depuis
   l'instant initial de prop.
- `ìnstantInitialJulian` : instant initiale en jours juliens. 
- `mats_ECI_ECEF` : Vecteur de matrice de rotation. La multiplication par mats_ECI_ECEF[i] 
   permet de passer des coordonnées ECI aux coordonnées ECEF à l'instant ts[i]. 
- `coordsCiblesECEF` : Coordonnées des cibles en mètres dans le référentiel ECEF PEF. 
- `nbCS` : Ce paramètre concerne l'inclinaison, la longitude du noeud ascendant et 
   l'anomalie moyenne des satellites appartenant aux constellations solution de 
   l'instance. C'est le nombre de chiffre significatif après la virgule que l'on doit avoir 
   sur ces variables pour être cohérent avec la discrétisation temporelle. Ce paramètre 
   intervient lors du codage en binaire des variables. 
   """
mutable struct DonneesSCDP
   # Les champs suivants sont les données d'une instance du problème.
   coordsCibles::AbstractVector{<:AbstractVector{<:Real}}
   angleOuverture::Real
   instantInitial::DateTime
   m_simu::Int 
   n_0::Integer

   # Les champs suivants se calculent à partir des champs ci-dessus. 
   nbCibles::Int
   nbMaxSat::Int 
   n_min::Int
   n_max::Int

   nbInstants::Int
   dt::Real
   ts::AbstractVector{<:Real}

   instantInitialJulian::Float64
   mats_ECI_ECEF::AbstractVector{<:AbstractMatrix{<:Real}}
   coordsCiblesECEF::AbstractVector{<:AbstractVector{<:Real}}

   nbCS::Int

   """
   # Arguments
   - `p` : le pas de temps pour la simulation est définit telque deux disques à la surface 
   de la Terre, correspondant aux zones vu par un satellite à deux instants successif, se 
   superposent de 100p %. 
   """
   function DonneesSCDP(coordsCibles, angleOuverture, nbMaxSat, instantInitial, m_simu, n_0, dt, p = 0.2)

	  nbCibles = length(coordsCibles)
	  # nbMaxSat = nbCibles*n_0 
	  n_min = ceil(Int, (m_simu*T_terre*√(μ/a_max^3))/(2π))
	  n_max = floor(Int, (m_simu*T_terre*√(μ/a_min^3))/(2π))
	  
	  nbInstants = Int(T_terre/(n_0*dt)) #calculNbIntants(m_simu, n_0, angleOuverture, n_max, p)
	  dt = (T_terre*m_simu)/(nbInstants*n_0)
	  ts = LinRange(0, m_simu*T_terre, Int(nbInstants*n_0 + 1))

	  instantInitialJulian = datetime2julian(instantInitial) 
	  mats_ECI_ECEF = 
		  @. r_eci_to_ecef(TOD(), PEF(), instantInitialJulian + ts/T_terre)
	  coordsCiblesECEF = @. geodetic_to_ecef(getindex(coordsCibles, 1),
										  getindex(coordsCibles, 2), 0)

	  ϵ = 2n_max*dt*π/(m_simu*T_terre)
	  nbCS = ceil(Int, -log10(ϵ))

	  return new(coordsCibles, angleOuverture, instantInitial, m_simu, n_0, 
				 nbCibles, 
				 nbMaxSat, 
				 n_min, 
				 n_max, 
				 nbInstants,
				 dt,
				 ts, 
				 instantInitialJulian, 
				 mats_ECI_ECEF, 
				 coordsCiblesECEF,
				 nbCS)
   end
end

"""
Générer une instance aléatoire du SCDP. Renvoie une variable de type `DonneesProb`.
"""
function instanceAleatoire(nbCibles::Int, angleOuverture::Real, nbMaxSat::Int, instantInitial::DateTime,  
      m_simu::Int, n_0::Int, p::Real = 0.2)
   # coordsCibles = [SA[π*rand() - π/2, 2π*rand() - π] for _ in 1:nbCibles]
   # coordsCibles = [SA[π*rand() - π/2, 2π*rand()] for _ in 1:nbCibles]
   phi  = 2. * pi * rand()
   csth = 1.0 - 2.0 * rand()
   snth = sqrt(1.0 - csth*csth)
   coordsCibles = [SA[asin(1.0 - 2.0 * rand()), 2. * pi * rand()] for _ in 1:nbCibles]
   return DonneesSCDP(coordsCibles, angleOuverture, nbMaxSat, instantInitial, m_simu, n_0, p)
end

"""
Les données définissant une constellation de satellites sont les éléments orbitaux des 
satellites qui composent cette constellation (rayon d'orbite, inclinaison, longitude du 
noeud ascendant, anomalie moyenne initiale) ainsi que la date à laquelle ces éléments 
sont considérés (on suppose que les éléments orbitaux de chaque satellite de la
constellation sont considérés à la même date). 

La structure suivante permet de représenter une constellation, indépendamment d'une instance
du SCDP. Les fonctions prenant en paramètre une constellation reçoivent générallement
une constellation représentée par ce type. 
"""
struct Conste 
   date::Float64 # Date des éléments orbitaux en jour julien.
   nbSats::Int # Nombre de satellites qui composent la constellation
   as::Vector{Float64} # Rayons d'orbite
   is::Vector{Float64} # Inclinaisons
   Ωs::Vector{Float64} # Longitudes des noeuds ascendants
   Ms::Vector{Float64} # Anomalies à l'instant date
end

"""
On définit d'autres représentations d'une constellation qui seront utilisées par NSGAII. Ces 
représentations sont associées à une instance du SCDP. 

   # Représentation de taille fixée avec des éléments orbitaux groupés par satellites. 
   
Une constellation est représentée par un array (de flottant) de la forme 
   [b₁, n₁, i₁, Ω₁, M₁, ..., bᵣ, nᵣ, iᵣ, Ωᵣ, Mᵣ] 
où 
- bᵢ dans {0, 1} est un booléen indiquant si le ième satellite est actif. 
- nᵢ dans {n_min, ..., n_max} est l'indice du rayon d'orbite du ième satellite. 
- iᵢ dans [0, π] est l'inclinaison du ième satellite. 
- Ωᵢ dans [-π, π] est la longitude du noeud ascendant du ième satellite. 
- Mᵢ dans [0, 2π] est l'anomalie moyenne initiale du ième satellite. 
- r est le nombre maximale de satellites d'une constellation. Ce nombre est définit par 
l'instance du SCDP à laquelle est associée cette constellation. 

   # Représentation de taille fixée avec des éléments orbitaux groupés par type. 
   
Une constellation est représentée par un array (de flottant) de la forme 
   [b₁, ..., bᵣ, n₁, ..., nᵣ, i₁, ..., iᵣ, Ω₁, ..., Ωᵣ, M₁, ..., Mᵣ] 
où 
- bᵢ dans {0, 1} est un booléen indiquant si le ième satellite est actif. 
- nᵢ dans {n_min, ..., n_max} est l'indice du rayon d'orbite du ième satellite. 
- iᵢ dans [0, π] est l'inclinaison du ième satellite. 
- Ωᵢ dans [-π, π] est la longitude du noeud ascendant du ième satellite. 
- Mᵢ dans [0, 2π] est l'anomalie moyenne initiale du ième satellite. 
- r est le nombre maximale de satellites d'une constellation. Ce nombre est définit par 
l'instance du problème à laquelle est associée cette constellation. 

   # Représentation de taille variable

Une constellation est représentée par un array (de flottant) de la forme 
   [n₁, i₁, Ω₁, M₁, ..., nᵣ, iᵣ, Ωᵣ, Mᵣ] 
où 
- nᵢ dans {n_min, ..., n_max} est l'indice du rayon d'orbite du ième satellite. 
- iᵢ dans [0, π] est l'inclinaison du ième satellite. 
- Ωᵢ dans [-π, π] est la longitude du noeud ascendant du ième satellite. 
- Mᵢ dans [0, 2π] est l'anomalie moyenne initiale du ième satellite. 
- r est le nombre de satellites de la constellation. 
"""

"""
Représente le groupement par satellite.
"""
struct GrpSat end 
"""
Représente le groupement par type. 
"""
struct GrpType end 

"""
Calculer le nombre de satellites actifs de la représentation d'une constellation de taille 
fixe groupé par satellite `conste`, solution de l'instance `inst`. 
"""
function calculNbSatsActifsConsteFixe(::Type{GrpSat}, conste::AbstractVector{<:Real}, 
      inst::DonneesSCDP)
   actifs = conste[1:5:5inst.nbMaxSat]
   return sum(actifs)
end

"""
Calculer le nombre de satellites actifs de la représentation d'une constellation de taille 
fixe groupé par type `conste`, solution de l'instance `inst`. 
"""
function calculNbSatsActifsConsteFixe(::Type{GrpType}, conste::AbstractVector{<:Real}, 
      inst::DonneesSCDP)
   actifs = conste[1:inst.nbMaxSat]
   return sum(actifs)
end

"""
Extraire les éléments orbitaux des satellites (actifs) d'une constellation `conste`,
solution de l'instance `inst`, représentée en taille fixe avec des éléments orbitaux
groupés par satellite.
"""
function extraireDeConsteFixe(::Type{GrpSat}, conste::AbstractVector{<:Real}, 
      inst::DonneesSCDP)
   actifs = Bool.(conste[1:5:5inst.nbMaxSat])
   nbSat = sum(actifs)
   ns = conste[2:5:5inst.nbMaxSat][actifs]
   is = conste[3:5:5inst.nbMaxSat][actifs]
   Ωs = conste[4:5:5inst.nbMaxSat][actifs]
   Ms = conste[5:5:5inst.nbMaxSat][actifs]
   return nbSat, ns, is, Ωs, Ms
end

"""
Extraire les éléments orbitaux des satellites (actifs) d'une constellation `conste`,
solution de l'instance `inst`, représentée en taille fixe avec des éléments orbitaux
groupés par type.
"""
function extraireDeConsteFixe(::Type{GrpType}, conste::AbstractVector{<:Real}, 
      inst::DonneesSCDP)
   actifs = Bool.(conste[1:inst.nbMaxSat])
   nbSat = sum(actifs)
   ns = conste[(inst.nbMaxSat+1):(2inst.nbMaxSat)][actifs]
   is = conste[(2inst.nbMaxSat+1):(3inst.nbMaxSat)][actifs]
   Ωs = conste[(3inst.nbMaxSat+1):(4inst.nbMaxSat)][actifs]
   Ms = conste[(4inst.nbMaxSat+1):(5inst.nbMaxSat)][actifs]
   return nbSat, ns, is, Ωs, Ms
end

"""
Extraire les éléments orbitaux des satellites d'une constellation `conste` représentée en
taille variable. 
"""
function extraireDeConsteVar(conste::AbstractVector{<:Real})
   ns = conste[1:4:end]
   is = conste[2:4:end]
   Ωs = conste[3:4:end]
   Ms = conste[4:4:end]
   return length(ns), ns, is, Ωs, Ms 
end

"""
Etant donnée une représentation de taille fixe groupée selon `Grp` d'une constellation
`conste`, solution de l'instance du problème `inst`, renvoie cette constellation convertit
en type 'Conste`.
"""
function consteFixe2Conste(Grp::Union{Type{GrpSat}, Type{GrpType}}, 
      conste::AbstractVector{<:Real}, inst::DonneesSCDP)
   nbSats, ns, is, Ωs, Ms = extraireDeConsteFixe(Grp, conste, inst)
   as = Float64[calculRayonOrbite(n, inst.m_simu) for n in ns]
   return Conste(inst.instantInitialJulian, nbSats, as, is, Ωs, Ms)
end

"""
Etant donnée une représentation de taille variable d'une constellation `conste`,  solution
de l'instance du problème `inst`, renvoie cette constellation convertie en type 'Conste`.
"""
function consteVar2Conste(conste::AbstractVector{<:Real}, inst::DonneesSCDP)
   nbSats, ns, is, Ωs, Ms = extraireDeConsteVar(conste)
   as = Float64[calculRayonOrbite(n, inst.m_simu) for n in ns]
   return Conste(inst.instantInitialJulian, nbSats , as, is, Ωs, Ms)
end
