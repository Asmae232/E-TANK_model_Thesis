# E-TANK: Carbon-Tax Recycling, Energy Poverty, and Green Innovation

E-TANK is a real, two-agent dynamic general-equilibrium (TANK) model of a
small open economy. It features hand-to-mouth and Ricardian households,
non-homothetic Stone-Geary preferences with a committed energy component,
a nested CES production structure, clean and imported dirty energy, and
endogenous clean-knowledge accumulation through green R&D.

The carbon tax enters as an additive price wedge on imported dirty energy:

    P_d,t = er_t P_oil,t + tau_t.

Carbon-tax revenue is recycled within each period through two policy levers.
The parameter theta determines the share allocated to clean R&D, with the
remainder returned through household transfers. The parameter phi determines
how that transfer envelope is allocated between hand-to-mouth and Ricardian
households.

The quarterly benchmark calibration is disciplined primarily by euro-area
macroeconomic moments, while household heterogeneity is informed by French
and broader European evidence. The model is used to examine how carbon-tax
revenue recycling affects clean innovation, the macroeconomic cost of the
energy transition, household welfare, and the distribution of transition
costs.

The included `e_tank_model.mod` file is the authoritative structural model
source used for my thesis production runs. The reported results are
obtained from deterministic perfect-foresight simulations using a numerical
horizon of 1,600 quarters (400 years).

